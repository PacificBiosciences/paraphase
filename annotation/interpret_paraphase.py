import urllib
import logging
import requests
import json
import yaml
import argparse
import os
import time
from argparse import RawTextHelpFormatter
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from typing import Dict, List, Tuple, Optional, Any
from dataclasses import dataclass

# ============================================================================
# CONSTANTS
# ============================================================================

# Default transcript file path (relative to script directory)
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DEFAULT_TRANSCRIPT_FILE = os.path.join(
    _SCRIPT_DIR, "data", "paraphase_genes_transcript_start_end.txt"
)

# Default config URL
DEFAULT_CONFIG_URL = (
    "https://raw.githubusercontent.com/PacificBiosciences/paraphase/main/"
    "paraphase/data/38/config.yaml"
)

# Variants to always include in ClinVar results
CLINVAR_WHITELIST = [
    # three variants in CYP21A2
    "6:g.32039081C>G",
    "6:g.32039810T>A",
    "6:g.32039816T>A",
    # c.840C>T in SMN1
    "5:g.70951946C>T",
]

# Medical regions for annotation
MEDICAL_REGIONS = [
    "rccx",
    "smn1",
    "pms2",
    "strc",
    "cfc1",
    "ikbkg",
    "ncf1",
    "neb",
    "f8",
    "hba",
    "TNXB",
    "OTOA",
    "GBA",
    "CFHclust",
]

# Regions that don't need VCF files
REGIONS_NO_VCF = ["f8", "CFHclust"]

# Genes that always need to be reported
ALWAYS_REPORT_GENES = ["SMN1", "SMN2", "PMS2", "STRC", "IKBKG", "NCF1"]

# CSV Headers
SUMMARY_HEADER = [
    "Region",
    "Region_median_depth",
    "Gene",
    "Total_haplotypes",
    "Num_haplotypes_with_pathogenic_variants",
    "Num_haplotypes_without_pathogenic_variants",
    "Useful_information_from_Paraphase_JSON",
]

DETAILED_HEADER = [
    "Region",
    "Gene",
    "Transcript",
    "Transcript_start",
    "Transcript_end",
    "Haplotype",
    "Haplotype_boundary",
    "Boundary_check",
    "Pathogenic_variants_annotated",
]

# ============================================================================
# API SESSION MANAGEMENT
# ============================================================================


class APISessionManager:
    """Manages HTTP session with rate limiting and retry logic."""

    def __init__(self, min_request_interval: float = 0.34):
        """
        Initialize API session manager.

        Args:
            min_request_interval: Minimum seconds between API requests
        """
        self.min_request_interval = min_request_interval
        self.last_request_time = 0
        self.session = self._create_session()

    def _create_session(self) -> requests.Session:
        """
        Create a requests session with retry logic for handling rate limits.

        Returns:
            requests.Session with retry adapter
        """
        session = requests.Session()

        # Define retry strategy
        retry_strategy = Retry(
            total=5,  # Total number of retries
            backoff_factor=2,  # Wait 2, 4, 8, 16, 32 seconds between retries
            status_forcelist=[429, 500, 502, 503, 504],  # Retry on these status codes
            allowed_methods=["GET", "POST"],  # Only retry on these methods
            respect_retry_after_header=True,  # Respect Retry-After header
        )

        adapter = HTTPAdapter(max_retries=retry_strategy)
        session.mount("http://", adapter)
        session.mount("https://", adapter)

        return session

    def request(self, method: str, url: str, **kwargs) -> requests.Response:
        """
        Make a rate-limited API request with automatic retry on 429 errors.

        Args:
            method: HTTP method (get, post, etc.)
            url: Request URL
            **kwargs: Additional arguments to pass to requests

        Returns:
            Response object
        """
        # Rate limiting: ensure minimum time between requests
        current_time = time.time()
        time_since_last_request = current_time - self.last_request_time
        if time_since_last_request < self.min_request_interval:
            sleep_time = self.min_request_interval - time_since_last_request
            logging.debug(f"Rate limiting: sleeping for {sleep_time:.2f} seconds")
            time.sleep(sleep_time)

        self.last_request_time = time.time()

        # Make the request with retry logic
        max_retries = 3
        for attempt in range(max_retries):
            try:
                response = self.session.request(method, url, **kwargs)

                # If we get a 429, wait and retry
                if response.status_code == 429:
                    retry_after = response.headers.get("Retry-After")
                    if retry_after:
                        try:
                            wait_time = int(retry_after)
                        except (ValueError, TypeError):
                            wait_time = 2**attempt
                    else:
                        # Exponential backoff: 2^attempt seconds
                        wait_time = 2**attempt

                    if attempt < max_retries - 1:
                        logging.warning(
                            f"Rate limit exceeded (429). Waiting {wait_time} seconds before retry "
                            f"(attempt {attempt + 1}/{max_retries})"
                        )
                        time.sleep(wait_time)
                        continue
                    else:
                        response.raise_for_status()

                # Success - return the response
                response.raise_for_status()
                return response

            except requests.exceptions.RequestException as e:
                if attempt < max_retries - 1:
                    wait_time = 2**attempt
                    logging.warning(
                        f"Request failed: {e}. Retrying in {wait_time} seconds "
                        f"(attempt {attempt + 1}/{max_retries})"
                    )
                    time.sleep(wait_time)
                else:
                    # Last attempt failed, re-raise the exception
                    raise


# ============================================================================
# DATA STRUCTURES
# ============================================================================


@dataclass
class RegionInfo:
    """Container for region-specific information from JSON."""

    # useful fields from paraphase json
    report_info: Optional[str]
    # region depth from paraphase json
    region_depth: Optional[float]
    # haplotype boundaries from paraphase json
    haplotype_boundary: Dict[str, str]
    # haplotypes with known variants that make them nonfunctional
    known_nonfunctional_haplotypes: List[str]


# ============================================================================
# CONFIGURATION AND SETUP
# ============================================================================


def load_parameters():
    parser = argparse.ArgumentParser(
        description="Proof-of-concept script to annotate Paraphase results",
        formatter_class=RawTextHelpFormatter,
    )
    inputp = parser.add_argument_group("Input Options")
    outputp = parser.add_argument_group("Output Options")
    inputp.add_argument(
        "-i",
        "--input",
        help="Path to Paraphase JSON file.\n"
        + f"Please make sure that the Paraphase VCF files are in the $sampleID_paraphase_vcfs directory in the same folder.",
        required=True,
    )
    outputp.add_argument(
        "-o",
        "--out",
        help="Output directory",
        required=True,
    )
    parser.add_argument(
        "-r",
        "--region",
        help="Specify the region(s) you want to annotate (separated by comma).\n"
        + "If not provided, the default is to annotate the predefined set of medical regions: rccx, smn1, pms2, strc, cfc1, ikbkg, ncf1, neb, f8, hba, TNXB, OTOA, GBA, CFHclust.\n"
        + "Use 'all' to annotate all regions in the config file (this could take a few hours).",
        required=False,
    )
    parser.add_argument(
        "-c",
        "--config",
        help="Paraphase config file. Using default paraphase config if not provided.",
        required=False,
    )
    parser.add_argument(
        "-l",
        "--log-level",
        help="Set the logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL). Default: INFO",
        required=False,
        type=str.upper,
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
    )
    parser.add_argument(
        "--api-delay",
        help="Minimum seconds between API requests (rate limiting). Default: 0.34 (~3 requests/sec). "
        "Increase this value if you encounter 'Too Many Requests' errors.",
        required=False,
        type=float,
        default=0.34,
    )
    return parser


def search_clinvar_by_variant(
    variant: str, api_session: APISessionManager
) -> List[str]:
    """
    Search ClinVar for a specific variant.

    Args:
        variant: Variant in HGVS format
        api_session: API session manager

    Returns:
        List of ClinVar IDs
    """
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
    esearch_url = f"{base_url}esearch.fcgi"

    params = {"db": "clinvar", "term": variant, "retmode": "json", "retmax": 10}

    response = api_session.request("get", esearch_url, params=params)
    esearch_results = response.json()

    clinvar_ids = esearch_results.get("esearchresult", {}).get("idlist", [])
    return clinvar_ids


def fetch_clinvar_details_multiple(
    clinvar_ids: List[str],
    api_session: APISessionManager,
    rettype: str = "vcv",
    retmode: str = "json",
) -> Any:
    """
    Retrieve detailed information for multiple ClinVar IDs.

    Args:
        clinvar_ids: List of ClinVar IDs
        api_session: API session manager
        rettype: Return type - "vcv" or "clinvarset"
        retmode: Return mode - "xml" or "json"

    Returns:
        ClinVar records data
    """
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
    efetch_url = f"{base_url}esummary.fcgi"

    params = {
        "db": "clinvar",
        "id": ",".join(clinvar_ids),  # Comma-separated list of IDs
        "rettype": rettype,
        "retmode": retmode,
    }

    response = api_session.request("get", efetch_url, params=params)

    if retmode == "json":
        try:
            return response.json()
        except json.JSONDecodeError:
            return response.text
    else:
        return response.text


def get_clinvar_result(
    variant: str, genes: List[str], api_session: APISessionManager
) -> List[List[str]]:
    """
    Get ClinVar result for a variant.

    Args:
        variant: Variant in HGVS format
        genes: List of genes in the region
        api_session: API session manager
    Returns:
        ClinVar result
    """
    clinvar_results = []
    # find clinvar ids for the variant
    clinvar_ids = search_clinvar_by_variant(variant, api_session)
    if clinvar_ids:
        # fetch clinvar details
        clinvar_details = fetch_clinvar_details_multiple(clinvar_ids, api_session)
        for clinvar_id in clinvar_ids:
            if "result" in clinvar_details and clinvar_id in clinvar_details["result"]:
                clinvar_detail = clinvar_details["result"][clinvar_id]
                if (
                    "obj_type" in clinvar_detail
                    and clinvar_detail["obj_type"] != "Haplotype"
                ):
                    classification = clinvar_detail.get("germline_classification", None)
                    description = classification.get("description", None)
                    review_status = classification.get("review_status", None)
                    clinvar_genes = clinvar_detail.get("genes", None)
                    if None not in [
                        classification,
                        description,
                        review_status,
                        clinvar_genes,
                    ]:
                        this_gene = [
                            a["symbol"] for a in clinvar_genes if a["symbol"] in genes
                        ]
                        if this_gene != []:
                            this_gene = this_gene[0]
                            clinvar_results.append(
                                [
                                    clinvar_id,
                                    this_gene,
                                    description,
                                    review_status,
                                ]
                            )
    return clinvar_results


def get_vep_result(variant: str, api_session: APISessionManager) -> Any:
    """
    Get VEP result for a variant.

    Args:
        variant: Variant in HGVS format
        api_session: API session manager

    Returns:
        VEP result
    """
    headers = {
        "Host": "rest.ensembl.org",
        "Accept": "application/json, text/javascript, */*; q=0.01",
        "Accept-Language": "en-US,en;q=0.5",
        "Connection": "keep-alive",
    }
    url_pattern = "http://rest.ensembl.org/vep/Homo_sapiens/hgvs/{variant}"
    variant_quoted = urllib.parse.quote(variant)
    url = url_pattern.format(variant=variant_quoted)
    response = api_session.request("get", url, headers=headers)
    json_text = response.json()
    return json_text


def filter_vep_result(vep_result, transcripts):
    """
    Filter VEP result to transcripts, find consequences.

    Args:
        vep_result: VEP result
        transcripts: List of transcripts

    Returns:
        List of consequences, boolean for high impact
    """
    consequence_terms = []
    filter_to_intranscript = []
    if "transcript_consequences" in vep_result[0]:
        filter_to_intranscript = [
            a
            for a in vep_result[0]["transcript_consequences"]
            if a.get("transcript_id", None) in transcripts
        ]
    high_impact = False
    for transcript_consequence in filter_to_intranscript:
        gene_symbol = transcript_consequence.get("gene_symbol", None)
        transcript_id = transcript_consequence.get("transcript_id", None)
        var_impact = transcript_consequence.get("impact", None)
        var_consequence = transcript_consequence.get("consequence_terms", None)
        if None not in [gene_symbol, transcript_id, var_consequence, var_impact]:
            consequence_terms.append(
                [
                    gene_symbol,
                    transcript_id,
                    var_consequence,
                    var_impact,
                ]
            )
            if var_impact == "HIGH":
                high_impact = True

    return consequence_terms, high_impact


def check_boundary(haplotype_boundary, transcript_boundary):
    """
    Check if a haplotype boundary overlaps with a transcript boundary.

    Args:
        haplotype_boundary: Haplotype boundary in the format of "start-end"
        transcript_boundary: Transcript boundary in the format of a tuple (transcript_id, start, end)

    Returns:
        "fully_contained", "partial", "truncated", or "no_overlap"
    """
    haplotype_boundary_split = haplotype_boundary.split("-")
    haplotype_start = int(haplotype_boundary_split[0].split("truncated")[0])
    haplotype_end = int(haplotype_boundary_split[1].split("truncated")[0])
    transcript_start = transcript_boundary[1]
    transcript_end = transcript_boundary[2]
    if haplotype_start <= transcript_start and haplotype_end >= transcript_end:
        return "fully_contained"
    if haplotype_start < transcript_end and transcript_start < haplotype_end:
        if (
            haplotype_start > transcript_start
            and "truncated" in haplotype_boundary_split[0]
        ) or (
            haplotype_end < transcript_end
            and "truncated" in haplotype_boundary_split[1]
        ):
            return "truncated"
        return "partial"
    return "no_overlap"


def modify_hapbound(bound1, bound2, truncated):
    """
    Get haplotype boundaries to appear in vcf

    Args:
        bound1: Start of the haplotype boundary
        bound2: End of the haplotype boundary
        truncated: truncated direction, 5 prime or 3 prime or both

    Returns:
        Haplotype boundary in the format of "start-end"
    """
    hap_bound = f"{bound1}-{bound2}"
    if truncated == ["5p"]:
        hap_bound = f"{bound1}truncated-{bound2}"
    elif truncated == ["3p"]:
        hap_bound = f"{bound1}-{bound2}truncated"
    elif truncated == ["5p", "3p"]:
        hap_bound = f"{bound1}truncated-{bound2}truncated"
    return hap_bound


def get_json_info(json_file: str, region_name: str) -> Optional[RegionInfo]:
    """
    Get information from the JSON file.

    Args:
        json_file: Path to the JSON file
        region_name: Name of the region

    Returns:
        RegionInfo object or None if region not found
    """
    report_info = None
    region_depth = None
    haplotype_boundary = {}
    known_nonfunctional_haplotypes = []

    with open(json_file, "r") as f:
        json_data = json.load(f)
        if region_name not in json_data:
            return None

        gene_call = json_data[region_name]
        if gene_call.get("final_haplotypes") is None:
            return None

        # Extract haplotype boundaries
        if gene_call.get("haplotype_details") is not None:
            for hap, hap_info in gene_call["haplotype_details"].items():
                bound = hap_info["boundary"]
                truncated = hap_info["is_truncated"]
                if bound is not None and len(bound) == 2:
                    bound1 = bound[0]
                    bound2 = bound[1]
                    haplotype_boundary.setdefault(
                        hap, modify_hapbound(bound1, bound2, truncated)
                    )
                    if (
                        gene_call["two_copy_haplotypes"] is not None
                        and hap in gene_call["two_copy_haplotypes"]
                    ):
                        haplotype_boundary.setdefault(
                            f"{hap}_cp2", modify_hapbound(bound1, bound2, truncated)
                        )

        # Get region depth
        region_depth = None
        if region_name != "CFHclust":
            region_depth = gene_call["region_depth"]["median"]
        elif "CFH" in json_data:
            region_depth = json_data["CFH"]["region_depth"]["median"]

        # Region-specific information extraction
        if region_name == "smn1":
            known_nonfunctional_haplotypes = list(
                gene_call["smn_del78_haplotypes"].values()
            )
            report_info = (
                f"smn1_cn: {gene_call['smn1_cn']}; "
                f"smn2_cn: {gene_call['smn2_cn']}; "
                f"smn_del78_cn: {gene_call['smn_del78_cn']}"
            )
        elif region_name == "rccx":
            annotated_alleles = "-".join(
                [
                    "'" + str(a).replace(",", "|") + "'"
                    for a in gene_call["annotated_alleles"]
                ]
            )
            report_info = (
                f"Phasing success: {gene_call['phasing_success']}; "
                f"annotated_alleles: {annotated_alleles}"
            )
        elif region_name == "hba":
            known_nonfunctional_haplotypes = [
                a for a in gene_call["final_haplotypes"].values() if "homology" in a
            ]
            report_info = (
                f"genotype: {gene_call['genotype']}; "
                f"sv_called: {gene_call['sv_called']}"
            )
        elif region_name == "ikbkg":
            known_nonfunctional_haplotypes = gene_call["deletion_haplotypes"]
            report_info = f"deletion_haplotypes: {gene_call['deletion_haplotypes']}"
        elif region_name in ["GBA", "CFHclust"]:
            fusion = gene_call["fusions_called"]
            if fusion == {}:
                report_info = f"fusions_called: {fusion}"
            else:
                fusion_info = ""
                for hap, hap_fusion_info in fusion.items():
                    breakpoint = hap_fusion_info["breakpoint"]
                    breakpoint1 = "-".join([str(a) for a in breakpoint[0]])
                    breakpoint2 = "-".join([str(a) for a in breakpoint[1]])
                    fusion_info += (
                        f"{hap}:{hap_fusion_info['type']}_{breakpoint1}_{breakpoint2};"
                    )
                report_info = f"fusions_called: {fusion_info}"
        elif region_name == "f8":
            report_info = f"sv_called: {gene_call['sv_called']}"
        elif region_name in ["strc", "ncf1"]:
            report_info = f"gene_cn: {gene_call['gene_cn']}"
        elif region_name == "opn1lw":
            report_info = (
                f"Phasing success: {gene_call['phasing_success']}; "
                f"annotated_alleles: {gene_call['annotated_alleles']}"
            )

    return RegionInfo(
        report_info=report_info,
        region_depth=region_depth,
        haplotype_boundary=haplotype_boundary,
        known_nonfunctional_haplotypes=known_nonfunctional_haplotypes,
    )


# ============================================================================
# HELPER FUNCTIONS FOR MAIN PROCESSING
# ============================================================================


def setup_logging(log_level: str) -> None:
    """Set up logging configuration."""
    log_level_attr = getattr(logging, log_level.upper(), logging.INFO)
    logging.basicConfig(
        level=log_level_attr,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    )

    # Disable verbose urllib3/requests connection pool logging
    logging.getLogger("urllib3.connectionpool").setLevel(logging.WARNING)
    logging.getLogger("urllib3").setLevel(logging.WARNING)
    logging.getLogger("requests").setLevel(logging.WARNING)


def load_config(
    config_file: Optional[str], api_session: APISessionManager
) -> Dict[str, Any]:
    """
    Load Paraphase configuration from file or default URL.

    Args:
        config_file: Path to config file, or None to use default
        api_session: API session manager

    Returns:
        Configuration dictionary
    """
    if config_file is None:
        logging.info("Using default paraphase config")
        response = api_session.request("get", DEFAULT_CONFIG_URL)
        config = yaml.safe_load(response.text)
    else:
        if not os.path.exists(config_file):
            raise FileNotFoundError(f"Config file {config_file} not found")
        with open(config_file, "r") as f:
            config = yaml.safe_load(f)
    return config


def load_transcript_mappings(
    transcript_file: str,
) -> Tuple[Dict[str, List], Dict[str, List]]:
    """
    Load transcript to gene mappings from file.

    Args:
        transcript_file: Path to transcript file

    Returns:
        Tuple of (gene_transcript, transcript_gene) dictionaries
    """
    gene_transcript = {}
    with open(transcript_file) as f:
        for line in f:
            parts = line.split()
            gene_transcript.setdefault(
                parts[1], [parts[0], int(parts[2]), int(parts[3])]
            )
    return gene_transcript


def format_variant(chrom: str, pos: int, ref: str, alt: str) -> Optional[str]:
    """
    Format variant in HGVS format.

    Args:
        chrom: Chromosome
        pos: Position
        ref: Reference allele
        alt: Alternate allele

    Returns:
        Variant in HGVS format or None if cannot be formatted
    """
    # SNV
    if len(ref) == 1 and len(alt) == 1:
        return f"{chrom}:g.{pos}{ref}>{alt}"
    # insertion
    elif len(ref) == 1 and len(alt) > 1:
        return f"{chrom}:g.{pos}_{pos+1}ins{alt[1:]}"
    # deletion
    elif len(ref) > 1 and len(alt) == 1:
        del_len = len(ref) - 1
        return f"{chrom}:g.{pos+1}_{pos+del_len}del"
    return None


def process_vcf_file(
    vcf_path: str,
    genes: List[str],
    transcripts: List[str],
    region_info: RegionInfo,
    api_session: APISessionManager,
) -> Tuple[Dict[str, List[str]], Dict[Tuple[str, str], List[str]], Dict[str, str]]:
    """
    Process VCF file to extract variants and annotations.

    Args:
        vcf_path: Path to VCF file
        genes: List of genes in region
        transcripts: List of transcripts
        region_info: Region information from JSON
        api_session: API session manager

    Returns:
        Tuple of (haplotype_variants, variant_info, haplotype_boundary)
    """
    haplotype_variants = {}
    variant_info = {}
    haplotype_boundary = region_info.haplotype_boundary.copy()

    if genes == []:
        return haplotype_variants, variant_info, haplotype_boundary

    with open(vcf_path) as f:
        header = None
        nhap = 0

        for line in f:
            if line.startswith("##"):
                continue

            if line.startswith("#"):
                header = line.split()
                nhap = len(header) - 9
                for i in range(nhap):
                    hap_name = header[i + 9]
                    haplotype_variants.setdefault(hap_name, [])
                continue

            parts = line.split()
            if parts[6] != "PASS":
                continue

            chrom = parts[0].replace("chr", "")
            pos = int(parts[1])
            ref = parts[3]
            alt = parts[4]

            if "<" in alt:
                logging.debug(f"Skipping large SV: {chrom}:{pos}_{ref}_{alt}")
                continue

            all_bounds = parts[7].split(";")[0].split("=")[1].split(",")
            variant = format_variant(chrom, pos, ref, alt)

            if variant is None:
                continue

            haps_with_variant = []
            for i in range(nhap):
                this_hap = parts[i + 9]
                gt = this_hap.split(":")[0]
                hap_name = header[i + 9]

                if hap_name in haplotype_boundary:
                    haplotype_boundary[hap_name] = all_bounds[i]
                else:
                    haplotype_boundary.setdefault(hap_name, all_bounds[i])

                if gt == "1":
                    haplotype_variants.setdefault(hap_name, []).append(variant)
                    haps_with_variant.append(hap_name)

            haps_with_variant_filtered = [
                a
                for a in haps_with_variant
                if a not in region_info.known_nonfunctional_haplotypes
            ]

            if not haps_with_variant_filtered:
                logging.debug(
                    f"Skipping variant: {variant} because it's only in known nonfunctional haplotypes"
                )
                continue

            logging.debug(f"Annotating variant: {variant}")
            # check VEP
            vep_result = get_vep_result(variant, api_session)
            consequence_terms, high_impact = filter_vep_result(vep_result, transcripts)

            # Get colocated variants at this position, find clinvar entries in VEP result
            clin_sig = []
            if "colocated_variants" in vep_result[0]:
                colocated_variants = vep_result[0]["colocated_variants"]
                for colocated_variant in colocated_variants:
                    if "clin_sig" in colocated_variant:
                        clin_sig.extend(colocated_variant["clin_sig"])

            pathogenic_clin_sig = [
                a for a in clin_sig if a in ["pathogenic", "likely_pathogenic"]
            ]

            clinvar_results = []
            pathogenic_clinvar_results = []
            if pathogenic_clin_sig or high_impact:
                # check ClinVar through ClinVar API
                clinvar_results = get_clinvar_result(variant, genes, api_session)
                for clinvar_result in clinvar_results:
                    if "athogenic" in clinvar_result[2] and (
                        (
                            "o conflicts" in clinvar_result[3]
                            and "multiple" in clinvar_result[3]
                        )
                        or variant in CLINVAR_WHITELIST
                    ):
                        pathogenic_clinvar_results.append(clinvar_result)

            if high_impact or pathogenic_clinvar_results != []:
                logging.debug(
                    f"Variant {variant} is high impact or has clinvar entries indicating pathogenicity"
                )

                # Process ClinVar results
                for clinvar_result in clinvar_results:
                    clinvar_id = clinvar_result[0]
                    clinvar_gene = clinvar_result[1]
                    clinvar_description = clinvar_result[2]

                    if clinvar_gene in genes:
                        filtered_consequence_terms = [
                            a for a in consequence_terms if a[0] == clinvar_gene
                        ]

                        variant_key = (variant, clinvar_gene)
                        if variant_key not in variant_info:
                            variant_info[variant_key] = []

                        this_info = f"ClinVar-id{clinvar_id};{clinvar_description};"
                        variant_info[variant_key].append(this_info)

                        for a in filtered_consequence_terms:
                            if a[0] == clinvar_gene:
                                this_transcript = a[1]
                                this_consequence = a[2]
                                this_impact = a[3]
                                this_info = (
                                    f"{this_transcript};{';'.join(this_consequence)};"
                                    f"{this_impact};"
                                )
                                if this_info not in variant_info[variant_key]:
                                    variant_info[variant_key].append(this_info)

                # Process high impact variants with no clinvar results
                if clinvar_results == []:
                    for a in consequence_terms:
                        this_gene = a[0]
                        this_transcript = a[1]
                        this_consequence = a[2]
                        this_impact = a[3]

                        if this_impact == "HIGH" and this_gene in genes:
                            variant_key = (variant, this_gene)
                            if variant_key not in variant_info:
                                variant_info[variant_key] = []

                            this_info = (
                                f"{this_transcript};{';'.join(this_consequence)};"
                                f"{this_impact};"
                            )
                            if this_info not in variant_info[variant_key]:
                                variant_info[variant_key].append(this_info)

    return haplotype_variants, variant_info, haplotype_boundary


def write_region_output(
    region: str,
    genes: List[str],
    region_info: RegionInfo,
    haplotype_variants: Dict[str, List[str]],
    variant_info: Dict[Tuple[str, str], List[str]],
    haplotype_boundary: Dict[str, str],
    gene_transcript: Dict[str, List],
    summary_out,
    detailed_out,
) -> None:
    """
    Write output for a region to summary and detailed files.

    Args:
        region: Region name
        genes: List of genes
        region_info: Region information
        haplotype_variants: Variants per haplotype
        variant_info: Variant annotation information
        haplotype_boundary: Haplotype boundaries
        gene_transcript: Gene to transcript mapping
        summary_out: Summary file handle
        detailed_out: Detailed file handle
    """
    gene_variants = {}
    genes_exclude = set()

    # Determine which genes to exclude
    for gene in genes:
        gene_variants.setdefault(gene, [])
        nhap = 0
        for hap in haplotype_boundary:
            boundary_check = check_boundary(
                haplotype_boundary[hap], gene_transcript[gene]
            )
            if boundary_check != "no_overlap":
                nhap += 1
        if nhap == 0 and gene not in ALWAYS_REPORT_GENES:
            genes_exclude.add(gene)
            logging.debug(
                f"Gene {gene} has no overlap with any haplotype and is excluded from reporting. This often means that its paralog is already included in reporting."
            )

    # Map variants to genes
    for variant, gene in variant_info:
        gene_variants.setdefault(gene, []).append(variant)

    # Write output for each gene
    for gene in gene_variants:
        if gene in genes_exclude:
            continue

        nhap = 0
        nonfunctional = 0

        for hap in haplotype_variants:
            boundary_check = check_boundary(
                haplotype_boundary[hap], gene_transcript[gene]
            )
            pathogenic_variants = list(
                set(haplotype_variants[hap]).intersection(set(gene_variants[gene]))
            )
            pathogenic_variants_annotated = [
                a + "=" + ":".join(variant_info[(a, gene)]) for a in pathogenic_variants
            ]

            if hap not in region_info.known_nonfunctional_haplotypes:
                detailed_line = [
                    str(a)
                    for a in [
                        region,
                        gene,
                        gene_transcript[gene][0],
                        gene_transcript[gene][1],
                        gene_transcript[gene][2],
                        hap,
                        haplotype_boundary[hap],
                        boundary_check,
                        ",".join(pathogenic_variants_annotated),
                    ]
                ]
            else:
                detailed_line = [
                    str(a)
                    for a in [
                        region,
                        gene,
                        gene_transcript[gene][0],
                        gene_transcript[gene][1],
                        gene_transcript[gene][2],
                        hap,
                        haplotype_boundary[hap],
                        boundary_check,
                        "known_nonfunctional",
                    ]
                ]

            detailed_out.write(",".join(detailed_line) + "\n")

            if boundary_check != "no_overlap":
                nhap += 1
                if (
                    boundary_check == "truncated"
                    or pathogenic_variants != []
                    or hap in region_info.known_nonfunctional_haplotypes
                ):
                    nonfunctional += 1

        report_info = region_info.report_info or ""
        summary_line = [
            str(a)
            for a in [
                region,
                region_info.region_depth,
                gene,
                nhap,
                nonfunctional,
                nhap - nonfunctional,
                report_info.replace(",", ";"),
            ]
        ]
        summary_out.write(",".join(summary_line) + "\n")


def main():
    parser = load_parameters()
    args = parser.parse_args()

    # Set up logging
    setup_logging(args.log_level)

    # Set up API session with rate limiting
    api_session = APISessionManager(min_request_interval=args.api_delay)
    if args.api_delay > 0:
        logging.info(
            f"Rate limiting enabled: minimum {args.api_delay:.2f} seconds between API requests"
        )

    # Parse input parameters
    json_file = args.input
    logging.info(f"Input JSON file: {json_file}")
    sample_id = os.path.basename(json_file).split(".paraphase.json")[0]

    # Load transcript mappings and config
    gene_transcript = load_transcript_mappings(DEFAULT_TRANSCRIPT_FILE)
    try:
        config = load_config(args.config, api_session)
    except FileNotFoundError as e:
        logging.warning(str(e))
        return

    regions_to_annotate = args.region
    if regions_to_annotate is None:
        logging.info("No region specified, annotating medical regions")
        regions_to_annotate = MEDICAL_REGIONS
    elif regions_to_annotate == "all":
        regions_to_annotate = [
            a for a in config.keys() if a not in ["CFH", "CFHR3"]
        ] + ["CFHclust"]
    else:
        regions_to_annotate = regions_to_annotate.split(",")

    # Set up output files
    outdir = args.out
    os.makedirs(outdir, exist_ok=True)
    output_summary_file = os.path.join(outdir, f"{sample_id}_summary.csv")
    output_detailed_file = os.path.join(outdir, f"{sample_id}_detailed.csv")

    for each_file in [output_summary_file, output_detailed_file]:
        if os.path.exists(each_file):
            logging.warning(f"File {each_file} already exists, removing it")
            os.remove(each_file)

    summary_out = open(output_summary_file, "w")
    summary_out.write(",".join(SUMMARY_HEADER) + "\n")
    detailed_out = open(output_detailed_file, "w")
    detailed_out.write(",".join(DETAILED_HEADER) + "\n")

    # Check VCF directory
    vcf_dir = os.path.join(
        os.path.dirname(json_file),
        f"{sample_id}_paraphase_vcfs",
    )
    if not os.path.exists(vcf_dir):
        raise FileNotFoundError(f"VCF directory {vcf_dir} does not exist")

    # Process each region
    for region in regions_to_annotate:
        if region not in config and region.upper() in config:
            region = region.upper()
        if region not in config and region.lower() in config:
            region = region.lower()

        logging.info(f"Annotating region: {region}")
        # Get region information from JSON
        region_info = get_json_info(json_file, region)
        if region_info is None:
            logging.warning(
                f"Skipping region {region} because it is not found or is no-call due to low coverage in JSON file"
            )
            continue

        if region not in REGIONS_NO_VCF:
            if region not in config:
                logging.warning(f"Region {region} not found in config file")
                continue

            # Get genes for this region
            genes = config[region]["genes"].split(",")
            if region == "rccx":
                # TNXB is annotated in its own region
                genes.remove("TNXB")

            logging.info(f"Genes in region {region}: {','.join(genes)}")

            # Get transcripts for genes
            transcripts = []
            for gene in genes:
                transcript = gene_transcript[gene][0]
                logging.info(
                    f"Gene: {gene}, transcript: {transcript}, "
                    f"start: {gene_transcript[gene][1]}, end: {gene_transcript[gene][2]}"
                )
                transcripts.append(transcript)

            # Process VCF file
            vcf = os.path.join(vcf_dir, f"{sample_id}_{region}.vcf")
            logging.info(f"Annotating VCF file: {vcf}")
            if not os.path.exists(vcf):
                logging.warning(f"VCF file {vcf} not found")
                continue
            haplotype_variants, variant_info, haplotype_boundary = process_vcf_file(
                vcf, genes, transcripts, region_info, api_session
            )

            # Write output for this region
            write_region_output(
                region,
                genes,
                region_info,
                haplotype_variants,
                variant_info,
                haplotype_boundary,
                gene_transcript,
                summary_out,
                detailed_out,
            )

        # Add summary for regions that do not need VCF
        else:
            summary_line = [
                str(a)
                for a in [
                    region,
                    region_info.region_depth,
                    "NA",
                    "NA",
                    "NA",
                    "NA",
                    region_info.report_info.replace(",", ";") or "",
                ]
            ]
            summary_out.write(",".join(summary_line) + "\n")

    summary_out.close()
    detailed_out.close()
    logging.info(
        f"Completed annotating all regions. Output files: {output_summary_file} and {output_detailed_file}"
    )


if __name__ == "__main__":
    main()
