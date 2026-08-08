use paraphase::io::json::ParsedParaphaseOutputJSON;
use paraphase::toolkit::low_complexity::LowConfidenceSites;
use paraphase::toolkit::util::DResult;
use paraphase::{
    config, depth, io,
    phaser::{self, Phaser},
    toolkit::util::{self, test_file},
};

use itertools::Itertools;

use std::collections::{BTreeMap, BTreeSet};

fn dict_matches_with_gaps(x: &BTreeMap<String, String>, y: &BTreeMap<String, String>) -> bool {
    use util::DeletionInsensitiveCompare;
    for (k, v) in x {
        let yv = y.get(k).unwrap();
        if !v.same_without_dels(yv) {
            return false;
        }
    }
    y.keys().all(|y| x.contains_key(y))
}

fn slice_matches_with_gaps(x: &[impl AsRef<[u8]>], y: &[impl AsRef<[u8]>]) -> bool {
    use util::DeletionInsensitiveCompare;
    x.iter()
        .all(|z| y.iter().any(|w| w.as_ref().same_without_dels(&z.as_ref())))
        && y.iter()
            .all(|z| x.iter().any(|w| w.as_ref().same_without_dels(&z.as_ref())))
}

#[test]
fn smn1_ok() -> DResult {
    util::init_log(log::LevelFilter::Info);

    let outdir = tempfile::TempDir::new()?;
    let genome_bam = test_file("bams/HG00733.smn1.bam");
    let genome_path = if let Ok(x) = std::env::var("HG38") {
        x.trim_end_matches(".mmi").to_string()
    } else {
        log::warn!("Skipping test: HG38 env is not set.");
        return Ok(());
    };
    let gene_name = "smn1";
    let depth = depth::Result {
        median: 35.0,
        median_absolute_difference: 0.1,
        sex: depth::Sex::Other,
    };
    let settings = phaser::Settings::new(
        "HG00733",
        (genome_path, genome_bam),
        outdir.path(),
        gene_name,
        &config::Region::try_load(None)?,
        /* genome depth= */ Some(depth),
        /* sex = */ None,
        String::from("38"),
        None,
        0.03,
        false,
    );

    let gene_config = config::Gene::try_load(None)?;
    let mut phaser = Phaser::new(
        settings,
        Some(gene_config),
        None, // Option<SiteSelectionSettings>
        None, // Option<RealignSettings>
    )?;
    let call = phaser.run()?;

    let json_to_match = test_file("jsons/HG00733.json.xz");
    let region_config = config::Region::try_load(None)?;
    let json_data = ParsedParaphaseOutputJSON::from_path(&json_to_match, Some(&region_config))?;
    let gene_data = json_data
        .gene_data
        .get(gene_name)
        .ok_or("Missing gene data in expected json")?;

    log::debug!("Depth: {:?}", phaser.region_avg_depth);
    log::debug!("Call: {call:?}");
    let sites: Vec<String> = gene_data
        .sites_for_phasing
        .iter()
        .map(std::string::ToString::to_string)
        .sorted()
        .collect::<Vec<_>>();
    let paraph_rs_sites = call
        .sites_for_phasing
        .iter()
        .cloned()
        .sorted()
        .collect::<Vec<_>>();
    log::trace!("Python sites: {sites:?}");
    log::trace!("Rust sites: {paraph_rs_sites:?}");

    assert_eq!(
        sites, paraph_rs_sites,
        "Python sites: {sites:?}\nRust sites: {paraph_rs_sites:?}\n"
    );

    Ok(())
}

// Same-input oracle: Paraphase v3.5.0 c8016dff40d105501868719bdae0defd9937abfe
// on the public HPRC HG00733 GRCh38 HiFi locus fixture (m54329U_2019).
// Genome depth 31 comes from the same source BAM's bundled depth-probe union.
#[test]
fn strc_ok() -> DResult {
    util::init_log(log::LevelFilter::Info);

    let outdir = tempfile::TempDir::new()?;
    let genome_bam = test_file("bams/HG00733.strc.bam");
    let genome_path = if let Ok(x) = std::env::var("HG38") {
        x.trim_end_matches(".mmi").to_string()
    } else {
        log::warn!("Skipping test: HG38 env is not set.");
        return Ok(());
    };
    let gene_name = "strc";
    let depth = depth::Result {
        median: 31.0,
        median_absolute_difference: 0.1,
        sex: depth::Sex::Other,
    };
    let region_config = config::Region::try_load(None)?;
    let settings = phaser::Settings::new(
        "HG00733",
        (genome_path, genome_bam),
        outdir.path(),
        gene_name,
        &region_config,
        /* genome depth= */ Some(depth),
        /* sex = */ None,
        String::from("38"),
        None,
        0.03,
        false,
    );

    let gene_config = config::Gene::try_load(None)?;
    let mut phaser = Phaser::new(
        settings,
        Some(gene_config),
        None, // Option<SiteSelectionSettings>
        None, // Option<RealignSettings>
    )?;
    let call = phaser.run()?;
    assert_eq!(phaser.region_avg_depth[0], (71.0f32, 82.0f32));

    let json_to_match = test_file("jsons/HG00733.strc.v3.5.0.json.xz");
    let json_data = ParsedParaphaseOutputJSON::from_path(&json_to_match, Some(&region_config))?;
    let expected_call = json_data
        .object
        .get(gene_name)
        .and_then(serde_json::Value::as_object)
        .ok_or("Missing STRC call in expected json")?;
    let gene_data = json_data
        .gene_data
        .get(gene_name)
        .ok_or("Missing STRC gene data in expected json")?;

    assert_eq!(call.total_cn, Some(4));
    assert_eq!(
        call.region_specific_info
            .get("gene_cn")
            .and_then(serde_json::Value::as_i64),
        Some(2)
    );
    assert_eq!(
        call.region_specific_info
            .get("intergenic_depth")
            .and_then(serde_json::Value::as_f64),
        Some(31.0)
    );

    let expected_final_haplotypes = serde_json::from_value::<BTreeMap<String, String>>(
        expected_call
            .get("final_haplotypes")
            .ok_or("Missing STRC final_haplotypes in expected json")?
            .clone(),
    )?;
    let haplotype_class = |name: &str| {
        if name.contains("strcp1") {
            "strcp1"
        } else {
            "strc"
        }
    };
    let final_classes_by_sequence = |haplotypes: &BTreeMap<String, String>| {
        haplotypes
            .iter()
            .map(|(sequence, name)| (sequence.clone(), haplotype_class(name)))
            .collect::<BTreeMap<_, _>>()
    };
    assert_eq!(
        final_classes_by_sequence(&call.final_haplotypes),
        final_classes_by_sequence(&expected_final_haplotypes)
    );

    let expected_two_copy_haplotypes = serde_json::from_value::<Vec<String>>(
        expected_call
            .get("two_copy_haplotypes")
            .ok_or("Missing STRC two_copy_haplotypes in expected json")?
            .clone(),
    )?;
    let two_copy_classes = |haplotypes: &[String]| {
        haplotypes
            .iter()
            .map(|name| haplotype_class(name))
            .sorted()
            .collect::<Vec<_>>()
    };
    assert_eq!(
        two_copy_classes(&call.two_copy_haplotypes),
        two_copy_classes(&expected_two_copy_haplotypes)
    );

    let expected_sites = gene_data
        .sites_for_phasing
        .iter()
        .map(std::string::ToString::to_string)
        .sorted()
        .collect::<Vec<_>>();
    let found_sites = call
        .sites_for_phasing
        .iter()
        .cloned()
        .sorted()
        .collect::<Vec<_>>();
    assert_eq!(expected_sites, found_sites);

    let expected_haplotypes = gene_data
        .assembled_haps
        .as_ref()
        .ok_or("Missing STRC assembled_haplotypes in expected json")?
        .iter()
        .map(std::string::ToString::to_string)
        .sorted()
        .collect::<Vec<_>>();
    let found_haplotypes = call
        .assembled_haplotypes
        .iter()
        .cloned()
        .sorted()
        .collect::<Vec<_>>();
    assert_eq!(expected_haplotypes, found_haplotypes);

    let expected_haplotype_details =
        serde_json::from_value::<BTreeMap<String, paraphase::phaser::HapInfoForJson>>(
            expected_call
                .get("haplotype_details")
                .ok_or("Missing STRC haplotype_details in expected json")?
                .clone(),
        )?;
    let variants_by_sequence =
        |haplotypes: &BTreeMap<String, String>,
         details: &BTreeMap<String, paraphase::phaser::HapInfoForJson>| {
            haplotypes
                .iter()
                .map(|(sequence, name)| {
                    let variants = details
                        .get(name)
                        .expect("Missing STRC haplotype details")
                        .variants
                        .iter()
                        .cloned()
                        .sorted()
                        .collect::<Vec<_>>();
                    (sequence.clone(), variants)
                })
                .collect::<BTreeMap<_, _>>()
        };
    assert_eq!(
        variants_by_sequence(&call.final_haplotypes, &call.haplotype_details),
        variants_by_sequence(&expected_final_haplotypes, &expected_haplotype_details)
    );

    Ok(())
}

#[test]
fn amy1_ok() -> DResult {
    util::init_log(log::LevelFilter::Info);
    let gene = "AMY1A";
    let sample = "HG001";
    let json_to_match = util::test_file("jsons/HG001.json.xz");
    let genome_bam = util::test_file("bams/HG001.amy1a.bam");
    let outdir = tempfile::TempDir::new()?;
    let genome_path = if let Ok(x) = std::env::var("HG38") {
        x.trim_end_matches(".mmi").to_string()
    } else {
        log::warn!("Skipping test: HG38 env is not set.");
        return Ok(());
    };
    let region_config = config::Region::try_load(None).unwrap();
    let json_data =
        io::json::ParsedParaphaseOutputJSON::from_path(&json_to_match, Some(&region_config))
            .unwrap();
    let amy1_data = json_data.gene_data.get(gene).unwrap();
    let depth = depth::Calculator::from_hg38(genome_bam.display().to_string(), None)?.compute();
    log::debug!("Gene data: {amy1_data:?}");
    log::debug!("depth: {depth:?}");
    assert!(amy1_data.pivot_site.is_none());
    log::debug!(
        "Building Settings struct for HG001 with paths {genome_path:?}/{genome_bam:?} with outdir = {outdir:?}"
    );
    let settings = phaser::Settings::new(
        sample,
        (genome_path, genome_bam),
        outdir.path(),
        gene,
        &region_config,
        /* genome depth= */ None,
        /* sex = */ None,
        String::from("38"),
        None,
        0.03,
        false,
    );
    log::debug!("Building gene-level config");
    let gene_config = config::Gene::try_load(None).unwrap();
    log::debug!("Building Phaser.");
    let mut phaser = Phaser::new(
        settings.clone(),
        Some(gene_config.clone()),
        None, // Option<SiteSelectionSettings>
        None, // Option<RealignSettings>
    )?;
    log::debug!("Built Phaser.");
    let call = phaser.run()?;

    assert_eq!(phaser.region_avg_depth[0], (71.0f32, 78.0f32));
    let sites: Vec<String> = amy1_data
        .sites_for_phasing
        .iter()
        .map(std::string::ToString::to_string)
        .sorted()
        .collect::<Vec<_>>();
    let found_sites = call
        .sites_for_phasing
        .iter()
        .cloned()
        .sorted()
        .collect::<Vec<_>>();
    assert_eq!(
        sites, found_sites,
        "Found sites: {found_sites:?}. Expected {sites:?}"
    );
    let read_to_hap = amy1_data
        .read_to_hap
        .iter()
        .map(|(k, v)| (k.read_name.clone(), v.to_string()))
        .collect::<BTreeMap<_, _>>();
    let read_details = call
        .read_details
        .clone()
        .into_iter()
        .map(|(k, v)| (k, vstr::VStr::from(&v).to_string()))
        .collect::<BTreeMap<String, String>>();
    assert!(
        dict_matches_with_gaps(&read_details, &read_to_hap),
        "Found {read_details:?} but expected {read_to_hap:?}"
    );
    let expected_haps = amy1_data
        .assembled_haps
        .as_ref()
        .unwrap()
        .iter()
        .map(std::string::ToString::to_string)
        .sorted()
        .collect::<Vec<String>>();
    let assembled_haps = call
        .assembled_haplotypes
        .iter()
        .map(std::string::ToString::to_string)
        .sorted()
        .collect::<Vec<_>>();
    assert!(
        slice_matches_with_gaps(&expected_haps, &assembled_haps),
        "Found {assembled_haps:?} but expected {expected_haps:?}"
    );
    Ok(())
}

/// Integration test for AGAP9
#[test]
fn agap9_ok() -> DResult {
    util::init_log(log::LevelFilter::Info);

    let json_to_match = util::test_file("jsons/HG00733.json.xz");
    let genome_bam = util::test_file("bams/HG00733.agap9.bam");
    let outdir = tempfile::TempDir::new()?;
    let genome_path = if let Ok(x) = std::env::var("HG38") {
        x.trim_end_matches(".mmi").to_string()
    } else {
        log::warn!("Skipping test: HG38 env is not set.");
        return Ok(());
    };
    let gene_name = "AGAP9";
    let region_config = config::Region::try_load(None)?;
    let json_data = ParsedParaphaseOutputJSON::from_path(&json_to_match, Some(&region_config))?;
    let gene_data = json_data
        .gene_data
        .get(gene_name)
        .ok_or("Missing AGAP9 data in expected json")?;
    let depth = depth::Calculator::from_hg38(genome_bam.display().to_string(), None)?.compute();
    log::debug!("Gene data: {gene_data:?}");
    log::debug!("depth: {depth:?}");
    assert!(gene_data.pivot_site.is_none());
    assert!(gene_data.pivot_site.is_none());

    log::debug!("Building Settings struct for BCH-35 with paths {genome_path:?}/{genome_bam:?} with outdir = {outdir:?}");
    let settings = phaser::Settings::new(
        "HG00733",
        (genome_path, genome_bam),
        outdir.path(),
        gene_name,
        &region_config,
        /* genome depth= */ None,
        /* sex = */ None,
        String::from("38"),
        None,
        0.03,
        false,
    );
    log::debug!("Building gene-level config");
    let gene_config = config::Gene::try_load(None)?;
    log::debug!("Building Phaser.");
    let mut phaser = Phaser::new(
        settings.clone(),
        Some(gene_config.clone()),
        None, // Option<SiteSelectionSettings>
        None, // Option<RealignSettings>
    )?;
    log::debug!("Built Phaser.");
    let expected_homopolymer_data = {
        let seqs = util::seq_name_pairs(&util::test_file("ref/AGAP9_ref.fa"), true)?;
        assert_eq!(seqs[0].0, b"chr10_47501354_47524138");
        let seq = seqs
            .into_iter()
            .map(|x| x.1)
            .next()
            .ok_or("Missing AGAP9 sequence")?;
        assert_eq!(22785, seq.len(), "Found seq: \"{seq}\"");
        LowConfidenceSites::new(&seq[..], phaser.offset(), None)
    };
    log::trace!("Expected {expected_homopolymer_data:?} for low_complexity_sites");
    let call = phaser.run()?;
    assert_eq!(
        phaser.low_complexity_sites, expected_homopolymer_data,
        "Found {:?} expected {expected_homopolymer_data:?} when running in/outside of Phaser",
        phaser.low_complexity_sites
    );
    let expected_homopolymer_data = util::parse_homopolymers(&test_file("agap9-hpol-expected.txt"));
    let found_homopol = phaser.low_complexity_sites.clone();
    let all_positions = phaser
        .low_complexity_sites
        .keys()
        .chain(expected_homopolymer_data.keys())
        .copied()
        .collect::<BTreeSet<_>>();
    let mut hpol_fails = (0usize, 0usize, 0usize);
    let mut hpol_success = 0;
    for pos in all_positions {
        if !found_homopol.contains_key(&pos) {
            hpol_fails.0 += 1;
            log::trace!("Found does not contain key {pos}");
            continue;
        }
        if !expected_homopolymer_data.contains_key(&pos) {
            hpol_fails.1 += 1;
            log::trace!("Expected does not contain key {pos}");
            continue;
        }
        if found_homopol.get(&pos) != expected_homopolymer_data.get(&pos) {
            hpol_fails.2 += 1;
            log::trace!(
                "Differing values at key {pos}: {:?}/{:?}",
                found_homopol.get(&pos),
                expected_homopolymer_data.get(&pos)
            );
        } else {
            log::trace!(
                "Found matching values at key {pos}: {:?}/{:?}",
                found_homopol.get(&pos),
                expected_homopolymer_data.get(&pos)
            );
            hpol_success += 1;
        }
    }
    let hpol_sum = hpol_fails.0 + hpol_fails.1 + hpol_fails.2;
    log::debug!("Call: {call:?}");
    assert_eq!(hpol_sum, 0, "Missing hp sites in python: {}. Missing in rust: {}. Mismatches: {}. Successes: {hpol_success}", hpol_fails.0, hpol_fails.1, hpol_fails.2);
    assert_eq!(
        phaser.low_complexity_sites, expected_homopolymer_data,
        "Found {:?} expected {expected_homopolymer_data:?} when running in Python vs in Rust",
        phaser.low_complexity_sites
    );

    // check depth
    assert_eq!(phaser.region_avg_depth[0], (79.0f32, 90.0f32));

    let sites: Vec<String> = gene_data
        .sites_for_phasing
        .iter()
        .map(std::string::ToString::to_string)
        .sorted()
        .collect::<Vec<_>>();
    let found_sites = call
        .sites_for_phasing
        .iter()
        .cloned()
        .sorted()
        .collect::<Vec<_>>();
    assert_eq!(
        sites, found_sites,
        "Found {found_sites:?} instead of {sites:?}",
    );

    log::debug!("Sites for phasing: {sites:?}/{found_sites:?}");

    let expected_reads = gene_data
        .read_to_hap
        .keys()
        .map(std::string::ToString::to_string)
        .collect::<BTreeSet<_>>();
    let found_reads = call.read_details.keys().cloned().collect::<BTreeSet<_>>();
    if expected_reads != found_reads {
        let not_found = expected_reads.difference(&found_reads).collect::<Vec<_>>();
        let not_expected = found_reads.difference(&expected_reads).collect::<Vec<_>>();
        log::debug!("Not found: {not_found:?}. Not expected: {not_expected:?}");
    }

    let read_to_hap = gene_data
        .read_to_hap
        .iter()
        .map(|(k, v)| (k.read_name.clone(), v.to_string()))
        .collect::<BTreeMap<_, _>>();
    assert_eq!(call.read_details, read_to_hap);
    let expected_haps = gene_data
        .assembled_haps
        .as_ref()
        .ok_or("Missing assembled_haps for HSFY1 expected data")?
        .iter()
        .map(std::string::ToString::to_string)
        .sorted()
        .collect::<Vec<String>>();
    let assembled_haps = call
        .assembled_haplotypes
        .iter()
        .map(std::borrow::ToOwned::to_owned)
        .collect::<Vec<_>>();
    assert_eq!(
        expected_haps, assembled_haps,
        "Found {assembled_haps:?} but expected {expected_haps:?}"
    );
    Ok(())
}
