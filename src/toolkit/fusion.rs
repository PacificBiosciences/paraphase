use crate::phaser;
use crate::toolkit::site_selection::CandidateSite;
use crate::toolkit::util::DError;
use serde_json::json;
use serde_json::Value;
use std::collections::BTreeMap;
use vstr::{VStr, VString};

impl phaser::Phaser {
    #[allow(clippy::type_complexity)]
    /// Identify fusion haplotypes and derive fusion annotations.
    ///
    /// Returns renamed hap labels, inferred two-copy haplotypes, and structured
    /// `fusions_called` metadata for output serialization.
    pub fn find_fusion<'a>(
        &mut self,
        assembled_haps: &BTreeMap<VStr<'a>, String>,
        fusion_direction: String,
    ) -> Result<
        (
            BTreeMap<VStr<'a>, String>,
            Vec<String>,
            BTreeMap<String, BTreeMap<String, Value>>,
        ),
        DError,
    > {
        // get PSVs
        let fusion_gene_def_variants = self.parse_psv()?;
        let psv_variants = fusion_gene_def_variants.get(self.gene_name());
        // update two-copy haplotypes
        let (assembled_haps_renamed, two_cp_haps) =
            self.update_twp_cp_in_fusion_cases(assembled_haps)?;
        let mut fusions_called: BTreeMap<String, BTreeMap<String, Value>> = BTreeMap::new();
        for (hap, hap_name) in &assembled_haps_renamed {
            let first_base = hap.first().ok_or_else(|| {
                phaser::Exception::new(format!(
                    "Haplotype '{}' is unexpectedly empty while calling fusion",
                    hap_name
                ))
            })?;
            let last_base = hap.last().ok_or_else(|| {
                phaser::Exception::new(format!(
                    "Haplotype '{}' is unexpectedly empty while calling fusion",
                    hap_name
                ))
            })?;
            let first_base0 = *first_base == b'0';
            let last_base0 = *last_base == b'0';
            if *first_base != b'x' && *last_base != b'x' {
                if (first_base0 && !last_base0) || (!first_base0 && last_base0) {
                    let (new_hap, all_sites) =
                        self.new_hap_for_breakpoint(hap.into(), psv_variants)?;
                    let fusion_breakpoint_index =
                        get_fusion_breakpoint_index(hap.into(), new_hap.clone())?;
                    log::debug!(
                        "Computed fusion breakpoint index for hap {hap_name}: {:?}",
                        fusion_breakpoint_index
                    );
                    if let Some(fusion_breakpoint_index) = fusion_breakpoint_index {
                        let bp1 = all_sites[fusion_breakpoint_index].pos + 1;
                        let bp2 = self.get_range_in_other_gene(bp1, Some(1000));
                        let bp3 = all_sites[fusion_breakpoint_index - 1].pos + 1;
                        let bp4 = self.get_range_in_other_gene(bp3, Some(1000));
                        if let (Some(bp2), Some(bp4)) = (bp2, bp4) {
                            let fusion_type = get_fusion_type(&fusion_direction, hap.into())?;
                            fusions_called.entry(hap_name.clone()).or_default().insert(
                                String::from("type"),
                                fusion_type.map_or(Value::Null, Value::String),
                            );
                            fusions_called
                                .entry(hap_name.clone())
                                .or_default()
                                .insert(String::from("sequence"), Value::String(new_hap.clone()));
                            if bp1 < bp2 {
                                let fusion_breakpoint = json!([[bp3, bp1], [bp4, bp2]]);
                                fusions_called
                                    .entry(hap_name.clone())
                                    .or_default()
                                    .insert(String::from("breakpoint"), fusion_breakpoint);
                            } else {
                                let fusion_breakpoint = json!([[bp4, bp2], [bp3, bp1]]);
                                fusions_called
                                    .entry(hap_name.clone())
                                    .or_default()
                                    .insert(String::from("breakpoint"), fusion_breakpoint);
                            }
                        }
                    }
                }
            }
        }
        Ok((assembled_haps_renamed, two_cp_haps, fusions_called))
    }

    /// Parse bundled fusion PSV definitions (`fusion_genes.json`).
    pub fn parse_psv(&mut self) -> Result<BTreeMap<String, Vec<CandidateSite>>, DError> {
        const PSV_FILE: &[u8] = std::include_bytes!(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/data/38/fusion_genes.json"
        ));
        let tmp = std::str::from_utf8(PSV_FILE)?;
        let psv_file_parsed: Value = serde_json::from_str(tmp)?;
        let mut fusion_gene_def_variants = BTreeMap::new();
        let psv_map = psv_file_parsed
            .as_object()
            .ok_or_else(|| phaser::Exception::new("fusion_genes.json root should be an object"))?;
        for (a, b) in psv_map {
            let gene_name = a.as_str().to_string();
            let variants = b.as_array().ok_or_else(|| {
                phaser::Exception::new(format!(
                    "fusion_genes.json['{}'] should be an array",
                    gene_name
                ))
            })?;
            let mut variant_sites = Vec::new();
            for var in variants {
                let var_str = var.as_str().ok_or_else(|| {
                    phaser::Exception::new(format!(
                        "fusion_genes.json['{}'] contains a non-string variant entry",
                        gene_name
                    ))
                })?;
                let fields = var_str.split_terminator('_').collect::<Vec<_>>();
                if fields.len() != 3 {
                    return Err(format!("Malformed fusion variant: {var_str}").into());
                }
                let var_pos = fields[0].parse::<i64>()? - 1;
                let ref_base = fields[1].to_string();
                let alt_base = fields[2].to_string();
                variant_sites.push(CandidateSite::new(var_pos, ref_base, alt_base));
            }
            fusion_gene_def_variants.insert(gene_name, variant_sites);
        }
        Ok(fusion_gene_def_variants)
    }

    /// Build a reduced haplotype sequence used to infer fusion breakpoints.
    ///
    /// Uses PSV sites when available; otherwise falls back to local het/hom sites.
    pub fn new_hap_for_breakpoint(
        &mut self,
        hap: VString,
        psv_variants: Option<&Vec<CandidateSite>>,
    ) -> Result<(String, Vec<CandidateSite>), DError> {
        let mut new_hap = String::from("");
        let mut all_sites: Vec<CandidateSite> = Vec::new();
        if psv_variants.is_none() {
            // no known psv sites, use existing variant sites
            all_sites.extend(self.het_sites.clone());
            all_sites.extend(self.hom_sites.clone());
            all_sites.sort_by(|a, b| a.pos.cmp(&b.pos));
            if !self.clip_5p_positions.is_empty() {
                if let Some(max_clip_5p) = self.clip_5p_positions.iter().max() {
                    all_sites = all_sites
                        .iter()
                        .filter(|x| x.pos > *max_clip_5p)
                        .cloned()
                        .collect::<Vec<_>>();
                }
            }
            if !self.clip_3p_positions.is_empty() {
                if let Some(min_clip_3p) = self.clip_3p_positions.iter().min() {
                    all_sites = all_sites
                        .iter()
                        .filter(|x| x.pos < *min_clip_3p)
                        .cloned()
                        .collect::<Vec<_>>();
                }
            }
            for var_site in &all_sites {
                if self.hom_sites.contains(var_site) {
                    new_hap.push('2');
                } else if self.het_sites.contains(var_site) {
                    if let Some(this_position) = self.het_sites.iter().position(|x| x == var_site) {
                        new_hap.push(hap[this_position] as char);
                    }
                }
            }
        } else {
            all_sites = psv_variants
                .ok_or_else(|| {
                    phaser::Exception::new(format!(
                        "Missing PSV variants for gene '{}' while inferring fusion breakpoint",
                        self.gene_name()
                    ))
                })?
                .clone();
            for var_site in &all_sites {
                let mut base = '1';
                if self.hom_sites.contains(var_site) {
                    base = '2'
                } else if self.het_sites.contains(var_site) {
                    if let Some(this_position) = self.het_sites.iter().position(|x| x == var_site) {
                        base = hap[this_position] as char;
                    }
                }
                new_hap.push(base);
            }
        }
        Ok((new_hap, all_sites))
    }

    /// Rename haplotypes into gene1/gene2/fusion classes and infer two-copy haps.
    pub fn update_twp_cp_in_fusion_cases<'a>(
        &mut self,
        assembled_haps: &BTreeMap<VStr<'a>, String>,
    ) -> Result<(BTreeMap<VStr<'a>, String>, Vec<String>), DError> {
        let mut assembled_haps_renamed = BTreeMap::new();
        let gene_name = self.gene_name();
        let mut unknown_ends = false;
        let mut two_cp_haps = Vec::new();
        let mut gene1s = Vec::new();
        let mut gene2s = Vec::new();
        let mut fusions = Vec::new();
        let mut counter_gene1 = 0;
        let mut counter_gene2 = 0;
        let mut counter_fusion = 0;
        let mut counter_unknown = 0;
        for hap in assembled_haps.keys() {
            let first_base = hap.first().ok_or_else(|| {
                phaser::Exception::new("Encountered empty haplotype while renaming fusion haps")
            })?;
            let last_base = hap.last().ok_or_else(|| {
                phaser::Exception::new("Encountered empty haplotype while renaming fusion haps")
            })?;
            let first_base0 = *first_base == b'0';
            let last_base0 = *last_base == b'0';
            if *first_base == b'x' || *last_base == b'x' {
                unknown_ends = true;
                counter_unknown += 1;
                assembled_haps_renamed
                    .insert(*hap, format!("{gene_name}_unknownhap{}", counter_unknown));
            }
            if !first_base0 && !last_base0 {
                gene1s.push(*hap);
                counter_gene1 += 1;
                assembled_haps_renamed
                    .insert(*hap, format!("{gene_name}_gene1hap{}", counter_gene1));
            } else if first_base0 && last_base0 {
                gene2s.push(*hap);
                counter_gene2 += 1;
                assembled_haps_renamed
                    .insert(*hap, format!("{gene_name}_gene2hap{}", counter_gene2));
            } else if (first_base0 && !last_base0) || (!first_base0 && last_base0) {
                fusions.push(*hap);
                counter_fusion += 1;
                assembled_haps_renamed
                    .insert(*hap, format!("{gene_name}_fusionhap{}", counter_fusion));
            }
        }
        if !unknown_ends {
            if fusions.is_empty() && assembled_haps_renamed.len() < 4 {
                if gene1s.len() == 1 {
                    if let Some(gene1_hap) = gene1s.first() {
                        if let Some(gene1_hap_name) = assembled_haps_renamed.get(gene1_hap) {
                            two_cp_haps.push(gene1_hap_name.to_string());
                        }
                    }
                }
                if gene2s.len() == 1 {
                    if let Some(gene2_hap) = gene2s.first() {
                        if let Some(gene2_hap_name) = assembled_haps_renamed.get(gene2_hap) {
                            two_cp_haps.push(gene2_hap_name.to_string());
                        }
                    }
                }
            } else if fusions.len() == 1 && assembled_haps_renamed.len() == 1 {
                // homozygous fusion
                if let Some(fusion_hap) = fusions.first() {
                    if let Some(fusion_hap_name) = assembled_haps_renamed.get(fusion_hap) {
                        two_cp_haps.push(fusion_hap_name.to_string());
                    }
                }
            }
        }
        Ok((assembled_haps_renamed, two_cp_haps))
    }
}

/// get fusion type: deletion or duplication
pub fn get_fusion_type(fusion_direction: &str, hap: VString) -> Result<Option<String>, DError> {
    let first_base = hap
        .first()
        .ok_or_else(|| phaser::Exception::new("Cannot infer fusion type from empty haplotype"))?;
    let last_base = hap
        .last()
        .ok_or_else(|| phaser::Exception::new("Cannot infer fusion type from empty haplotype"))?;
    let first_base0 = *first_base == b'0';
    let last_base0 = *last_base == b'0';
    if fusion_direction == "5p" {
        if !last_base0 && first_base0 {
            return Ok(Some(String::from("duplication")));
        }
        if last_base0 && !first_base0 {
            return Ok(Some(String::from("deletion")));
        }
    }
    if fusion_direction == "3p" {
        if !last_base0 && first_base0 {
            return Ok(Some(String::from("deletion")));
        }
        if last_base0 && !first_base0 {
            return Ok(Some(String::from("duplication")));
        }
    }
    Ok(None)
}

/// Infer the switch from gene1 sequence to gene2 sequence
pub fn get_fusion_breakpoint_index(hap: VString, new_hap: String) -> Result<Option<usize>, DError> {
    let first_base = hap.first().ok_or_else(|| {
        phaser::Exception::new("Cannot infer fusion breakpoint from empty source haplotype")
    })?;
    let last_base = hap.last().ok_or_else(|| {
        phaser::Exception::new("Cannot infer fusion breakpoint from empty source haplotype")
    })?;
    let first_base0 = *first_base == b'0';
    let last_base0 = *last_base == b'0';
    let hap_len = new_hap.len();
    // 2s to 1s
    if first_base0 && !last_base0 {
        let mut counts = Vec::new();
        for i in 0..hap_len {
            let first_seg = new_hap[..i].to_string();
            let second_seg = new_hap[i..].to_string();
            let n = first_seg
                .into_bytes()
                .iter()
                .filter(|x| **x == b'2')
                .count()
                + second_seg
                    .into_bytes()
                    .iter()
                    .filter(|x| **x == b'1')
                    .count();
            counts.push(n);
        }
        let Some(max_count) = counts.iter().max() else {
            return Ok(None);
        };
        let Some(bp_index) = counts.iter().position(|x| x == max_count) else {
            return Ok(None);
        };
        if bp_index == 0 || bp_index == counts.len() - 1 {
            return Ok(None);
        }
        return Ok(Some(bp_index));
    }
    // 1s to 2s
    if !first_base0 && last_base0 {
        let mut counts = Vec::new();
        for i in 0..hap_len {
            let first_seg = new_hap[..i].to_string();
            let second_seg = new_hap[i..].to_string();
            let n = first_seg
                .into_bytes()
                .iter()
                .filter(|x| **x == b'1')
                .count()
                + second_seg
                    .into_bytes()
                    .iter()
                    .filter(|x| **x == b'2')
                    .count();
            counts.push(n);
        }
        let Some(max_count) = counts.iter().max() else {
            return Ok(None);
        };
        let Some(bp_index) = counts.iter().position(|x| x == max_count) else {
            return Ok(None);
        };
        if bp_index == 0 || bp_index == counts.len() - 1 {
            return Ok(None);
        }
        return Ok(Some(bp_index));
    }
    Ok(None)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config;
    use crate::phaser;
    use crate::phaser::Phaser;
    use crate::toolkit::site_selection::CandidateSite;
    use crate::toolkit::util;
    use std::str::FromStr;

    fn build_test_phaser() -> Option<Phaser> {
        let outdir = tempfile::TempDir::new().expect("tempdir should build");
        let Some(genome_path) = std::env::var("HG38")
            .ok()
            .map(|x| x.trim_end_matches(".mmi").to_string())
        else {
            log::warn!(
                "Skipping fusion test because the HG38 environment variable is not configured."
            );
            return None;
        };
        let settings = phaser::Settings::new(
            "HG00733",
            (genome_path, util::test_file("bams/HG00733.smn1.bam")),
            outdir.path(),
            "smn1",
            &config::Region::try_load(None).expect("region config should load"),
            None,
            None,
            String::from("38"),
            None,
            0.03,
            false,
        );
        let gene_config = config::Gene::try_load(None).expect("gene config should load");
        Some(Phaser::new(settings, Some(gene_config), None, None).expect("phaser should build"))
    }

    #[test]
    fn update_two_cp_in_fusion_cases_matches_python_cases() {
        let Some(mut phaser) = build_test_phaser() else {
            return;
        };
        let haplotypes = BTreeMap::from([
            (VStr::from("12121212"), String::from("hap1")),
            (VStr::from("01212120"), String::from("hap2")),
            (VStr::from("21212121"), String::from("hap3")),
            (VStr::from("02121210"), String::from("hap4")),
        ]);
        let (renamed_haps, two_cp_haps) = phaser
            .update_twp_cp_in_fusion_cases(&haplotypes)
            .expect("fusion rename should succeed");
        assert_eq!(
            renamed_haps
                .iter()
                .map(|(k, v)| (k.to_string(), v.clone()))
                .collect::<BTreeMap<_, _>>(),
            BTreeMap::from([
                (String::from("12121212"), String::from("smn1_gene1hap1")),
                (String::from("01212120"), String::from("smn1_gene2hap1")),
                (String::from("21212121"), String::from("smn1_gene1hap2")),
                (String::from("02121210"), String::from("smn1_gene2hap2")),
            ])
        );
        assert!(two_cp_haps.is_empty());

        let haplotypes = BTreeMap::from([
            (VStr::from("12121212"), String::from("hap1")),
            (VStr::from("01212120"), String::from("hap2")),
            (VStr::from("21212121"), String::from("hap3")),
        ]);
        let (_renamed_haps, two_cp_haps) = phaser
            .update_twp_cp_in_fusion_cases(&haplotypes)
            .expect("fusion rename should succeed");
        assert_eq!(two_cp_haps, vec![String::from("smn1_gene2hap1")]);

        let haplotypes = BTreeMap::from([
            (VStr::from("01212120"), String::from("hap1")),
            (VStr::from("02121210"), String::from("hap2")),
            (VStr::from("21212121"), String::from("hap3")),
        ]);
        let (_renamed_haps, two_cp_haps) = phaser
            .update_twp_cp_in_fusion_cases(&haplotypes)
            .expect("fusion rename should succeed");
        assert_eq!(two_cp_haps, vec![String::from("smn1_gene1hap1")]);

        let haplotypes = BTreeMap::from([
            (VStr::from("0121212x"), String::from("hap1")),
            (VStr::from("21212121"), String::from("hap2")),
            (VStr::from("02121210"), String::from("hap3")),
        ]);
        let (_renamed_haps, two_cp_haps) = phaser
            .update_twp_cp_in_fusion_cases(&haplotypes)
            .expect("fusion rename should succeed");
        assert!(two_cp_haps.is_empty());
    }

    #[test]
    fn get_fusion_type_matches_python_cases() {
        assert_eq!(
            get_fusion_type("5p", VString::from("012121")).unwrap(),
            Some(String::from("duplication"))
        );
        assert_eq!(
            get_fusion_type("5p", VString::from("121210")).unwrap(),
            Some(String::from("deletion"))
        );
        assert_eq!(
            get_fusion_type("5p", VString::from("121211")).unwrap(),
            None
        );
    }

    #[test]
    fn get_fusion_breakpoint_index_matches_python_cases() {
        assert_eq!(
            get_fusion_breakpoint_index(VString::from("121210"), String::from("111111122222"))
                .unwrap(),
            Some(7)
        );
        assert_eq!(
            get_fusion_breakpoint_index(VString::from("121210"), String::from("2222211111111"))
                .unwrap(),
            None
        );
        assert_eq!(
            get_fusion_breakpoint_index(VString::from("012121"), String::from("2222211111111"))
                .unwrap(),
            Some(5)
        );
        assert_eq!(
            get_fusion_breakpoint_index(VString::from("012121"), String::from("111111122222"))
                .unwrap(),
            None
        );
        assert_eq!(
            get_fusion_breakpoint_index(VString::from("112121"), String::from("111111122222"))
                .unwrap(),
            None
        );
    }

    #[test]
    fn new_hap_for_breakpoint_matches_python_cases() {
        let Some(mut phaser) = build_test_phaser() else {
            return;
        };
        phaser.hom_sites = vec![CandidateSite::from_str("7_C_T").unwrap()];
        phaser.het_sites = vec![
            CandidateSite::from_str("1_A_T").unwrap(),
            CandidateSite::from_str("3_C_T").unwrap(),
            CandidateSite::from_str("11_C_T").unwrap(),
        ];
        let hap = VString::from("212");

        let psv_variants = vec![
            CandidateSite::from_str("1_A_T").unwrap(),
            CandidateSite::from_str("3_C_T").unwrap(),
            CandidateSite::from_str("5_A_T").unwrap(),
            CandidateSite::from_str("7_C_T").unwrap(),
            CandidateSite::from_str("9_A_T").unwrap(),
            CandidateSite::from_str("11_C_T").unwrap(),
        ];
        let (new_hap, all_sites) = phaser
            .new_hap_for_breakpoint(hap.clone(), Some(&psv_variants))
            .unwrap();
        assert_eq!(new_hap, String::from("211212"));
        assert_eq!(all_sites, psv_variants);

        let (new_hap, all_sites) = phaser.new_hap_for_breakpoint(hap.clone(), None).unwrap();
        assert_eq!(new_hap, String::from("2122"));
        assert_eq!(
            all_sites
                .iter()
                .map(std::string::ToString::to_string)
                .collect::<Vec<_>>(),
            vec![
                String::from("1_A_T"),
                String::from("3_C_T"),
                String::from("7_C_T"),
                String::from("11_C_T")
            ]
        );

        phaser.clip_3p_positions = vec![10];
        let (new_hap, all_sites) = phaser.new_hap_for_breakpoint(hap.clone(), None).unwrap();
        assert_eq!(new_hap, String::from("212"));
        assert_eq!(
            all_sites
                .iter()
                .map(std::string::ToString::to_string)
                .collect::<Vec<_>>(),
            vec![
                String::from("1_A_T"),
                String::from("3_C_T"),
                String::from("7_C_T")
            ]
        );

        phaser.clip_5p_positions = vec![1];
        let (new_hap, all_sites) = phaser.new_hap_for_breakpoint(hap, None).unwrap();
        assert_eq!(new_hap, String::from("12"));
        assert_eq!(
            all_sites
                .iter()
                .map(std::string::ToString::to_string)
                .collect::<Vec<_>>(),
            vec![String::from("3_C_T"), String::from("7_C_T")]
        );
    }
}
