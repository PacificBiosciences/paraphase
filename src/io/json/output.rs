use crate::phaser::HapInfoForJson;
use crate::toolkit::util::DError;

use std::collections::BTreeMap;

///
/// Primary per-gene output payload emitted by the phasing pipeline.
///
/// Most fields mirror Python output keys directly for parity. Additional,
/// gene-specific metadata can be attached via `region_specific_info`.
#[must_use]
#[derive(Clone, Debug, Default, serde::Serialize)]
pub struct GeneCall {
    #[serde(rename = "region_name")]
    pub gene_name: String,
    pub phase_region: String,
    pub genes_in_region: Option<String>,
    pub sample_sex: Option<String>,
    pub genome_depth: Option<f32>,
    pub region_depth: BTreeMap<String, f32>,
    pub failed_for_coverage: bool,

    pub total_cn: Option<i32>,
    pub final_haplotypes: BTreeMap<String, String>,
    pub two_copy_haplotypes: Vec<String>,
    pub region_specific_info: BTreeMap<String, serde_json::Value>, // additional key-value metadata

    pub sites_for_phasing: Vec<String>,
    pub assembled_haplotypes: Vec<String>,
    pub unique_supporting_reads: BTreeMap<String, Vec<String>>,

    pub highest_total_cn: Option<i32>,
    pub heterozygous_sites: Vec<String>,
    pub het_sites_not_used_in_phasing: Vec<String>,
    pub homozygous_sites: Vec<String>,
    pub haplotype_details: BTreeMap<String, HapInfoForJson>,
    pub nonunique_supporting_reads: BTreeMap<String, Vec<String>>,
    pub read_details: BTreeMap<String, String>,
}

/// Serialize per-gene calls as pretty JSON and append a trailing newline.
///
/// This is the canonical JSON output writer used by the pipeline entrypoint.
pub fn write_outputs(
    x: &BTreeMap<String, GeneCall>,
    writer: &mut impl std::io::Write,
) -> Result<(), DError> {
    writeln!(writer, "{}", serde_json::to_string_pretty(x)?)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn write_outputs_emits_json() -> Result<(), DError> {
        let mut calls = BTreeMap::new();
        calls.insert(
            String::from("GENE1"),
            GeneCall {
                gene_name: String::from("GENE1"),
                ..GeneCall::default()
            },
        );
        let mut buf = Vec::<u8>::new();
        write_outputs(&calls, &mut buf)?;
        let rendered = String::from_utf8(buf).expect("json writer should emit utf8");
        assert!(rendered.contains("\"GENE1\""), "output was: {rendered}");
        assert!(rendered.ends_with('\n'), "output should end with newline");
        let parsed: serde_json::Value =
            serde_json::from_str(rendered.trim_end()).expect("must parse as json");
        assert!(parsed.get("GENE1").is_some(), "missing GENE1 key");
        Ok(())
    }
}
