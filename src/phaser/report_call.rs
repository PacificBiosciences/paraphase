//! Fill final report fields on [`GeneCall`] from phasing outputs.

use crate::io::json::GeneCall;
use crate::phaser::Phaser;
use crate::phaser::{PhasedResult, PhasedResultForJson};

use std::collections::BTreeMap;

impl Phaser {
    /// Build a baseline `GeneCall` populated with run metadata and region depth.
    pub fn get_default_call(&self) -> GeneCall {
        let (region_median_depth, region_80percentile_depth) = self
            .region_avg_depth
            .first()
            .copied()
            .unwrap_or((f32::NAN, f32::NAN));
        let mut region_depth = BTreeMap::new();
        region_depth.insert(String::from("median"), region_median_depth);
        region_depth.insert(String::from("percentile80"), region_80percentile_depth);
        let genes_in_region = self
            .locus_config()
            .get("genes")
            .and_then(|x| x.as_str())
            .map(ToString::to_string);
        GeneCall {
            gene_name: self.settings.gene_name.clone(),
            sample_sex: match self.settings.sample_sex {
                crate::depth::Sex::Other => None,
                sex => Some(sex.to_string()),
            },
            genome_depth: self.settings.depth.map(|x| x.median as f32),
            region_depth,
            genes_in_region,
            phase_region: format!(
                "{}:{}:{}-{}",
                self.settings.genome,
                self.chr().unwrap_or("unknown"),
                self.left_boundary(),
                self.right_boundary()
            ),
            ..Default::default()
        }
    }

    /// Populate call-site arrays from the current phaser state.
    ///
    /// This is shared by both normal reporting and early-return failure paths.
    pub fn populate_call_site_fields(&mut self, call: &mut GeneCall) {
        call.sites_for_phasing = self
            .het_sites
            .iter()
            .map(ToString::to_string)
            .collect::<Vec<_>>();
        self.init_het_sites.sort_by_key(|a| a.pos);
        call.heterozygous_sites = self
            .init_het_sites
            .iter()
            .map(ToString::to_string)
            .collect::<Vec<_>>();
        call.homozygous_sites = self
            .hom_sites
            .iter()
            .map(ToString::to_string)
            .collect::<Vec<_>>();
        self.het_sites_no_phasing.sort_by_key(|a| a.pos);
        call.het_sites_not_used_in_phasing = self
            .het_sites_no_phasing
            .iter()
            .map(ToString::to_string)
            .collect::<Vec<_>>();
    }

    /// Fill in fields for the final report.
    pub fn fill_in_call(&mut self, phase_results: PhasedResult, call: &mut GeneCall) {
        let phase_result = PhasedResultForJson::new(&phase_results);
        self.populate_call_site_fields(call);
        call.assembled_haplotypes = phase_result
            .assemblies
            .final_haps
            .iter()
            .map(std::string::ToString::to_string)
            .collect::<Vec<_>>();
        call.highest_total_cn = Some(phase_result.assemblies.highest_cn as i32);
        call.unique_supporting_reads = phase_result
            .uniquely_supporting_reads
            .iter()
            .map(|(k, v)| (k.to_string(), v.clone()))
            .collect::<BTreeMap<_, _>>();
        call.nonunique_supporting_reads = phase_result
            .nonuniquely_supporting_reads
            .iter()
            .map(|(k, v)| {
                (
                    k.clone(),
                    v.iter()
                        .map(std::borrow::ToOwned::to_owned)
                        .collect::<Vec<_>>(),
                )
            })
            .collect::<BTreeMap<_, _>>();
    }
}
