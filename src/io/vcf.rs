//! VCF I/O module.
//!
//! This top-level module follows the same layout as `io::json`:
//! a facade file (`vcf.rs`) with focused submodules under `io/vcf/`.
//!
//! Structure:
//! - [`writer`]: `VcfWriter` orchestration (header building, variant collection, record writing).
//! - [`helpers`]: pure/helper logic and VCF constants used by the writer.
//! - [`types`]: shared data structures passed between collection and merge steps.

mod helpers;
mod types;
mod writer;

pub use self::types::*;
pub use self::writer::*;
#[cfg(test)]
mod tests {
    use super::helpers::*;
    use super::*;
    use crate::config;
    use crate::io::json::GeneCall;
    use crate::phaser;
    use crate::phaser::Phaser;
    use crate::toolkit::util;
    use rust_htslib::bcf::record::GenotypeAllele;
    use serde_json::json;
    use std::collections::BTreeMap;

    fn hg38_reference() -> Option<std::path::PathBuf> {
        std::env::var("HG38").ok().map(std::path::PathBuf::from)
    }

    fn build_test_phaser(gene: &str, outdir: &std::path::Path) -> Option<Phaser> {
        let Some(hg38_reference) = hg38_reference() else {
            log::warn!(
                "Skipping VCF test because the HG38 environment variable is not configured."
            );
            return None;
        };
        let settings = phaser::Settings::new(
            "TEST",
            (hg38_reference, util::test_file("bams/HG00733.smn1.bam")),
            outdir,
            gene,
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

    fn stage_vcf_fixtures(outdir: &std::path::Path, sample_id: &str, gene: &str) {
        std::fs::create_dir_all(outdir).expect("test outdir should be creatable");
        let dest_bam = outdir.join(format!("{sample_id}_{gene}_realigned_tagged.bam"));
        let dest_bai = outdir.join(format!("{sample_id}_{gene}_realigned_tagged.bam.bai"));
        let src_bam = util::test_file(&format!("bams/HG00733.{gene}.bam"));
        let src_bai = util::test_file(&format!("bams/HG00733.{gene}.bam.bai"));
        if src_bam.exists() && src_bai.exists() {
            if !dest_bam.exists() {
                std::fs::copy(&src_bam, &dest_bam).expect("tagged BAM fixture should copy");
            }
            if !dest_bai.exists() {
                std::fs::copy(&src_bai, &dest_bai).expect("tagged BAM index should copy");
            }
        }
        let src_local_ref = util::test_file(&format!("ref/{gene}_ref.fa"));
        if src_local_ref.exists() {
            let dest_local_ref = outdir.join(format!("{gene}_ref.fa"));
            if !dest_local_ref.exists() {
                std::fs::copy(&src_local_ref, &dest_local_ref)
                    .expect("local reference fixture should copy");
            }
            crate::phaser::build_faidx(&dest_local_ref)
                .expect("local reference index should build");
        }
    }

    fn write_empty_tagged_bam_for_region(phaser: &Phaser) {
        use rust_htslib::bam;
        let path = phaser.realigned_tagged_bam_path();
        if path.exists() {
            return;
        }
        let mut header = bam::Header::new();
        let chr = phaser
            .chr()
            .expect("gene should have a chromosome for empty tagged BAM header");
        let ln = phaser.right_boundary_0based() + 1;
        let mut sq = bam::header::HeaderRecord::new(b"SQ");
        sq.push_tag(b"SN", chr);
        sq.push_tag(b"LN", ln);
        header.push_record(&sq);
        let writer = bam::Writer::from_path(&path, &header, bam::Format::Bam)
            .expect("empty tagged BAM should be creatable");
        drop(writer);
        bam::index::build(&path, None, bam::index::Type::Bai, 1)
            .expect("empty tagged BAM index should build");
    }

    #[test]
    fn convert_alt_record_matches_python_cases() {
        assert_eq!(
            convert_alt_record(String::from("T"), String::from("A")),
            "A"
        );
        assert_eq!(
            convert_alt_record(String::from("T"), String::from("TACG")),
            "T+3ACG"
        );
        assert_eq!(
            convert_alt_record(String::from("TGC"), String::from("T")),
            "T-2GC"
        );
    }

    #[test]
    fn format_hap_bound_matches_python_cases() {
        assert_eq!(format_hap_bound(1, 2, &[]), "1-2");
        assert_eq!(
            format_hap_bound(1, 2, &[String::from("5p")]),
            "1truncated-2"
        );
        assert_eq!(
            format_hap_bound(1, 2, &[String::from("3p")]),
            "1-2truncated"
        );
        assert_eq!(
            format_hap_bound(1, 2, &[String::from("5p"), String::from("3p")]),
            "1truncated-2truncated"
        );
    }

    #[test]
    fn consensus_and_gt_match_key_python_get_var_cases() {
        let ref_seq = b"AAAAAAA";

        let call = get_consensus_var(vec![b"A".to_vec(); 10], 0, ref_seq).expect("consensus");
        assert_eq!(call.base.as_deref(), Some("A"));
        assert_eq!(call.depth, 10);
        assert_eq!(call.nread, 10);
        assert_eq!(
            get_gt(&call, &String::from("A"), &String::from("T")),
            GenotypeAllele::Phased(0)
        );

        let call = get_consensus_var(vec![], 0, ref_seq).expect("empty consensus");
        assert_eq!(call.base, None);
        assert_eq!(call.depth, 0);
        assert_eq!(call.nread, 0);
        assert_eq!(
            get_gt(&call, &String::from("A"), &String::from("T")),
            GenotypeAllele::UnphasedMissing
        );

        let call = get_consensus_var(vec![b"T".to_vec(); 4], 0, ref_seq).expect("alt consensus");
        assert_eq!(call.base.as_deref(), Some("T"));
        assert_eq!(call.depth, 4);
        assert_eq!(call.nread, 4);
        assert_eq!(
            get_gt(&call, &String::from("A"), &String::from("T")),
            GenotypeAllele::Phased(1)
        );

        let call = get_consensus_var(vec![b"*".to_vec(); 4], 0, ref_seq).expect("star consensus");
        assert_eq!(call.base.as_deref(), Some("*"));
        assert_eq!(call.depth, 4);
        assert_eq!(call.nread, 4);
        assert_eq!(
            get_gt(&call, &String::from("A"), &String::from("*")),
            GenotypeAllele::Phased(1)
        );

        let call = get_consensus_var(
            vec![
                b"A+2GT".to_vec(),
                b"A+2GT".to_vec(),
                b"A+2GT".to_vec(),
                b"A+2GT".to_vec(),
                b"A".to_vec(),
            ],
            0,
            ref_seq,
        )
        .expect("ins consensus");
        assert_eq!(call.base.as_deref(), Some("AGT"));
        assert_eq!(call.ref_base, "A");
        assert_eq!(call.original_base.as_deref(), Some("A+2GT"));
        assert_eq!(call.depth, 5);
        assert_eq!(call.nread, 4);
        assert_eq!(
            get_gt(&call, &String::from("A"), &String::from("AGT")),
            GenotypeAllele::Phased(1)
        );

        let call = get_consensus_var(
            vec![
                b"A".to_vec(),
                b"A".to_vec(),
                b"A".to_vec(),
                b"T".to_vec(),
                b"T".to_vec(),
            ],
            0,
            ref_seq,
        )
        .expect("mixed consensus");
        assert_eq!(call.base.as_deref(), Some("A"));
        assert_eq!(call.depth, 5);
        assert_eq!(call.nread, 3);
        assert_eq!(
            get_gt(&call, &String::from("A"), &String::from("T")),
            GenotypeAllele::UnphasedMissing
        );
    }

    #[test]
    fn low_depth_alt_consensus_stays_missing_in_rust() {
        let ref_seq = b"AAAAAAA";
        let call = get_consensus_var(vec![b"T".to_vec(); 2], 0, ref_seq).expect("alt consensus");
        assert_eq!(call.base.as_deref(), Some("T"));
        assert_eq!(call.depth, 2);
        assert_eq!(call.nread, 2);
        assert_eq!(
            get_gt(&call, &String::from("A"), &String::from("T")),
            GenotypeAllele::UnphasedMissing
        );
    }

    #[test]
    fn python_low_depth_get_var_examples_are_intentionally_missing_in_rust() {
        let ref_seq = b"AAAAAAA";

        // Python test_get_var expects GT=1 for this low-depth SNV case.
        // Rust intentionally emits missing due to MIN_DEPTH gating.
        let snv = get_consensus_var(vec![b"T".to_vec(); 2], 0, ref_seq).expect("snv consensus");
        assert_eq!(snv.base.as_deref(), Some("T"));
        assert_eq!(snv.depth, 2);
        assert_eq!(snv.nread, 2);
        assert_eq!(
            get_gt(&snv, &String::from("A"), &String::from("T")),
            GenotypeAllele::UnphasedMissing
        );

        // Python test_get_var expects GT=1 for this low-depth deletion case.
        // Rust intentionally emits missing for the same reason.
        let del =
            get_consensus_var(vec![b"*".to_vec(); 3], 0, ref_seq).expect("deletion consensus");
        assert_eq!(del.base.as_deref(), Some("*"));
        assert_eq!(del.depth, 3);
        assert_eq!(del.nread, 3);
        assert_eq!(
            get_gt(&del, &String::from("A"), &String::from("*")),
            GenotypeAllele::UnphasedMissing
        );
    }

    #[test]
    fn homozygous_case_hap_info_matches_python_shape() {
        let outdir = tempfile::TempDir::new().expect("tempdir should build");
        stage_vcf_fixtures(outdir.path(), "TEST", "smn1");
        let Some(phaser) = build_test_phaser("smn1", outdir.path()) else {
            return;
        };
        let call = GeneCall {
            heterozygous_sites: vec![],
            final_haplotypes: BTreeMap::new(),
            ..Default::default()
        };
        let writer = VcfWriter::new(&phaser, &call, false, false);

        let hap_variant_info = writer
            .get_variants_for_vcf(&call.final_haplotypes, false, false)
            .expect("homozygous no-realign VCF path should succeed");

        assert_eq!(hap_variant_info.hap_info.len(), 2);
        assert_eq!(
            hap_variant_info.hap_info[0].hap_name,
            String::from("smn1_homozygous_hap1")
        );
        assert_eq!(
            hap_variant_info.hap_info[1].hap_name,
            String::from("smn1_homozygous_hap1_cp2")
        );
        assert!(hap_variant_info.hap_info[0].is_truncated.is_empty());
        assert!(hap_variant_info.hap_info[1].is_truncated.is_empty());
    }

    #[test]
    fn truncated_two_copy_hap_info_matches_python_shape() {
        let outdir = tempfile::TempDir::new().expect("tempdir should build");
        stage_vcf_fixtures(outdir.path(), "TEST", "smn1");
        let Some(phaser) = build_test_phaser("smn1", outdir.path()) else {
            return;
        };
        let left = phaser.left_boundary();
        let right = phaser.right_boundary();

        let call = GeneCall {
            final_haplotypes: BTreeMap::from([
                (String::from("111"), String::from("smn1_hap1")),
                (String::from("222"), String::from("smn1_hap2")),
                (String::from("333"), String::from("smn1_hap3")),
            ]),
            two_copy_haplotypes: vec![String::from("smn1_hap1")],
            haplotype_details: BTreeMap::from([
                (
                    String::from("smn1_hap1"),
                    phaser::HapInfoForJson {
                        variants: vec![],
                        boundary: [left, right],
                        boundary_gene2: None,
                        is_truncated: vec![],
                    },
                ),
                (
                    String::from("smn1_hap2"),
                    phaser::HapInfoForJson {
                        variants: vec![],
                        boundary: [left + 100, right],
                        boundary_gene2: None,
                        is_truncated: vec![String::from("5p")],
                    },
                ),
                (
                    String::from("smn1_hap3"),
                    phaser::HapInfoForJson {
                        variants: vec![],
                        boundary: [left + 100, right],
                        boundary_gene2: None,
                        is_truncated: vec![String::from("5p")],
                    },
                ),
            ]),
            ..Default::default()
        };
        let writer = VcfWriter::new(&phaser, &call, false, false);

        let hap_variant_info = writer
            .get_variants_for_vcf(&call.final_haplotypes, false, false)
            .expect("two-copy/truncated no-realign VCF path should succeed");

        assert_eq!(hap_variant_info.hap_info.len(), 4);
        assert_eq!(
            hap_variant_info
                .hap_info
                .iter()
                .map(|x| x.hap_name.clone())
                .collect::<Vec<_>>(),
            vec![
                String::from("smn1_hap1"),
                String::from("smn1_hap1_cp2"),
                String::from("smn1_hap2"),
                String::from("smn1_hap3"),
            ]
        );
        assert_eq!(hap_variant_info.hap_info[0].start, left - 1);
        assert_eq!(hap_variant_info.hap_info[0].end, right - 1);
        assert!(hap_variant_info.hap_info[0].is_truncated.is_empty());
        assert_eq!(hap_variant_info.hap_info[1].start, left - 1);
        assert_eq!(hap_variant_info.hap_info[1].end, right - 1);
        assert!(hap_variant_info.hap_info[1].is_truncated.is_empty());
        assert_eq!(hap_variant_info.hap_info[2].start, left + 99);
        assert_eq!(hap_variant_info.hap_info[2].end, right - 1);
        assert_eq!(
            hap_variant_info.hap_info[2].is_truncated,
            vec![String::from("5p")]
        );
        assert_eq!(hap_variant_info.hap_info[3].start, left + 99);
        assert_eq!(hap_variant_info.hap_info[3].end, right - 1);
        assert_eq!(
            hap_variant_info.hap_info[3].is_truncated,
            vec![String::from("5p")]
        );
    }

    #[test]
    fn ikbkg_symbolic_deletion_is_present_in_hap_variant_info() {
        let outdir = tempfile::TempDir::new().expect("tempdir should build");
        stage_vcf_fixtures(outdir.path(), "TEST", "ikbkg");
        let Some(phaser) = build_test_phaser("ikbkg", outdir.path()) else {
            return;
        };
        write_empty_tagged_bam_for_region(&phaser);
        let call = GeneCall {
            final_haplotypes: BTreeMap::from([
                (String::from("111"), String::from("ikbkg_hap1")),
                (String::from("222"), String::from("ikbkg_hap2")),
            ]),
            haplotype_details: BTreeMap::from([
                (
                    String::from("ikbkg_hap1"),
                    phaser::HapInfoForJson {
                        variants: vec![],
                        boundary: [phaser.left_boundary(), phaser.right_boundary()],
                        boundary_gene2: None,
                        is_truncated: vec![],
                    },
                ),
                (
                    String::from("ikbkg_hap2"),
                    phaser::HapInfoForJson {
                        variants: vec![],
                        boundary: [phaser.left_boundary(), phaser.right_boundary()],
                        boundary_gene2: None,
                        is_truncated: vec![],
                    },
                ),
            ]),
            region_specific_info: BTreeMap::from([(
                String::from("deletion_haplotypes"),
                json!(["ikbkg_hap1"]),
            )]),
            ..Default::default()
        };
        let writer = VcfWriter::new(&phaser, &call, false, false);

        let hap_variant_info = writer
            .get_variants_for_vcf(&call.final_haplotypes, false, false)
            .expect("ikbkg no-realign VCF path should succeed");

        let deletion_name = phaser
            .locus_config()
            .get("deletion1_in_gene1")
            .and_then(|value| value.as_str())
            .expect("ikbkg config should include deletion1_in_gene1")
            .to_string();
        let (start_1based, _, _) =
            parse_symbolic_sv(&deletion_name).expect("deletion name should parse");
        let start_0based = start_1based - 1;
        let pos_calls = hap_variant_info
            .sv_variants
            .get(&start_0based)
            .expect("symbolic deletion should be present at the deletion anchor");
        assert_eq!(pos_calls.len(), hap_variant_info.hap_info.len());
        assert_eq!(pos_calls[0], Some(deletion_name));
        assert_eq!(pos_calls[1], None);
    }

    #[test]
    fn special_variants_include_ikbkg_deletions() {
        let outdir = tempfile::TempDir::new().expect("tempdir should build");
        let Some(phaser) = build_test_phaser("ikbkg", outdir.path()) else {
            return;
        };
        let deletion_name_gene1 = phaser
            .locus_config()
            .get("deletion1_in_gene1")
            .and_then(|value| value.as_str())
            .expect("ikbkg config should include deletion1_in_gene1")
            .to_string();
        let deletion_name_gene2 = phaser
            .locus_config()
            .get("deletion1_in_gene2")
            .and_then(|value| value.as_str())
            .expect("ikbkg config should include deletion1_in_gene2")
            .to_string();
        let mut call = GeneCall::default();
        call.region_specific_info.insert(
            String::from("deletion_haplotypes"),
            json!(["ikbkg_hap1", "ikbkg_pseudohap1"]),
        );
        let writer = VcfWriter::new(&phaser, &call, false, false);

        let gene1 = writer
            .get_special_variants(false)
            .expect("special variants should be collected");
        assert_eq!(gene1.get("ikbkg_hap1"), Some(&deletion_name_gene1));
        assert_eq!(gene1.get("ikbkg_pseudohap1"), Some(&deletion_name_gene1));

        let gene2 = writer
            .get_special_variants(true)
            .expect("special variants should be collected");
        assert_eq!(gene2.get("ikbkg_hap1"), Some(&deletion_name_gene1));
        assert_eq!(gene2.get("ikbkg_pseudohap1"), Some(&deletion_name_gene2));
    }

    #[test]
    fn special_variants_include_f8_symbolic_events() {
        let outdir = tempfile::TempDir::new().expect("tempdir should build");
        let Some(phaser) = build_test_phaser("f8", outdir.path()) else {
            return;
        };
        let mut call = GeneCall::default();
        call.region_specific_info.insert(
            String::from("sv_called"),
            json!({
                "f8_int22invhap1": "inversion",
                "f8_int22delhap1": "deletion"
            }),
        );
        let writer = VcfWriter::new(&phaser, &call, false, false);
        let variants = writer
            .get_special_variants(false)
            .expect("special variants should be collected");

        assert_eq!(
            variants.get("f8_int22invhap1"),
            Some(&String::from("154890327_INV_155454650"))
        );
        assert_eq!(
            variants.get("f8_int22delhap1"),
            Some(&String::from("154890327_DEL_155376007"))
        );
    }
}
