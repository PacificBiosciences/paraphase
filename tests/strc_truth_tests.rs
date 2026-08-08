use paraphase::toolkit::util::DResult;
use paraphase::{
    config, depth,
    phaser::{self, Phaser},
    toolkit::util::{self, test_file},
};

use std::collections::BTreeSet;

fn strc_haplotype_class(name: &str) -> &'static str {
    let identity = name
        .rsplit('_')
        .next()
        .expect("STRC haplotype name should contain an identity label");
    let has_numeric_suffix = |prefix: &str| {
        identity.strip_prefix(prefix).is_some_and(|suffix| {
            !suffix.is_empty() && suffix.chars().all(|digit| digit.is_ascii_digit())
        })
    };

    if has_numeric_suffix("strchap") {
        "strc"
    } else if has_numeric_suffix("strcp1hap") {
        "strcp1"
    } else {
        panic!("Unexpected STRC haplotype identity label in {name}");
    }
}

// Public NIST GIAB HG005 HiFi locus fixture (524 reads), with genome depth 34
// injected. The identity/CN truth target is described by Tsai et al. (2024).
#[test]
fn hg005_strc_truth_ok() -> DResult {
    util::init_log(log::LevelFilter::Info);

    let outdir = tempfile::TempDir::new()?;
    let genome_bam = test_file("bams/HG005.strc.bam");
    let genome_path = if let Ok(x) = std::env::var("HG38") {
        x.trim_end_matches(".mmi").to_string()
    } else {
        log::warn!("Skipping test: HG38 env is not set.");
        return Ok(());
    };
    let depth = depth::Result {
        median: 34.0,
        median_absolute_difference: 0.1,
        sex: depth::Sex::Male,
    };
    let region_config = config::Region::try_load(None)?;
    let settings = phaser::Settings::new(
        "HG005",
        (genome_path, genome_bam),
        outdir.path(),
        "strc",
        &region_config,
        Some(depth),
        None,
        String::from("38"),
        None,
        0.03,
        false,
    );
    let mut phaser = Phaser::new(settings, Some(config::Gene::try_load(None)?), None, None)?;
    let call = phaser.run()?;

    assert!(!call.failed_for_coverage);
    assert_eq!(phaser.region_avg_depth[0], (81.0, 85.0));
    assert_eq!(call.total_cn, Some(4));
    assert_eq!(
        call.region_specific_info
            .get("gene_cn")
            .and_then(serde_json::Value::as_i64),
        Some(3)
    );
    assert_eq!(
        call.region_specific_info
            .get("intergenic_depth")
            .and_then(serde_json::Value::as_f64),
        Some(34.0)
    );

    let two_copy = call
        .two_copy_haplotypes
        .iter()
        .cloned()
        .collect::<BTreeSet<_>>();
    let copies = |name: &str| 1 + usize::from(two_copy.contains(name));
    let mut copy_classes = Vec::new();
    for name in call.final_haplotypes.values() {
        for _ in 0..copies(name) {
            copy_classes.push(strc_haplotype_class(name));
        }
    }
    copy_classes.sort_unstable();
    assert_eq!(copy_classes, vec!["strc", "strc", "strc", "strcp1"]);

    let marker_index = |site: &str| {
        let indices = call
            .sites_for_phasing
            .iter()
            .enumerate()
            .filter_map(|(index, found_site)| (found_site == site).then_some(index))
            .collect::<Vec<_>>();
        assert_eq!(
            indices.len(),
            1,
            "Expected exactly one {site} site for STRC identity markers"
        );
        indices[0]
    };
    let strc_snv_index = marker_index("43602487_C_G");
    let strc_deletion_index = marker_index("43602630_del_314");
    let mut identity_marker_pairs = BTreeSet::new();
    for (sequence, name) in &call.final_haplotypes {
        let marker = |index: usize| {
            sequence
                .as_bytes()
                .get(index)
                .copied()
                .map(char::from)
                .unwrap_or_else(|| panic!("Missing STRC identity marker in haplotype {sequence}"))
        };
        let found_marker_pair = (marker(strc_snv_index), marker(strc_deletion_index));
        let expected_marker_pair = match strc_haplotype_class(name) {
            "strc" => ('1', '1'),
            "strcp1" => ('2', '3'),
            identity => unreachable!("Unhandled STRC haplotype identity {identity}"),
        };
        assert_eq!(
            found_marker_pair, expected_marker_pair,
            "Unexpected identity markers for {name} on haplotype {sequence}"
        );
        identity_marker_pairs.insert(found_marker_pair);
    }
    assert_eq!(
        identity_marker_pairs,
        BTreeSet::from([('1', '1'), ('2', '3')])
    );

    let c5125_strc_copy_count = call
        .final_haplotypes
        .values()
        .filter(|name| strc_haplotype_class(name) == "strc")
        .filter(|name| {
            call.haplotype_details.get(*name).is_some_and(|detail| {
                detail
                    .variants
                    .iter()
                    .any(|variant| variant == "43600074_T_C")
            })
        })
        .map(|name| copies(name))
        .sum::<usize>();
    assert_eq!(c5125_strc_copy_count, 1);

    Ok(())
}
