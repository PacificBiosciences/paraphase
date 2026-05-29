use crate::depth::Sex;
use crate::io::json::{ReadAlignmentId, ReadFingerprintMap};
use crate::phaser;
use crate::phaser::PhasedResult;
use crate::toolkit::util::DError;
use itertools::{iproduct, Itertools};
use std::collections::{BTreeMap, BTreeSet, HashSet};
use vstr::{VStr, VString};

#[derive(Debug, Clone)]
pub struct AlleleResult {
    pub alleles: Vec<Vec<String>>,
    pub raw_alleles: Vec<Vec<String>>,
    pub haplotype_links: BTreeMap<String, Vec<String>>,
    pub directed_links: BTreeMap<(String, String), usize>,
    pub directed_links_loose: BTreeMap<String, BTreeMap<String, i32>>,
}

/// Links between haplotypes.
/// Used for linking haplotypes in tandem duplications.
#[derive(Debug, Clone)]
pub struct DirectedLinks<'a>(
    pub BTreeMap<(&'a str, &'a str), Vec<i32>>,
    pub BTreeMap<(&'a str, &'a str), Vec<i32>>,
    pub BTreeMap<&'a str, Vec<&'a str>>,
);

impl phaser::Phaser {
    /// Check if we should connect haplotypes based on reads that link them together.
    fn check_linking_read(aln1: VStr<'_>, aln2: VStr<'_>, is_reverse: bool) -> Option<[u8; 2]> {
        // Check first/last bases for 'x'
        let (Some(start1), Some(start2), Some(end1), Some(end2)) = (
            aln1.first().copied(),
            aln2.first().copied(),
            aln1.last().copied(),
            aln2.last().copied(),
        ) else {
            return None;
        };
        let start1 = start1 != b'x';
        let start2 = start2 != b'x';
        let end1 = end1 != b'x';
        let end2 = end2 != b'x';
        if is_reverse {
            if (end1 && end2) || (start1 && start2) {
                return Some(*b"00");
            }
        } else {
            if end1 && start2 {
                return Some(*b"12");
            }
            if end2 && start1 {
                return Some(*b"21");
            }
        }
        None
    }

    /// Find links between haplotypes. Used to resolve tandem duplication ordering.
    #[must_use]
    pub fn get_directed_links<'a>(
        new_reads: &'a BTreeMap<String, Vec<BTreeMap<&str, Vec<VString>>>>,
        raw_read_haps: &ReadFingerprintMap,
        assembled_haps: &'a BTreeMap<VStr<'_>, String>,
        is_reverse: bool,
    ) -> DirectedLinks<'a> {
        let mut nondirected_links = BTreeMap::<(&str, &str), Vec<i32>>::new();
        let mut directed_links = BTreeMap::<(&str, &str), Vec<i32>>::new();
        let mut directed_links_loose = BTreeMap::<&str, Vec<&str>>::new();
        for hap_info in new_reads.values().filter(|x| x.len() >= 2) {
            let mut links_found = BTreeSet::<(&str, &str)>::new();
            // All-pairs
            let nseg = hap_info.len();
            for i in 0..nseg {
                for j in (i + 1)..nseg {
                    let left_seg = &hap_info[i];
                    let right_seg = &hap_info[j];
                    let Some(info1) = left_seg.iter().next() else {
                        continue;
                    };
                    let hap1_len_is_1 = info1.1.len() == 1;
                    let Some(aln1) = raw_read_haps
                        .get(&ReadAlignmentId::from_name(*info1.0))
                        .map(|x| x.vstr())
                    else {
                        continue;
                    };
                    log::trace!(
                        "[get_directed_links] Segment 1 candidate: info={:?}, alignment={:?}",
                        info1,
                        aln1,
                    );

                    let Some(info2) = right_seg.iter().next() else {
                        continue;
                    };
                    let hap2_len_is_1 = info2.1.len() == 1;
                    let Some(aln2) = raw_read_haps
                        .get(&ReadAlignmentId::from_name(*info2.0))
                        .map(|x| x.vstr())
                    else {
                        continue;
                    };
                    log::trace!(
                        "[get_directed_links] Segment 2 candidate: info={:?}, alignment={:?}",
                        info2,
                        aln2,
                    );
                    let check_link = Self::check_linking_read(aln1, aln2, is_reverse);
                    if hap1_len_is_1 && hap2_len_is_1 {
                        let hap1 = &info1.1[0];
                        let hap2 = &info2.1[0];
                        if hap1 == hap2 {
                            continue;
                        }
                        let hap1_renamed = &assembled_haps[&hap1.vstr()];
                        let hap2_renamed = &assembled_haps[&hap2.vstr()];
                        let link_to_add = (&hap1_renamed[..], &hap2_renamed[..]);
                        if !links_found.contains(&link_to_add) {
                            log::debug!(
                                "[get_directed_links] Adding non-directed link from read {:?}: {:?}",
                                info1.0,
                                link_to_add
                            );
                            links_found.insert(link_to_add);
                            nondirected_links.entry(link_to_add).or_default().push(1);
                        }
                        if let Some(check_link) = check_link {
                            match check_link {
                            [b'1', b'2'] | [b'0', b'0'] => {
                                directed_links
                                    .entry((hap1_renamed, hap2_renamed))
                                    .or_default()
                                    .push(1);
                                directed_links_loose.entry(hap1_renamed).or_default()
                                    .push(hap2_renamed);
                            }
                            [b'2', b'1'] /*| [b'0', b'0'] */ => { // 00 case in paraphase is unreachable because of prior case.
                                directed_links
                                    .entry((hap2_renamed, hap1_renamed))
                                    .or_default()
                                    .push(1);
                                directed_links_loose.entry(hap2_renamed).or_default()
                                    .push(hap1_renamed);
                            }
                            _ => {}
                            }
                        }
                    } else if check_link.is_some() && (hap1_len_is_1 || hap2_len_is_1) {
                        // all these steps are no-ops if check link is None
                        for (hap1, hap2) in
                            iproduct!(info1.1.iter(), info2.1.iter()).filter(|(x, y)| x != y)
                        {
                            let hap1_renamed = &assembled_haps[&hap1.vstr()];
                            let hap2_renamed = &assembled_haps[&hap2.vstr()];
                            match check_link {
                            Some([b'1', b'2'] | [b'0', b'0']) => {
                                directed_links_loose
                                    .entry(hap1_renamed)
                                    .or_default()
                                    .push(hap2_renamed);
                            }
                            Some([b'2', b'1']) /*| Some([b'0', b'0']) */ => { // 00 case in paraphase is unreachable because of prior case.
                                directed_links_loose
                                    .entry(hap2_renamed)
                                    .or_default()
                                    .push(hap1_renamed);
                            }
                            _ => {}
                        }
                        }
                    }
                }
            }
        }

        DirectedLinks(nondirected_links, directed_links, directed_links_loose)
    }

    /// Generate alleles from linked haplotypes.
    fn get_alleles_from_links<'a>(
        links: &'a [(&'a str, Vec<&'a str>)],
        assembled_haps: &BTreeMap<VStr<'_>, String>,
    ) -> Vec<Vec<&'a str>> {
        let mut alleles: Vec<Vec<&str>> = Vec::new();
        if !links.is_empty() {
            let mut first_allele = vec![links[0].0];
            first_allele.extend_from_slice(&links[0].1[..]);
            for (hap1, hap2) in links
                .iter()
                .flat_map(|(hap1, hap2s)| hap2s.iter().map(move |hap2| (hap1, hap2)))
            {
                let hap1_contained = alleles
                    .iter()
                    .filter(|x| x.contains(hap1))
                    .collect::<Vec<_>>();
                let hap2_contained = alleles
                    .iter()
                    .filter(|x| x.contains(hap2))
                    .collect::<Vec<_>>();
                match (hap1_contained.is_empty(), hap2_contained.is_empty()) {
                    (true, true) => {
                        alleles.push(vec![hap1, hap2]);
                    }
                    (true, false) => {
                        for allele in alleles.iter_mut().filter(|x| x.contains(hap2)) {
                            if !allele.contains(hap1) {
                                allele.push(hap1);
                            }
                            if !allele.contains(hap2) {
                                allele.push(hap2);
                            }
                        }
                    }

                    _ => {
                        // (false, true) + (false, false)
                        for allele in alleles.iter_mut().filter(|x| x.contains(hap1)) {
                            if !allele.contains(hap2) {
                                allele.push(hap2);
                            }
                            if !allele.contains(hap1) {
                                allele.push(hap1);
                            }
                        }
                    }
                }
            }
        }
        // merge alleles
        log::debug!(
            "[get_alleles_from_links] Allele groups before merge pass: {:?}",
            alleles
        );
        loop {
            let mut to_merge = Vec::new();
            for hap in assembled_haps.values() {
                let mut hap_found_in_alleles = Vec::new();
                for allele in &alleles {
                    hap_found_in_alleles.push(allele.contains(&&hap[..]));
                }
                if hap_found_in_alleles
                    .iter()
                    .filter(|a| **a)
                    .collect::<Vec<_>>()
                    .len()
                    > 1
                {
                    to_merge.push(hap);
                    break;
                }
            }
            log::debug!(
                "[get_alleles_from_links] Haplotypes that trigger merging this pass: {:?}",
                to_merge
            );
            if to_merge.is_empty() {
                break;
            }
            let mut new_alleles = Vec::with_capacity(alleles.len());
            let hap = &to_merge[0][..];
            let mut merged = Vec::new();
            for allele_set in &alleles {
                if !allele_set.contains(&hap) {
                    new_alleles.push(allele_set.clone());
                } else {
                    for a in allele_set {
                        merged.push(*a);
                    }
                }
            }
            merged.sort_unstable();
            merged.dedup();
            log::debug!(
                "[get_alleles_from_links] Merged allele group for current pass: {:?}",
                merged
            );
            new_alleles.push(merged);
            alleles = new_alleles;
        }

        alleles
    }

    /// Identify cases where all haplotypes are phased onto one allele
    /// Exclude sex chromosomes in males
    pub fn all_haps_phased_onto_one_allele(
        &self,
        alleles: &Vec<Vec<String>>,
        assembled_haps: &BTreeMap<VStr<'_>, String>,
    ) -> Result<bool, DError> {
        let nchr = self.chr().ok_or_else(|| {
            crate::phaser::Exception::new(format!(
                "Failed to determine chromosome while evaluating allele phasing for gene '{}'",
                self.gene_name()
            ))
        })?;
        let sample_sex = self.settings.sample_sex;
        if alleles.len() == 1
            && !(sample_sex == Sex::Male && (nchr.contains("X") || nchr.contains("Y")))
        {
            let mut first_allele = alleles[0].clone();
            first_allele.sort();
            let mut all_haps = assembled_haps.values().cloned().collect::<Vec<_>>();
            all_haps.sort();
            let matching = first_allele
                .iter()
                .zip(&all_haps)
                .filter(|&(a, b)| a == b)
                .count();
            if first_allele.len() == all_haps.len() && matching == first_allele.len() {
                return Ok(true);
            }
        }
        Ok(false)
    }

    /// Filter potentially inconsistent allele groupings by heuristic rules.
    ///
    /// Handles over-clustered alleles, weakly supported links, and explicit
    /// haplotype exclusions from read-level evidence.
    fn filter_alleles(
        &self,
        raw_alleles: &Vec<Vec<String>>,
        assembled_haps: &BTreeMap<VStr<'_>, String>,
        haps_to_exclude: Vec<VString>,
    ) -> Vec<Vec<String>> {
        let mut alleles = raw_alleles.clone();
        // 1. more than two alleles
        let allele_len = raw_alleles.len();
        if allele_len > 2 {
            return vec![];
        }
        // # 2. two alleles sharing haplotypes
        if allele_len == 2 {
            let first_allele = &raw_alleles[0];
            let second_allele = &raw_alleles[1];
            let hap_overlap = first_allele
                .iter()
                .filter(|x| second_allele.contains(*x))
                .count();
            if hap_overlap > 0 {
                return vec![];
            }
        }
        // 3. not all haplotypes are included in alleles
        let mut haps_in_alleles = Vec::new();
        for allele in raw_alleles {
            for hap in allele {
                haps_in_alleles.push(hap);
            }
        }
        let haps_in_alleles_len = haps_in_alleles.len();
        let mut haps_considered = Vec::new();
        for (hap_seq, hap_name) in assembled_haps {
            let hap_seq_vstring: VString = (*hap_seq).into();
            if !haps_to_exclude.contains(&hap_seq_vstring) {
                haps_considered.push(hap_name.clone());
            }
        }
        let haps_considered_len = haps_considered.len();
        if allele_len == 2 && haps_in_alleles_len < haps_considered_len {
            return vec![];
        }
        if allele_len == 1 && haps_in_alleles_len < haps_considered_len {
            if haps_in_alleles_len == haps_considered_len - 1
                && (haps_in_alleles_len == 2 || haps_in_alleles_len == 3)
            {
                // 1/2 or 1/3 are okay, add the single haplotype as the second allele
                let last_hap = haps_considered
                    .iter()
                    .filter(|x| !haps_in_alleles.contains(x))
                    .collect::<Vec<_>>()[0]
                    .to_string();
                log::debug!(
                    "[filter_alleles] Promoting singleton haplotype {:?} into a second allele",
                    last_hap.clone()
                );
                alleles.push(vec![last_hap]);
            } else {
                return vec![];
            }
        }
        alleles
    }

    /// Returns alleles + read links.
    /// Phase assembled haplotypes into allele groups from read-link evidence.
    ///
    /// Returns both strict and raw allele groupings plus link diagnostics.
    pub fn phase_alleles(
        &self,
        result: &mut PhasedResult,
        assembled_haps: &BTreeMap<VStr<'_>, String>,
        haps_to_exclude: Option<Vec<VString>>,
    ) -> AlleleResult {
        let min_read = 2i32;
        let is_reverse = self.is_reverse();
        let haps_to_exclude = haps_to_exclude.unwrap_or_default();
        let mut new_reads = BTreeMap::<String, Vec<BTreeMap<&str, Vec<VString>>>>::new();
        // unique
        for (hap, reads) in &result.uniquely_supporting_reads {
            if !haps_to_exclude.contains(hap) {
                for read in reads {
                    let short_name = read.split_terminator("_sup").next().unwrap_or(read);
                    new_reads.entry(short_name.into()).or_default().push(
                        [(&read[..], &[hap])]
                            .into_iter()
                            .map(|(name, haplist)| {
                                (name, haplist.iter().copied().cloned().collect::<Vec<_>>())
                            })
                            .collect::<BTreeMap<_, _>>(),
                    );
                }
            }
        }
        for (read_name, haps) in &result.nonuniquely_supporting_reads {
            let short_name = read_name
                .split_terminator("_sup")
                .next()
                .unwrap_or(read_name);
            let mut sub = BTreeMap::new();
            let new_supported_haps = haps
                .iter()
                .filter(|x| !haps_to_exclude.contains(*x))
                .cloned()
                .collect::<Vec<_>>();
            sub.insert(&read_name[..], new_supported_haps);
            new_reads.entry(short_name.into()).or_default().push(sub);
        }

        let links = Self::get_directed_links(
            &new_reads,
            &result.raw_read_haps,
            assembled_haps,
            is_reverse,
        );
        log::debug!(
            "[phase_alleles] Directed-link summary before normalization: {:?}",
            links
        );
        let DirectedLinks(nondir, dir_links, dir_loose) = links;
        log::debug!(
            "[phase_alleles] Non-directed links that passed initial extraction: {:?}",
            nondir
        );

        let mut directed_links = BTreeMap::new();
        for ((a1, a2), b) in dir_links {
            directed_links.insert((a1.to_string(), a2.to_string()), b.len());
        }
        let mut directed_links_loose = BTreeMap::new();
        for (a, b) in dir_loose {
            let counter = b
                .iter()
                .map(|x| x.to_string())
                .collect::<counter::Counter<String, i32>>()
                .into_iter()
                .map(|(k, v)| (k.to_string(), v))
                .collect::<BTreeMap<String, i32>>();
            directed_links_loose.insert(a.to_string(), counter);
        }

        let mut haplotype_links = BTreeMap::<&str, HashSet<&str>>::new();
        for (hap1, hap2) in nondir.iter().filter_map(|x| {
            if x.1.len() as i32 >= min_read {
                Some(x.0)
            } else {
                None
            }
        }) {
            haplotype_links.entry(hap1).or_default().insert(hap2);
            haplotype_links.entry(hap2).or_default().insert(hap1);
        }

        let haplotype_links = haplotype_links
            .into_iter()
            .map(|(x, y)| (x, y.into_iter().collect::<Vec<_>>()))
            .sorted_by_key(|x| std::cmp::Reverse(x.1.len()))
            .collect::<Vec<(&str, Vec<&str>)>>();
        let raw_alleles = Self::get_alleles_from_links(&haplotype_links, assembled_haps)
            .into_iter()
            .map(|v| {
                v.into_iter()
                    .map(std::borrow::ToOwned::to_owned)
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        let haplotype_links = haplotype_links
            .into_iter()
            .map(|(k, v)| {
                (
                    k.to_owned(),
                    v.into_iter()
                        .map(std::borrow::ToOwned::to_owned)
                        .collect::<Vec<_>>(),
                )
            })
            .collect::<BTreeMap<_, _>>();

        log::debug!(
            "[phase_alleles] Raw allele groups before final filtering: {raw_alleles:?}. Haplotype read-link graph: {haplotype_links:?}"
        );
        let alleles = self.filter_alleles(&raw_alleles, assembled_haps, haps_to_exclude);
        AlleleResult {
            alleles,
            raw_alleles,
            haplotype_links,
            directed_links,
            directed_links_loose,
        }
    }
} // impl Phaser

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::json::ReadAlignmentId;

    #[test]
    fn check_linking_read_matches_python_cases() {
        assert_eq!(
            phaser::Phaser::check_linking_read(VStr::from("xxx12"), VStr::from("11xxxx"), false),
            Some(*b"12")
        );
        assert_eq!(
            phaser::Phaser::check_linking_read(VStr::from("xxx12"), VStr::from("11xxxx"), true),
            None
        );

        assert_eq!(
            phaser::Phaser::check_linking_read(VStr::from("11xxxx"), VStr::from("xxx12"), false),
            Some(*b"21")
        );
        assert_eq!(
            phaser::Phaser::check_linking_read(VStr::from("11xxxx"), VStr::from("xxx12"), true),
            None
        );

        assert_eq!(
            phaser::Phaser::check_linking_read(VStr::from("xxx12"), VStr::from("xxxxxx11"), false),
            None
        );
        assert_eq!(
            phaser::Phaser::check_linking_read(VStr::from("xxx12"), VStr::from("xxxxxx11"), true),
            Some(*b"00")
        );

        assert_eq!(
            phaser::Phaser::check_linking_read(VStr::from("xxx1x"), VStr::from("1xxxxx11"), false),
            None
        );
        assert_eq!(
            phaser::Phaser::check_linking_read(VStr::from("xxx1x"), VStr::from("1xxxxx11"), true),
            None
        );
    }

    #[test]
    fn get_alleles_from_links_matches_python_cases() {
        let assembled_haps = BTreeMap::from([
            (VStr::from("1"), String::from("1")),
            (VStr::from("2"), String::from("2")),
            (VStr::from("3"), String::from("3")),
            (VStr::from("4"), String::from("4")),
            (VStr::from("5"), String::from("5")),
            (VStr::from("6"), String::from("6")),
            (VStr::from("7"), String::from("7")),
        ]);

        let links = vec![("1", vec!["2"]), ("2", vec!["1", "3"]), ("3", vec!["2"])];
        let alleles = phaser::Phaser::get_alleles_from_links(&links, &assembled_haps);
        assert_eq!(alleles, vec![vec!["1", "2", "3"]]);

        let links = vec![
            ("1", vec!["2"]),
            ("2", vec!["1", "3"]),
            ("4", vec!["5"]),
            ("5", vec!["4", "3"]),
            ("3", vec!["5"]),
            ("6", vec!["7"]),
            ("7", vec!["6"]),
        ];
        let alleles = phaser::Phaser::get_alleles_from_links(&links, &assembled_haps);
        assert_eq!(alleles.len(), 2);
        assert!(alleles.iter().any(|a| {
            let mut x = a.clone();
            x.sort_unstable();
            x == vec!["1", "2", "3", "4", "5"]
        }));
    }

    #[test]
    fn get_directed_links_matches_python_case() {
        let raw_read_haps = ReadFingerprintMap::from([
            (ReadAlignmentId::from_name("r1"), VString::from("xx1")),
            (ReadAlignmentId::from_name("r1_sup"), VString::from("2xx")),
            (ReadAlignmentId::from_name("r2"), VString::from("xx1")),
            (ReadAlignmentId::from_name("r2_sup"), VString::from("2xx")),
            (ReadAlignmentId::from_name("r3"), VString::from("x1x")),
            (ReadAlignmentId::from_name("r3_sup"), VString::from("2xx")),
        ]);
        let new_reads = BTreeMap::from([
            (
                String::from("r1"),
                vec![
                    BTreeMap::from([("r1", vec![VString::from("111")])]),
                    BTreeMap::from([("r1_sup", vec![VString::from("211")])]),
                ],
            ),
            (
                String::from("r2"),
                vec![
                    BTreeMap::from([("r2", vec![VString::from("111")])]),
                    BTreeMap::from([("r2_sup", vec![VString::from("211")])]),
                ],
            ),
            (
                String::from("r3"),
                vec![
                    BTreeMap::from([("r3", vec![VString::from("111")])]),
                    BTreeMap::from([("r3_sup", vec![VString::from("211")])]),
                ],
            ),
        ]);
        let assembled_haps = BTreeMap::from([
            (VStr::from("111"), String::from("hap1")),
            (VStr::from("211"), String::from("hap2")),
        ]);

        let DirectedLinks(nondirected_links, directed_links, _loose) =
            phaser::Phaser::get_directed_links(&new_reads, &raw_read_haps, &assembled_haps, false);

        assert_eq!(
            directed_links
                .iter()
                .map(|((a, b), v)| (format!("{a}-{b}"), v.clone()))
                .collect::<BTreeMap<_, _>>(),
            BTreeMap::from([(String::from("hap1-hap2"), vec![1, 1])])
        );
        assert_eq!(
            nondirected_links
                .iter()
                .map(|((a, b), v)| (format!("{a}-{b}"), v.clone()))
                .collect::<BTreeMap<_, _>>(),
            BTreeMap::from([(String::from("hap1-hap2"), vec![1, 1, 1])])
        );
    }
}
