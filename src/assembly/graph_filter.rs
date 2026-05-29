use crate::assembly::assembly_result::AssembledPaths;
use crate::assembly::graph_stats::compute_median;
use crate::assembly::variant_graph::{Graph, ReadHapPair};
use crate::io::json::ReadFingerprintMap;
use crate::toolkit::hapcmp::HapCompare;
use crate::toolkit::util::DError;
use vstr::{VStr, VString};

use itertools::Itertools;

use std::collections::{BTreeMap, BTreeSet};

fn trimmed_x_len(seq: &[u8]) -> usize {
    let start = seq.iter().position(|&x| x != b'x');
    let end = seq.iter().rposition(|&x| x != b'x');
    match (start, end) {
        (Some(s), Some(e)) if e >= s => e - s + 1,
        _ => 0,
    }
}

impl Graph {
    /// Choose the haplotypes with the best support.
    #[must_use]
    pub fn select_main_haps(
        &self,
        final_haps: &AssembledPaths,
        mut candidates: Vec<AssembledPaths>,
    ) -> AssembledPaths {
        log::debug!(
            "Selecting best assembled paths from final_haps={final_haps:?} using pivot candidates={candidates:?}; pivot_index={}",
            self.settings.pivot_index
        );

        let ret = if self.settings.pivot_index < 0 {
            assert!(!candidates.is_empty(), "No extended candidates");
            std::mem::take(&mut candidates[0])
        } else {
            let best_scores = candidates
                .into_iter()
                .map(|cand| {
                    log::trace!(
                        "Pivot-nonmissing haplotypes for candidate {cand:?}: {:?}. pivot_index={}",
                        cand.clone()
                            .into_inner()
                            .into_iter()
                            .filter(|hap| hap[self.settings.pivot_index as usize] != b'x')
                            .collect::<Vec<_>>(),
                        self.settings.pivot_index
                    );
                    let pivot_ones = cand
                        .into_inner()
                        .into_iter()
                        .filter(|hap| hap[self.settings.pivot_index as usize] != b'x')
                        .collect::<Vec<_>>();
                    let nhap = pivot_ones.len();
                    let total_length = pivot_ones
                        .iter()
                        .map(|seq| trimmed_x_len(seq))
                        .sum::<usize>();
                    let num_matching_fhaps = final_haps
                        .iter()
                        .filter(|fhap| {
                            pivot_ones
                                .iter()
                                .filter_map(|x| {
                                    HapCompare::from_haps(fhap.vstr(), x.vstr())
                                        .ok()
                                        .map(|cmp| cmp.mismatches)
                                })
                                .min()
                                .unwrap_or(1)
                                == 0
                        })
                        .count();
                    ((nhap, num_matching_fhaps, total_length), pivot_ones)
                })
                .sorted_by_key(|x| std::cmp::Reverse(x.0))
                .collect::<Vec<_>>();
            log::trace!("Sorted candidate scores for pivot selection: {best_scores:?}");
            if let Some(best) = best_scores.first().map(|x| x.0) {
                AssembledPaths::from_seqs(
                    best_scores
                        .into_iter()
                        .filter_map(|x| if x.0 == best { Some(x.1) } else { None })
                        .min_by_key(std::vec::Vec::len)
                        .unwrap_or_default(),
                )
            } else {
                final_haps.clone()
            }
        };
        if ret.is_empty() {
            final_haps.clone()
        } else {
            ret
        }
    }

    /// Return the highest possible copy number given the assembled haplotypes
    /// Returns (highest count, `candidate_hap_sets`)
    pub(crate) fn get_highest_cn<T>(
        &self,
        paths: impl IntoIterator<Item = T>,
    ) -> Result<(usize, Vec<Vec<VString>>), DError>
    where
        T: std::convert::Into<VString>,
    {
        let seqs = paths
            .into_iter()
            .map(std::convert::Into::into)
            .collect::<Vec<_>>();
        let mut cand = Vec::with_capacity(self.nvar / 2);
        for set in (0..self.nvar)
            .rev()
            .map(|idx| seqs.iter().filter(|x| x[idx] != b'x').collect::<Vec<_>>())
        {
            if !cand.contains(&set) {
                cand.push(set);
            }
        }
        let cand = cand
            .into_iter()
            .sorted_by_key(|x| {
                (
                    std::cmp::Reverse(x.len()),
                    x.iter()
                        .map(|hap| hap.iter().filter(|&n| *n == b'x').collect::<Vec<_>>().len())
                        .sum::<usize>(),
                )
            })
            .map(|x| x.into_iter().map(VString::from).collect::<Vec<_>>())
            .collect::<Vec<_>>();
        log::trace!("[get_highest_cn] candidate haplotype sets: {cand:?}");
        Ok(cand
            .first()
            .map(std::vec::Vec::len)
            .map(|x| (x, cand))
            .ok_or_else(|| {
                anyhow::anyhow!("get_highest_cn failed because there were no passing haps")
            })?)
    }

    /// Choose the best candidates from an input set.
    #[must_use]
    pub(crate) fn filter_candidates(
        &self,
        hap: VStr<'_>,
        candidates: &BTreeSet<VString>,
        len_x: usize,
        candidate_first: bool,
    ) -> Vec<VString> {
        let mut read_support = BTreeMap::<VStr<'_>, Vec<i32>>::new();
        let mut ret = vec![];
        for c in candidates {
            let len = c.len();
            let mut hap2 = VString::from(vec![
                b'x';
                if candidate_first {
                    len_x - len
                } else {
                    self.nvar - len_x
                }
            ]);
            hap2.extend_from_slice(&c[..]);
            hap2.resize(self.nvar, b'x');
            let hap2 = hap2;
            for read_seq in self.reads.values() {
                let Ok(cmp1) = HapCompare::from_haps(hap, read_seq) else {
                    continue;
                };
                let Ok(cmp2) = HapCompare::from_haps(&hap2, read_seq) else {
                    continue;
                };
                if cmp1.mismatches == 0
                    && cmp2.mismatches == 0
                    && cmp1.matches > 0
                    && cmp2.matches > 0
                {
                    read_support
                        .entry(c.vstr())
                        .or_default()
                        .push(std::cmp::min(cmp1.matches, cmp2.matches));
                }
            }
        }
        log::trace!("read_support: {read_support:?}");
        match read_support.len() {
            0 => {
                return ret;
            }
            1 => {
                if let Some((first_hap, _)) = read_support.iter().next() {
                    ret.push((*first_hap).into());
                }
                return ret;
            }
            _ => {}
        }
        let medians_and_maxes = read_support
            .values()
            .map(|x| {
                (
                    compute_median(x),
                    x.iter().copied().max().map_or(f32::NAN, |x| x as f32),
                )
            })
            .collect::<Vec<(f32, f32)>>();
        for ((max_other, median, id), (hap, read_ids)) in medians_and_maxes
            .iter()
            .map(|x| x.0)
            .enumerate()
            .map(|(id, median)| {
                let max_other = medians_and_maxes
                    .iter()
                    .map(|(_median, max)| max)
                    .copied()
                    .enumerate()
                    .filter_map(
                        |(other_id, max)| {
                            if other_id == id {
                                None
                            } else {
                                Some(max)
                            }
                        },
                    )
                    .max_by(|x, y| x.partial_cmp(y).unwrap_or(std::cmp::Ordering::Equal))
                    .unwrap_or(f32::NAN);
                (max_other, median, id)
            })
            .zip(read_support.iter())
        {
            log::trace!("median: {median}. max other: {max_other} for id {id}");
            if median > max_other && median >= 4. && read_ids.len() >= 4 {
                ret.push(hap.into());
                break;
            }
        }
        if ret.is_empty() {
            ret = candidates.iter().cloned().collect::<Vec<_>>();
        }
        ret
    }

    /// Filter low-support haplotypes.
    /// Takes an `AssembledPaths` object and returns a `Result<AssembledPaths, DError>`
    /// where `DError` is an alias for `Box<dyn std::error::Error>`.
    pub fn filter_low_support_haps(
        &self,
        init_haps: &AssembledPaths,
    ) -> Result<AssembledPaths, DError> {
        Self::filter_low_support_haps_detail(
            init_haps,
            self.settings.min_hap_support as usize,
            &self.reads_original,
        )
    }

    /// Iteratively remove haplotypes without sufficient unique read support.
    ///
    /// Recomputes support after each filtering round until the surviving set
    /// stabilizes.
    pub fn filter_low_support_haps_detail(
        init_haps: &AssembledPaths,
        min_count: usize,
        reads: &ReadFingerprintMap,
    ) -> Result<AssembledPaths, DError> {
        log::trace!(
            "Starting low-support haplotype filtering. min_count={min_count}, haps={init_haps:?}, reads={reads:?}"
        );
        let flat_init_haps = init_haps
            .iter()
            .map(|x| x.vstr())
            .collect::<Vec<VStr<'_>>>();
        let read_hap_support =
            Self::match_reads_and_haps(reads, &flat_init_haps[..], /* min_match= */ None)?;
        let good_reads = read_hap_support.support_reads;
        log::trace!("Initial read-to-haplotype support map: {good_reads:?}");
        let mut filtered_ass_haps = good_reads
            .into_iter()
            .filter(|(_id, support)| !support.is_empty())
            .collect::<Vec<_>>();
        log::trace!(
            "{} haplotypes passed initial support screening out of {} candidates (min_count={min_count})",
            filtered_ass_haps.len(),
            flat_init_haps.len()
        );
        let mut iternum = 0;
        loop {
            let haps_to_assess = filtered_ass_haps
                .iter()
                .map(|x| flat_init_haps[x.0 as usize])
                .collect::<Vec<_>>();
            log::trace!("Iteration {iternum}: assessing haplotypes {haps_to_assess:?}");
            let read_hap_support =
                Self::match_reads_and_haps(reads, &haps_to_assess[..], /* min_match= */ None)?;
            log::trace!(
                "Number of uniquely-matching reads: {}. Support: {:?}",
                read_hap_support.support_reads.len(),
                read_hap_support.support_reads
            );
            let filtered_indices = filtered_ass_haps
                .iter()
                .map(|x| x.0 as usize)
                .collect::<Vec<_>>();
            let human_readable_support = read_hap_support
                .support_reads
                .iter()
                .map(|(k, v)| (flat_init_haps[filtered_indices[*k as usize]], &v[..]))
                .collect::<BTreeMap<VStr<'_>, &[ReadHapPair<'_>]>>();
            log::trace!(
                "Iteration {iternum}: support by haplotype sequence: {human_readable_support:?}"
            );
            let local_passing_haps = read_hap_support.support_reads.iter().filter_map(|(k, v)| {
                if v.len() >= min_count {
                    Some(k)
                } else {
                    log::trace!(
                        "Dropping hap index {k} ({}) due to insufficient support: count={}, supports={v:?}",
                        v.len(),
                        flat_init_haps[filtered_indices[*k as usize]]
                    );
                    None
                }
            });
            log::trace!(
                "Iteration {iternum}: passing hap indices after support cutoff: {local_passing_haps:?}"
            );
            let offset_passing_haps = local_passing_haps
                .map(|x| std::mem::take(&mut filtered_ass_haps[*x as usize]))
                .collect::<Vec<_>>();
            filtered_ass_haps = offset_passing_haps;
            log::trace!(
                "Iteration {iternum}: filtered hap entries after applying support cutoff: {filtered_ass_haps:?}"
            );
            iternum += 1;
            log::trace!(
                "{} haplotypes remain after {iternum} iterations",
                filtered_ass_haps.len()
            );
            if filtered_ass_haps.len() == haps_to_assess.len() {
                log::trace!(
                    "Filtered set stabilized at iteration {iternum}. Filtered={filtered_ass_haps:?}, assessed={haps_to_assess:?}. Stopping."
                );
                break;
            }
        }
        let filtered_ass_haps = AssembledPaths::from_seqs(
            filtered_ass_haps
                .into_iter()
                .map(|x| flat_init_haps[x.0 as usize]),
        );
        log::debug!(
            "{} haplotypes remain after filtering. Output={filtered_ass_haps:?}, input_count={}, input={init_haps:?}",
            filtered_ass_haps.len(),
            init_haps.len(),
        );
        Ok(filtered_ass_haps)
    }
}

#[cfg(test)]
mod tests {
    use super::trimmed_x_len;

    #[test]
    fn trimmed_x_len_handles_edge_and_mixed_cases() {
        assert_eq!(trimmed_x_len(b""), 0);
        assert_eq!(trimmed_x_len(b"x"), 0);
        assert_eq!(trimmed_x_len(b"xxxx"), 0);
        assert_eq!(trimmed_x_len(b"abc"), 3);
        assert_eq!(trimmed_x_len(b"xabc"), 3);
        assert_eq!(trimmed_x_len(b"abcx"), 3);
        assert_eq!(trimmed_x_len(b"xxabcxx"), 3);
        assert_eq!(trimmed_x_len(b"axxxb"), 5);
    }
}
