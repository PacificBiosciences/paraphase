//! Deletion discovery and labeling helpers for [`Phaser`].

use crate::phaser::{check_del, fiveprime_clip_length, threeprime_clip_length};
use crate::phaser::{Exception, Phaser};
use crate::toolkit::deletion::Datum as DeletionDatum;
use crate::toolkit::range;
use crate::toolkit::util::{DError, DResult};

use rust_htslib::{bam, bam::ext::BamRecordExtensions, bam::Read};

use std::collections::BTreeSet;

impl Phaser {
    /// Parse predefined deletions from region config.
    pub fn parse_deletions_from_config(&mut self) -> DResult {
        let deletion1_name = self.locus_config().get("deletion1_name");
        if deletion1_name.is_some() {
            let deletion1_name = deletion1_name
                .and_then(|x| x.as_str())
                .ok_or_else(|| {
                    Exception::new(format!(
                        "Missing deletion1_name in region config for gene '{}'",
                        self.gene_name()
                    ))
                })?
                .to_string();
            let deletion1_start = deletion1_name
                .split("_")
                .map(std::borrow::ToOwned::to_owned)
                .collect::<Vec<_>>()
                .first()
                .ok_or_else(|| {
                    Exception::new(format!(
                        "Malformed deletion1_name '{}'; expected 'start_*'",
                        deletion1_name
                    ))
                })?
                .parse::<i64>()?
                - 1;
            let deletion1_size = self
                .locus_config()
                .get("deletion1_size")
                .and_then(|x| {
                    x.as_str()
                        .and_then(|s| s.parse::<i64>().ok())
                        .or(x.as_i64())
                })
                .ok_or_else(|| {
                    Exception::new(format!(
                        "Missing deletion1_size in region config for gene '{}'",
                        self.gene_name()
                    ))
                })?;
            let del1_3p_pos1 = self
                .locus_config()
                .get("del1_3p_pos1")
                .and_then(|x| {
                    x.as_str()
                        .and_then(|s| s.parse::<i64>().ok())
                        .or(x.as_i64())
                })
                .ok_or_else(|| {
                    Exception::new(format!(
                        "Missing del1_3p_pos1 in region config for gene '{}'",
                        self.gene_name()
                    ))
                })?;
            let del1_3p_pos2 = self
                .locus_config()
                .get("del1_3p_pos2")
                .and_then(|x| {
                    x.as_str()
                        .and_then(|s| s.parse::<i64>().ok())
                        .or(x.as_i64())
                })
                .ok_or_else(|| {
                    Exception::new(format!(
                        "Missing del1_3p_pos2 in region config for gene '{}'",
                        self.gene_name()
                    ))
                })?;
            let del1_5p_pos1 = self
                .locus_config()
                .get("del1_5p_pos1")
                .and_then(|x| {
                    x.as_str()
                        .and_then(|s| s.parse::<i64>().ok())
                        .or(x.as_i64())
                })
                .ok_or_else(|| {
                    Exception::new(format!(
                        "Missing del1_5p_pos1 in region config for gene '{}'",
                        self.gene_name()
                    ))
                })?;
            let del1_5p_pos2 = self
                .locus_config()
                .get("del1_5p_pos2")
                .and_then(|x| {
                    x.as_str()
                        .and_then(|s| s.parse::<i64>().ok())
                        .or(x.as_i64())
                })
                .ok_or_else(|| {
                    Exception::new(format!(
                        "Missing del1_5p_pos2 in region config for gene '{}'",
                        self.gene_name()
                    ))
                })?;

            let deletion_end = deletion1_start + deletion1_size;
            let new_del_data = DeletionDatum::new(
                range::I64::new(deletion1_start, deletion_end),
                None,
                Some(del1_3p_pos1),
                Some(del1_3p_pos2),
                Some(del1_5p_pos1),
                Some(del1_5p_pos2),
            );
            self.del_data.push(new_del_data);
        }
        let deletion2_name = self.locus_config().get("deletion2_name");
        if deletion2_name.is_some() {
            let deletion2_name = deletion2_name
                .and_then(|x| x.as_str())
                .ok_or_else(|| {
                    Exception::new(format!(
                        "Missing deletion2_name in region config for gene '{}'",
                        self.gene_name()
                    ))
                })?
                .to_string();
            let deletion2_start = deletion2_name
                .split("_")
                .map(std::borrow::ToOwned::to_owned)
                .collect::<Vec<_>>()
                .first()
                .ok_or_else(|| {
                    Exception::new(format!(
                        "Malformed deletion2_name '{}'; expected 'start_*'",
                        deletion2_name
                    ))
                })?
                .parse::<i64>()?
                - 1;
            let deletion2_size = self
                .locus_config()
                .get("deletion2_size")
                .and_then(|x| {
                    x.as_str()
                        .and_then(|s| s.parse::<i64>().ok())
                        .or(x.as_i64())
                })
                .ok_or_else(|| {
                    Exception::new(format!(
                        "Missing deletion2_size in region config for gene '{}'",
                        self.gene_name()
                    ))
                })?;
            let del2_3p_pos1 = self
                .locus_config()
                .get("del2_3p_pos1")
                .and_then(|x| {
                    x.as_str()
                        .and_then(|s| s.parse::<i64>().ok())
                        .or(x.as_i64())
                })
                .ok_or_else(|| {
                    Exception::new(format!(
                        "Missing del2_3p_pos1 in region config for gene '{}'",
                        self.gene_name()
                    ))
                })?;
            let del2_3p_pos2 = self
                .locus_config()
                .get("del2_3p_pos2")
                .and_then(|x| {
                    x.as_str()
                        .and_then(|s| s.parse::<i64>().ok())
                        .or(x.as_i64())
                })
                .ok_or_else(|| {
                    Exception::new(format!(
                        "Missing del2_3p_pos2 in region config for gene '{}'",
                        self.gene_name()
                    ))
                })?;
            let del2_5p_pos1 = self
                .locus_config()
                .get("del2_5p_pos1")
                .and_then(|x| {
                    x.as_str()
                        .and_then(|s| s.parse::<i64>().ok())
                        .or(x.as_i64())
                })
                .ok_or_else(|| {
                    Exception::new(format!(
                        "Missing del2_5p_pos1 in region config for gene '{}'",
                        self.gene_name()
                    ))
                })?;
            let del2_5p_pos2 = self
                .locus_config()
                .get("del2_5p_pos2")
                .and_then(|x| {
                    x.as_str()
                        .and_then(|s| s.parse::<i64>().ok())
                        .or(x.as_i64())
                })
                .ok_or_else(|| {
                    Exception::new(format!(
                        "Missing del2_5p_pos2 in region config for gene '{}'",
                        self.gene_name()
                    ))
                })?;

            let deletion_end = deletion2_start + deletion2_size;
            let new_del_data = DeletionDatum::new(
                range::I64::new(deletion2_start, deletion_end),
                None,
                Some(del2_3p_pos1),
                Some(del2_3p_pos2),
                Some(del2_5p_pos1),
                Some(del2_5p_pos2),
            );
            self.del_data.push(new_del_data);
        }
        log::debug!(
            "Configured deletion templates parsed from region config: del_data={:?}",
            self.del_data
        );
        Ok(())
    }

    /// Get records overlapping pos1..pos2 and store as a `Vec<bam::Record>`.
    fn get_records_for_range(&mut self, range: &range::I64) -> Result<Vec<bam::Record>, DError> {
        self.get_records(range.start, range.end)
    }

    /// Get records overlapping pos1..pos2 and store as a `Vec<bam::Record>`.
    fn get_records(&mut self, pos1: i64, pos2: i64) -> Result<Vec<bam::Record>, DError> {
        let mut realigned_bam = self.try_realigned_bam()?;
        let tid = self
            .genome_tid()
            .map(|x| x as i32)
            .ok_or_else(|| Exception::new("get_records: missing chr tid in BAM header"))?;
        realigned_bam.fetch((tid, pos1, pos2))?;
        Ok(realigned_bam
            .records()
            .filter_map(Result::ok)
            .collect::<Vec<_>>())
    }

    /// Find big deletions in data and seed deletion candidates.
    pub(crate) fn discover_big_dels(&mut self) -> DResult {
        let mut realigned_bam = self.try_realigned_bam()?;
        let tid = self
            .genome_tid()
            .map(|x| x as i32)
            .ok_or_else(|| Exception::new("discover_big_dels: missing chr tid in BAM header"))?;
        realigned_bam.fetch((
            tid,
            self.left_boundary_0based(),
            self.right_boundary_0based(),
        ))?;
        let del_reads = realigned_bam
            .rc_records()
            .filter_map(Result::ok)
            .map(|x| {
                let max_len = i64::from(
                    x.cigar()
                        .iter()
                        .map(|x| match x {
                            bam::record::Cigar::Del(x) => *x,
                            _ => 0,
                        })
                        .max()
                        .unwrap_or(0),
                );
                log::trace!(
                    "[discover_big_dels] raw deletion candidate read: read_name={}, deletion_length={}",
                    String::from_utf8_lossy(x.qname()),
                    max_len
                );
                (x, max_len)
            })
            .filter(|(_x, max_len)| *max_len >= self.settings.big_deletion_settings.min_size)
            .filter_map(|(x, max_len)| {
                let del_position = x
                    .cigar()
                    .iter()
                    .position(|&x| x == bam::record::Cigar::Del(max_len as u32));
                let Some(del_position) = del_position else {
                    log::warn!(
                        "[discover_big_dels] expected deletion CIGAR op was not found in candidate read: read_name={}",
                        String::from_utf8_lossy(x.qname())
                    );
                    return None;
                };
                let del_pos = x.pos() - 1
                    + x.cigar()
                        .iter()
                        .take(del_position)
                        .map(|x| {
                            if let bam::record::Cigar::Match(_len) = x {
                                i64::from(x.len())
                            } else if let bam::record::Cigar::Del(_len) = x {
                                i64::from(x.len())
                            } else if let bam::record::Cigar::Equal(_len) = x {
                                i64::from(x.len())
                            } else if let bam::record::Cigar::Diff(_len) = x {
                                i64::from(x.len())
                            } else {
                                0
                            }
                        })
                        .sum::<i64>();
                log::trace!(
                    "[discover_big_dels] long deletion coordinates from read: read_name={}, start={}, end={}",
                    String::from_utf8_lossy(x.qname()),
                    del_pos,
                    del_pos + max_len
                );
                Some((del_pos, del_pos + max_len))
            })
            .collect::<Vec<(i64, i64)>>();
        log::debug!(
            "[discover_big_dels] candidate long-deletion intervals observed in reads: del_reads={del_reads:?}"
        );
        let common = del_reads
            .iter()
            .copied()
            .collect::<counter::Counter<(i64, i64), i64>>()
            .k_most_common_ordered(self.settings.max_number_deletions as usize);

        let padding = self.settings.big_deletion_settings.padding;
        self.del_data = common
            .into_iter()
            .filter_map(|x| {
                if x.1 >= self.settings.big_deletion_settings.min_count {
                    Some(x.0)
                } else {
                    None
                }
            })
            .map(|item| {
                let raw = item.into();
                DeletionDatum::new(raw, Some(padding), None, None, None, None)
            })
            .collect::<Vec<_>>();
        log::debug!(
            "[discover_big_dels] retained deletion candidates after frequency/size filters: del_data={:?}",
            self.del_data
        );
        Ok(())
    }

    /// Call recurrent clip sites from data.
    pub fn find_clip_site(
        &mut self,
        min_clip_length: Option<u32>,
        min_count: Option<usize>,
        padding: Option<i64>,
    ) -> DResult {
        let min_clip_length = min_clip_length.unwrap_or(800);
        let min_count = min_count.unwrap_or(6);
        let padding = padding.unwrap_or(1000);
        let mut clip_reads = Vec::new();
        let mut realigned_bam = self.try_realigned_bam()?;
        let tid = self
            .genome_tid()
            .map(|x| x as i32)
            .ok_or_else(|| Exception::new("find_clip_site: missing chr tid in BAM header"))?;
        realigned_bam.fetch((
            tid,
            self.left_boundary_0based(),
            self.right_boundary_0based(),
        ))?;
        let records = realigned_bam
            .records()
            .filter_map(Result::ok)
            .collect::<Vec<_>>();
        for record in records {
            let cigar = record.cigar();
            let clip_len_5p = fiveprime_clip_length(&cigar);
            if clip_len_5p >= min_clip_length {
                let clip_pos = record.reference_start() - 1;
                if clip_pos > self.left_boundary_0based() + padding
                    && clip_pos < self.right_boundary_0based() - padding
                {
                    clip_reads.push((clip_pos, "5p"));
                }
            }
            let clip_len_3p = threeprime_clip_length(&cigar);
            if clip_len_3p >= min_clip_length {
                let clip_pos = record.reference_end() - 1;
                if clip_pos > self.left_boundary_0based() + padding
                    && clip_pos < self.right_boundary_0based() - padding
                {
                    clip_reads.push((clip_pos, "3p"));
                }
            }
        }
        let clip_counter = clip_reads.iter().copied().collect::<counter::Counter<_>>();
        for ((pos, clip_direction), this_count) in clip_counter {
            if this_count >= min_count {
                if clip_direction == "5p" {
                    let too_close = self
                        .clip_5p_positions
                        .iter()
                        .any(|&existing_pos| (pos - 100..pos + 100).contains(&existing_pos));
                    if !too_close {
                        let mut in_deletion = false;
                        for each_known_deletion in &self.del_data {
                            if ranges_overlap(
                                pos - 100,
                                pos + 100,
                                each_known_deletion.fivep().start,
                                each_known_deletion.fivep().end,
                            ) {
                                in_deletion = true;
                            }
                        }
                        if !in_deletion {
                            self.clip_5p_positions.push(pos);
                        }
                    }
                }
                if clip_direction == "3p" {
                    let too_close = self
                        .clip_3p_positions
                        .iter()
                        .any(|&existing_pos| (pos - 100..pos + 100).contains(&existing_pos));
                    if !too_close {
                        let mut in_deletion = false;
                        for each_known_deletion in &self.del_data {
                            if ranges_overlap(
                                pos - 100,
                                pos + 100,
                                each_known_deletion.threep().start,
                                each_known_deletion.threep().end,
                            ) {
                                in_deletion = true;
                            }
                        }
                        if !in_deletion {
                            self.clip_3p_positions.push(pos);
                        }
                    }
                }
            }
        }
        self.clip_5p_positions.sort();
        self.clip_3p_positions.sort();

        Ok(())
    }

    /// Updates `self.del_data` vector with relevant read names per deletion.
    pub(crate) fn label_big_dels(&mut self) -> DResult {
        log::debug!("Labeling reads by deletion-support category (full/partial/negative)");
        let mut fetched_reads = Vec::with_capacity(self.del_data.len());
        for (threep, fivep) in self.del_data.clone().iter().map(|x| {
            let threep = self.get_records_for_range(&x.threep());
            let fivep = self.get_records_for_range(&x.fivep());
            (threep, fivep)
        }) {
            fetched_reads.push((threep?, fivep?));
        }
        let min_clip_len = self.settings.big_deletion_settings.min_clip_len;
        let min_extend = self.settings.big_deletion_settings.min_extend;
        let padding_negative_reads = self.settings.big_deletion_settings.padding_negative_reads;
        let use_supplementary = self.use_supplementary();
        for ((threep, fivep), del_data) in fetched_reads.iter_mut().zip(self.del_data.iter_mut()) {
            let mut p3_reads = BTreeSet::<String>::new();
            let mut p5_reads = BTreeSet::<String>::new();
            let mut del_reads = BTreeSet::<String>::new();
            for read in threep {
                let read_name = Self::get_read_name_free(read, use_supplementary);
                read.cache_cigar();
                if read.reference_start() < del_data.threep().start - padding_negative_reads
                    && read.reference_end() > del_data.threep().end + padding_negative_reads
                {
                    del_data.del_negative_reads.insert(read_name.clone());
                }

                let reference_start_cutoff = del_data.threep().start - min_extend;
                let threep_length = threeprime_clip_length(&read.cigar());
                let end = read.reference_end();

                if (i64::from(threep_length) >= min_clip_len)
                    && (del_data.threep().start < end)
                    && (end < del_data.threep().end)
                    && read.pos() < reference_start_cutoff
                {
                    p3_reads.insert(read_name.clone());
                }
                let has_deletion_in_cigar = check_del(read, del_data.size(), del_data.threep());
                if has_deletion_in_cigar {
                    del_reads.insert(read_name.clone());
                }
            }

            for read in fivep {
                let read_name = Self::get_read_name_free(read, use_supplementary);
                read.cache_cigar();
                if read.reference_start() < del_data.fivep().start - padding_negative_reads
                    && read.reference_end() > del_data.fivep().end + padding_negative_reads
                {
                    del_data.del_negative_reads.insert(read_name.clone());
                }

                let reference_end_cutoff = del_data.fivep().end + min_extend;
                let fivep_length = fiveprime_clip_length(&read.cigar());
                let pos = read.pos();

                if (i64::from(fivep_length) >= min_clip_len)
                    && (del_data.fivep().start < pos)
                    && (pos < del_data.fivep().end)
                    && read.reference_end() > reference_end_cutoff
                {
                    p5_reads.insert(read_name.clone());
                }
            }
            log::debug!("Supporting reads for this deletion: del_reads={del_reads:?}, p3_reads={p3_reads:?}, p5_reads={p5_reads:?}");

            let p3_p5_del_union;
            let full_union;
            if !del_reads.is_empty() || (!p3_reads.is_empty() && !p5_reads.is_empty()) {
                let p3_p5_intersect = p3_reads
                    .intersection(&p5_reads)
                    .cloned()
                    .collect::<BTreeSet<_>>();
                p3_p5_del_union = p3_p5_intersect
                    .union(&del_reads)
                    .cloned()
                    .collect::<BTreeSet<_>>();
                full_union = del_reads
                    .union(&p3_reads)
                    .chain(p5_reads.iter())
                    .cloned()
                    .collect::<BTreeSet<_>>();
            } else {
                p3_p5_del_union = BTreeSet::new();
                full_union = BTreeSet::new();
            }
            del_data.del_reads = p3_p5_del_union;
            del_data.del_reads_partial = full_union;
            log::debug!(
                "Deletion support labels finalized for one candidate: del_data={del_data:?}"
            );
        }
        Ok(())
    }
}

/// Helper function to check if two ranges overlap
fn ranges_overlap(start1: i64, end1: i64, start2: i64, end2: i64) -> bool {
    start1 < end2 && start2 < end1
}

#[cfg(test)]
mod tests {
    use crate::phaser::check_del;
    use crate::toolkit::range::I64 as Range64;
    use rust_htslib::bam::record::{Cigar, CigarString, Record};

    #[test]
    fn test_ranges_overlap() {
        assert!(super::ranges_overlap(0, 100, 50, 150));
        assert!(!super::ranges_overlap(0, 100, 100, 200));
        assert!(super::ranges_overlap(0, 100, 50, 75));
    }

    #[test]
    fn test_check_del() {
        let mut test_record = Record::new();
        test_record.set_pos(200);
        let test_cigar = CigarString(vec![Cigar::Match(151), Cigar::Del(50), Cigar::Match(100)]);
        test_record.set(
            "test_reads".as_bytes(),
            Some(&test_cigar),
            "AAA".as_bytes(),
            "~~~".as_bytes(),
        );
        let test_range = Range64::new(350, 360);
        let deletion_size = 50;
        let deletion_present = check_del(&test_record, deletion_size, test_range);
        assert!(deletion_present);

        // deletion position off
        let mut test_record = Record::new();
        test_record.set_pos(200);
        let test_cigar = CigarString(vec![
            Cigar::Match(200),
            Cigar::SoftClip(50),
            Cigar::Match(100),
        ]);
        test_record.set(
            "test_reads".as_bytes(),
            Some(&test_cigar),
            "AAA".as_bytes(),
            "~~~".as_bytes(),
        );
        let test_range = Range64::new(350, 360);
        let deletion_size = 50;
        let deletion_present = check_del(&test_record, deletion_size, test_range);
        assert!(!deletion_present);

        // deletion size different
        let mut test_record = Record::new();
        test_record.set_pos(200);
        let test_cigar = CigarString(vec![Cigar::Match(151), Cigar::Del(100), Cigar::Match(100)]);
        test_record.set(
            "test_reads".as_bytes(),
            Some(&test_cigar),
            "AAA".as_bytes(),
            "~~~".as_bytes(),
        );
        let test_range = Range64::new(350, 360);
        let deletion_size = 50;
        let deletion_present = check_del(&test_record, deletion_size, test_range);
        assert!(!deletion_present);
    }
}
