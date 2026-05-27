use crate::io::json::{ReadAlignmentId, ReadFingerprintMap};
use crate::phaser::{Exception, Phaser};
use crate::toolkit::deletion::Datum as DeletionDatum;
use crate::toolkit::site_selection::CandidateSite;
use crate::toolkit::util::DError;

use std::collections::{BTreeMap, BTreeSet};

const CLIP_OFFSET_HAPS_FROM_READS: i64 = 30;

/// Internal implementation for variant-site filtering from read fingerprints.
///
/// Updates `phaser.het_sites` in-place and returns per-site absent-base hints
/// used downstream for homozygous augmentation logic.
fn remove_var_impl(
    phaser: &mut Phaser,
    raw_read_haps: &ReadFingerprintMap,
    kept_sites: &[CandidateSite],
    hom_sites: &[CandidateSite],
) -> BTreeMap<CandidateSite, u8> {
    let mut absent_base_per_site = BTreeMap::new();
    if phaser.het_sites.is_empty() {
        log::debug!("No heterozygous variant sites available; skipping remove_var.");
        return absent_base_per_site;
    }

    assert_eq!(
        raw_read_haps
            .values()
            .map(|x| x.len())
            .collect::<BTreeSet<_>>()
            .len(),
        1,
        "All read fingerprints should be of uniform length."
    );

    let hap_len = raw_read_haps.values().next().map_or(0, |x| x.len());
    let mut bases_per_site = vec![counter::Counter::<u8, i32>::new(); hap_len];
    log::debug!(
        "Initial read fingerprints before variant filtering: count={}",
        raw_read_haps.len()
    );
    for hap in raw_read_haps.values() {
        assert_eq!(hap.len(), bases_per_site.len());
        for (base, counter) in hap.iter().copied().zip(bases_per_site.iter_mut()) {
            *counter.entry(base).or_default() += 1;
        }
    }

    log::trace!("Counts for each variant site: {bases_per_site:?}");

    let mut sites_to_remove = Vec::new();
    for (pos, base_counter) in bases_per_site
        .iter()
        .enumerate()
        .filter(|x| !x.1.is_empty())
    {
        let ref_count = base_counter.get(&b'1').copied().unwrap_or(0);
        let alt_count = base_counter.get(&b'2').copied().unwrap_or(0);
        let base0_count = base_counter.get(&b'0').copied().unwrap_or(0);
        let x_count = base_counter.get(&b'x').copied().unwrap_or(0);
        let ref_plus_alt = ref_count + alt_count;
        let total = base_counter.total::<i32>();
        let in_kept_sites = kept_sites.contains(&phaser.het_sites[pos]);
        let mut to_remove = false;
        if x_count == total - base0_count {
            to_remove = true;
        } else if ref_plus_alt == (total - x_count - base0_count)
            && (alt_count <= 3 || ref_count <= 3)
        {
            if alt_count <= 3 && ref_count <= 3 {
                to_remove = true;
            } else {
                if !in_kept_sites {
                    to_remove = true;
                }
                if hom_sites.contains(&phaser.het_sites[pos]) {
                    if alt_count <= 3 {
                        absent_base_per_site.insert(phaser.het_sites[pos].clone(), b'2');
                    } else if ref_count <= 3 {
                        absent_base_per_site.insert(phaser.het_sites[pos].clone(), b'1');
                    }
                }
            }
        }
        if to_remove {
            log::debug!(
                "Filtering heterozygous site index={pos}: x_count={x_count}, ref_count={ref_count}, alt_count={alt_count}, ref_plus_alt={ref_plus_alt}, kept_site={in_kept_sites}"
            );
            sites_to_remove.push(pos);
        }
    }
    log::debug!("Filtered heterozygous site indices: {sites_to_remove:?}");
    for idx in &sites_to_remove {
        let var_to_remove = &phaser.het_sites[*idx];
        if phaser.init_het_sites.contains(var_to_remove) {
            if let Some(idx_init) = phaser
                .init_het_sites
                .iter()
                .position(|x| x == var_to_remove)
            {
                phaser.init_het_sites.swap_remove(idx_init);
            }
        }
    }
    sites_to_remove.into_iter().rev().for_each(|idx| {
        phaser.het_sites.swap_remove(idx);
    });
    absent_base_per_site
}

/// Apply one deletion event update to haplotype fingerprints and site list.
///
/// Inserts/decorates deletion marker symbols and updates deletion category map.
fn sub_update_for_deletions(
    phaser: &mut Phaser,
    hap_map: &mut ReadFingerprintMap,
    data: &DeletionDatum,
    index: usize,
    categories: &mut BTreeMap<char, String>,
) -> Result<(), DError> {
    let initial_het_len = phaser.het_sites.len();
    for read in data.del_reads_partial.iter() {
        if hap_map.get(&(read.clone().into())).is_none() {
            hap_map.insert(read.into(), vec![b'x'; initial_het_len].into());
        }
    }
    for read in data.del_negative_reads.iter() {
        if hap_map.get(&(read.clone().into())).is_none() {
            hap_map.insert(read.into(), vec![b'x'; initial_het_len].into());
        }
    }

    let mut het_sites = std::mem::take(&mut phaser.het_sites);
    let signifier = match index {
        0 => '3',
        1 => '4',
        _ => ((65 + index) as u8) as char,
    };
    let base = match u8::try_from(signifier) {
        Ok(v) => v,
        Err(_) => return Err("char out of u8 bounds. Too many deletions?".into()),
    };

    let del_range_start = data.threep().start;
    let del_range_end = data.fivep().end;
    let pos1 = het_sites.iter().position(|x| x.pos > del_range_start);
    let pos2 = het_sites.iter().position(|x| x.pos > del_range_end);
    if let (Some(pos1), Some(pos2)) = (pos1, pos2) {
        let range = pos1..pos2;
        match pos1.cmp(&pos2) {
            std::cmp::Ordering::Less => {
                for hap in data.del_reads_partial.iter().map(ReadAlignmentId::from) {
                    if let Some(hap) = hap_map.get_mut(&hap) {
                        hap[range.clone()].fill(base);
                    }
                }
            }
            std::cmp::Ordering::Equal => {
                let cand = CandidateSite::from(data);
                het_sites.insert(pos1, cand);
                if pos1 == 0 {
                    for (name, hap) in hap_map.iter_mut() {
                        hap.insert(pos1, b'x');
                        if data.del_reads_partial.contains(&name.read_name) {
                            hap[pos1] = base;
                        } else if pos1 + 1 < hap.len() && hap.get(pos1 + 1) == Some(&b'0') {
                            hap[pos1] = b'0';
                        } else if data.del_negative_reads.contains(&name.read_name) {
                            hap[pos1] = b'1';
                        }
                    }
                } else {
                    for (name, hap) in hap_map.iter_mut() {
                        hap.insert(pos1, b'x');
                        if data.del_reads_partial.contains(&name.read_name) {
                            hap[pos1] = base;
                        } else if pos1 >= 1
                            && hap.get(pos1 - 1) == Some(&b'0')
                            && pos1 + 1 < hap.len()
                            && hap.get(pos1 + 1) == Some(&b'0')
                        {
                            hap[pos1] = b'0';
                        } else {
                            let hap_len = hap.len();
                            let flanking_left_start = pos1.saturating_sub(2);
                            let flanking_left_has_x =
                                hap[flanking_left_start..pos1].contains(&b'x');
                            let flanking_right_has_x = hap[std::cmp::min(pos1 + 1, hap_len)
                                ..std::cmp::min(pos1 + 3, hap_len)]
                                .contains(&b'x');
                            if !flanking_left_has_x && !flanking_right_has_x {
                                hap[pos1] = b'1';
                            }
                        }
                    }
                }
            }
            std::cmp::Ordering::Greater => {}
        }
    } else if pos1.is_some() && pos2.is_none() {
        let nvar = het_sites.len();
        let Some(pos1) = pos1 else {
            phaser.het_sites = het_sites;
            return Ok(());
        };
        let range = pos1..nvar;
        for hap in data.del_reads_partial.iter().map(ReadAlignmentId::from) {
            if let Some(hap) = hap_map.get_mut(&hap) {
                hap[range.clone()].fill(base);
            }
        }
    } else if pos1.is_none() && pos2.is_none() {
        let cand = CandidateSite::from(data);
        het_sites.push(cand);
        let pos1 = het_sites.len() - 1;
        for (name, hap) in hap_map.iter_mut() {
            hap.push(b'x');
            if data.del_reads_partial.contains(&name.read_name) {
                hap[pos1] = base;
            } else if pos1 >= 1 && hap.get(pos1 - 1) == Some(&b'0') {
                hap[pos1] = b'0';
            } else if data.del_negative_reads.contains(&name.read_name) {
                hap[pos1] = b'1';
            }
        }
    }

    categories.insert(signifier, data.name());
    phaser.het_sites = het_sites;
    log::debug!(
        "Deletion update adjusted heterozygous-site count: before={initial_het_len}, after={}",
        phaser.het_sites.len()
    );
    Ok(())
}

/// Internal implementation that applies all configured deletion updates.
fn update_for_deletions_impl(
    phaser: &mut Phaser,
    hap_map: &mut ReadFingerprintMap,
) -> Result<BTreeMap<char, String>, DError> {
    let mut ret = BTreeMap::new();
    for (deletion_index, deletion) in phaser.del_data.clone().into_iter().enumerate() {
        if phaser.del_data[deletion_index].del_reads_partial.len() > 1 {
            sub_update_for_deletions(phaser, hap_map, &deletion, deletion_index, &mut ret)?;
        }
    }

    Ok(ret)
}

use crate::realign::reference_length;
use crate::toolkit::range::I64 as Range64;
use crate::toolkit::util::{DResult, HashSet};

use vstr::VStr;

use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::{self, pileup, Read};

use std::ops::Index;

/// Create a unique identifier for the read alignment.
/// This lets us distinguish multiple alignments per input read.
///
/// This is a free function, so we can call it without borrowing the `Phaser` struct.
pub(super) fn get_read_name_free(record: &bam::Record, use_supplementary: bool) -> String {
    let qname = VStr::from(record.qname());
    if use_supplementary && record.is_supplementary() {
        let ref_start = record.pos();
        let ref_len = reference_length(record);
        format!("{qname}_sup_{ref_start}_{ref_len}")
    } else {
        qname.to_string()
    }
}

/// Create a unique identifier for the read alignment.
/// This lets us distinguish multiple alignments per input read.
///
/// This requires borrowing `Phaser`.
/// To avoid borrowing, use `get_read_name_free`.
pub(super) fn get_read_name(phaser: &Phaser, record: &bam::Record) -> String {
    get_read_name_free(record, phaser.use_supplementary())
}

/// Yields a vector of names. First the normal read name, and,
/// # Panics
/// Panics if read name is not utf-encoded.
pub(super) fn get_read_names(
    phaser: &Phaser,
    record: &bam::Record,
    partial_deletion_reads: Option<&BTreeSet<String>>,
) -> Vec<String> {
    let qname = VStr::from(record.qname()).to_string();
    let mut ret = vec![phaser.get_read_name(record)];
    if phaser.locus_config().use_supplementary() && record.is_supplementary() {
        if let Some(partial_deletion_reads) = partial_deletion_reads {
            if partial_deletion_reads.contains(&qname) && partial_deletion_reads.contains(&ret[0]) {
                ret.push(qname);
            }
        }
    }
    ret
}

/// Create a unique identifier for an alignment.
/// This lets us disambiguate primary and supplementary alignments.
#[must_use]
pub(super) fn labeled_read_name(record: &bam::Record) -> String {
    let qname = vstr::VStr::from(record.qname());
    if record.is_supplementary() {
        let pos = record.pos();
        let reference_end = record.reference_end();
        let length = reference_end - pos;
        format!("{qname}_sup_{pos}_{length}")
    } else {
        qname.to_string()
    }
}

/// Compute the length of a cigar operation.
#[inline]
/// Return clipping length for one CIGAR op (soft/hard clips only).
fn clip_size(x: bam::record::Cigar) -> u32 {
    match x {
        bam::record::Cigar::SoftClip(size) | bam::record::Cigar::HardClip(size) => size,
        _ => 0,
    }
}

/// 3' clip length.
/// Only checks last cigar operation.
/// Return 3' clipping length from a CIGAR vector.
///
/// Returns `0` when the terminal operation is not a clip.
pub fn threeprime_clip_length(x: &[bam::record::Cigar]) -> u32 {
    x.last().copied().map_or(0u32, clip_size)
}

/// 5' clip length.
/// Only checks first cigar operation.
/// Return 5' clipping length from a CIGAR vector.
///
/// Returns `0` when the leading operation is not a clip.
pub fn fiveprime_clip_length(x: &[bam::record::Cigar]) -> u32 {
    x.first().copied().map_or(0u32, clip_size)
}

/// Checks if a known deletion is present in the read.
/// # Arguments
/// * `record` - bam record
/// * `size` - deletion size
/// * `del_threeprime_range` - coordinate range for left of the deletion (reads with three prime clips)
/// # Returns
/// * a bool indicating the presence/absence of the deletion in the record
/// Check whether a read supports a target deletion by overlap and size match.
pub fn check_del(record: &bam::Record, size: i64, del_threeprime_range: Range64) -> bool {
    let mut starting_pos = record.pos();
    for cigar in record.cigar().iter() {
        let cigar_len = i64::from(cigar.len());
        let cigar_char = cigar.char();
        if cigar_char == 'D' {
            let len_diff = (cigar_len - size).abs() as f64;
            let diff_cutoff = (size as f64 * 0.1).min(50.0);
            if len_diff < diff_cutoff {
                let padding = size / 10;
                if starting_pos >= del_threeprime_range.start - padding
                    && starting_pos <= del_threeprime_range.end + padding
                {
                    return true;
                }
            }
        }
        if cigar_char == 'M' || cigar_char == 'D' || cigar_char == '=' || cigar_char == 'X' {
            starting_pos += cigar_len;
        }
    }
    false
}

/// Handle clipped positions within `haplotypes_from_reads_step`.
///
/// # Errors
/// Propagates `DError` from bam open/fetch/read operations.
/// Update clip-based hap marker symbols around configured clip positions.
fn handle_clip_step(
    phaser: &mut Phaser,
    read_haps: &mut ReadFingerprintMap,
    min_clip_len: u32,
    het_sites: &[CandidateSite],
    tid: i32,
    clip_buffer: Option<i32>,
) -> DResult {
    log::debug!("Adding clips to read abstractions");
    let mut reader = phaser.try_realigned_bam()?;
    let mut record = bam::Record::new();
    let nvar = phaser.het_sites.len();
    let clip_buffer_value: i64 = clip_buffer.unwrap_or(20).into();

    for (fingerprint_index, allele_site) in het_sites.iter().enumerate() {
        for clip_position in phaser
            .clip_3p_positions
            .iter()
            .copied()
            .filter(|&clip_position| allele_site.pos > clip_position)
        {
            reader.fetch((
                tid,
                std::cmp::max(0, clip_position - clip_buffer_value),
                clip_position + clip_buffer_value,
            ))?;
            while let Some(code) = reader.read(&mut record) {
                if let Err(e) = code {
                    log::warn!(
                        "Failed to read record while applying 3' clip mask; skipping record: {e:?}"
                    );
                    continue;
                }
                let read_name = phaser.get_read_name(&record);
                let entry = read_haps
                    .entry(read_name.into())
                    .or_insert_with(|| vec![b'x'; nvar].into());
                if (clip_position + 1 - record.reference_end()).abs() < clip_buffer_value {
                    let cigar = record.cigar();
                    let threep_clip_len = threeprime_clip_length(&cigar);
                    if threep_clip_len >= min_clip_len {
                        entry[fingerprint_index] = b'0';
                    }
                }
            }
        }
    }

    for (fingerprint_index, allele_site) in het_sites.iter().enumerate() {
        for clip_position in phaser
            .clip_5p_positions
            .iter()
            .rev()
            .copied()
            .filter(|&clip_position| allele_site.pos < clip_position)
        {
            reader.fetch((
                tid,
                std::cmp::max(0, clip_position - clip_buffer_value),
                clip_position + clip_buffer_value,
            ))?;
            while let Some(code) = reader.read(&mut record) {
                if let Err(e) = code {
                    log::warn!(
                        "Failed to read record while applying 5' clip mask; skipping record: {e:?}"
                    );
                    continue;
                }
                let entry = read_haps
                    .entry(phaser.get_read_name(&record).into())
                    .or_insert_with(|| vec![b'x'; nvar].into());
                if (clip_position + 1 - record.reference_start()).abs() < clip_buffer_value {
                    let cigar = record.cigar();
                    if fiveprime_clip_length(&cigar) >= min_clip_len {
                        entry[fingerprint_index] = b'0';
                    }
                }
            }
        }
    }
    Ok(())
}

/// For each variant site, performs updates.
///
/// # Errors
/// 1. bam index query failure.
/// 2. pileup iteration failure. (htslib error propagation.)
#[allow(clippy::too_many_arguments)]
/// Update one read's fingerprint symbol at one heterozygous-site index.
fn update_fingerprint_map(
    phaser: &mut Phaser,
    read_haps: &mut ReadFingerprintMap,
    flanking_indel_reads: &mut HashSet<String>,
    exclude_reads: Option<&BTreeSet<String>>,
    partial_deletion_reads: Option<&BTreeSet<String>>,
    site_data: (&CandidateSite, usize, u8),
    tid: i32,
    absent_base_per_site: Option<&BTreeMap<CandidateSite, u8>>,
) -> DResult {
    let (allele_site, fingerprint_index, min_mapq) = site_data;
    log::trace!(
        "Fetching reads from aligned-locus BAM for fingerprint index {fingerprint_index}: bam={:?}",
        phaser.realigned_bam_path(),
    );
    let mut bam = phaser.try_realigned_bam()?;
    bam.fetch((tid, allele_site.pos - 1, allele_site.pos + 1))?;

    let empty_map = BTreeMap::new();
    let absent_base_per_site = absent_base_per_site.unwrap_or(&empty_map);

    let mut read_name_inclusion = BTreeMap::<String, i32>::new();
    let mut read_name_exclusion = BTreeMap::<String, i32>::new();

    let nvar = phaser.het_sites.len();
    let mut num_pileups = 0usize;

    let mut in_exclude_count = 0usize;
    let mut in_flanking_indel_count = 0usize;
    let mut num_included = 0usize;

    flanking_indel_reads.clear();
    for p in bam.pileup() {
        let p = p?;
        num_pileups += 1;
        let diff = allele_site.pos - i64::from(p.pos());
        match diff {
            1 => {
                log::trace!("Inspecting upstream pileup position before target site: allele_site={allele_site}, pileup_pos={}", p.pos());
                for name in p
                    .alignments()
                    .filter(|x| {
                        (x.indel() != pileup::Indel::None || x.is_del())
                            && super::base_qual(x)
                                >= phaser.settings.site_selection_settings.min_base_quality
                    })
                    .flat_map(|aln| phaser.get_read_names(&aln.record(), partial_deletion_reads))
                {
                    log::trace!(
                        "Marking read as flanking-indel support: read={name}, allele_pos={}",
                        allele_site.pos
                    );
                    flanking_indel_reads.insert(name);
                }
            }
            0 => {
                log::trace!(
                    "Processing target site pileup: allele_site={allele_site}, flanking_indel_reads={flanking_indel_reads:?}"
                );
                for pileup_read in p.alignments() {
                    let base_qual = super::base_qual(&pileup_read);
                    let record = pileup_read.record();
                    let read_is_included = !pileup_read.is_del()
                        && !pileup_read.is_refskip()
                        && !record.is_secondary()
                        && record.mapq() >= min_mapq
                        && pileup_read.indel() == pileup::Indel::None
                        && base_qual >= phaser.settings.site_selection_settings.min_base_quality;
                    if read_is_included {
                        /*
                        log::trace!(
                            "Read included: {}. Pos: {}. Is del: {}. is refskip: {}. is secondary: {}. mapq {} vs min {min_mapq}, and indel {:?}. Base qual {base_qual} vs {}",
                            VStr::from(record.qname()),
                            p.pos(),
                            pileup_read.is_del(),
                            pileup_read.is_refskip(),
                            record.is_secondary(),
                            record.mapq(),
                            pileup_read.indel(),
                            phaser.settings.site_selection_settings.min_base_quality
                        );
                        */
                        let names = phaser.get_read_names(&record, partial_deletion_reads);
                        let read_name = VStr::from(record.qname());
                        for name in names {
                            let in_exclude =
                                exclude_reads.is_some_and(|exclude| exclude.contains(&name));
                            let in_flanking_indel = flanking_indel_reads.contains(&name);
                            if !in_exclude && !in_flanking_indel {
                                //log::trace!(
                                //    "Read passes exclusion/flanking-indel filters: read={read_name}"
                                //);
                                *read_name_inclusion.entry(name.clone()).or_default() += 1;
                                num_included += 1;
                                let Some(qpos) = pileup_read.qpos() else {
                                    continue;
                                };
                                let threep_clip_len =
                                    threeprime_clip_length(&record.cigar()) as usize;
                                if qpos < record.seq_len() - threep_clip_len - 1 {
                                    let prohibited_base = absent_base_per_site.get(allele_site);
                                    let base1 = prohibited_base != Some(&b'1');
                                    let base2 = prohibited_base != Some(&b'2');
                                    let entry = read_haps
                                        .entry((&name).into())
                                        .or_insert_with(|| vec![b'x'; nvar].into());
                                    let base = record.seq().index(qpos).to_ascii_uppercase();
                                    let base = if base1
                                        && (allele_site.ref_seq == base
                                            || phaser
                                                .settings
                                                .site_selection_settings
                                                .permit_list
                                                .get(&allele_site.pos)
                                                .is_some_and(|seq| seq.vstr() == base))
                                    {
                                        Some(b'1')
                                    } else if allele_site.var_seq == base && base2 {
                                        Some(b'2')
                                    } else {
                                        None
                                    };
                                    if let Some(base) = base {
                                        entry[fingerprint_index] = base;
                                        log::trace!(
                                            "Updated read fingerprint base call: read={name}, index={fingerprint_index}, encoded_base={base}, char_base={}",
                                            base as char
                                        );
                                    }
                                } else {
                                    let seq_len = record.seq_len();
                                    assert!(
                                        qpos < seq_len,
                                        "qpos out of bounds... qpos >= record seq len?? {qpos}, {seq_len}",
                                    );
                                }
                            } else {
                                log::trace!(
                                    "Read excluded from fingerprint update: read={read_name}, in_exclude={in_exclude}, in_flanking_indel={in_flanking_indel}"
                                );
                                *read_name_exclusion.entry(name.clone()).or_default() += 1;
                                if in_exclude {
                                    in_exclude_count += 1;
                                }
                                if in_flanking_indel {
                                    in_flanking_indel_count += 1;
                                }
                            }
                        }
                    } else {
                        for name in phaser.get_read_names(&record, partial_deletion_reads) {
                            *read_name_exclusion.entry(name.clone()).or_default() += 1;
                            log::trace!(
                                "Read filtered out at target pileup: read={}, pos={}, is_del={}, is_refskip={}, is_secondary={}, mapq={} (min={min_mapq}), indel={:?}, base_qual={} (min={})",
                                VStr::from(record.qname()),
                                p.pos(),
                                pileup_read.is_del(),
                                pileup_read.is_refskip(),
                                record.is_secondary(),
                                record.mapq(),
                                pileup_read.indel(),
                                base_qual,
                                phaser.settings.site_selection_settings.min_base_quality
                            );
                        }
                    }
                }
                break;
            }
            _ => {
                assert!(
                    diff >= 0,
                    "Expected loop termination upon reaching variant site. diff: {diff}"
                );
                continue;
            }
        }
    }

    log::trace!("Processed pileup rows for site update: num_pileups={num_pileups}");
    log::trace!(
        "Fingerprint update inclusion summary: included={num_included}, excluded_by_list={in_exclude_count}, excluded_by_flanking_indel={in_flanking_indel_count}"
    );
    log::trace!(
        "Per-read inclusion/exclusion counters: included={read_name_inclusion:?}, excluded={read_name_exclusion:?}"
    );
    Ok(())
}

/// Build read haplotypes by scanning reads at each selected variant site.
///
/// # Errors
/// 1. Updating fingerprint map.
/// 2. Failure to handle clips.
#[allow(clippy::too_many_arguments)]
/// Internal step: project pileup observations into read fingerprint map.
fn haplotypes_from_reads_step_impl(
    phaser: &mut Phaser,
    exclude_reads: Option<&BTreeSet<String>>,
    min_clip_len: u32,
    check_clip: bool,
    partial_deletion_reads: Option<&BTreeSet<String>>,
    min_mapq: u8,
    tid: i32,
    clip_buffer: Option<i32>,
    absent_base_per_site: Option<&BTreeMap<CandidateSite, u8>>,
) -> Result<ReadFingerprintMap, DError> {
    log::debug!(
        "Starting haplotypes_from_reads step: gene={}, sample={}, het_site_count={}",
        phaser.gene_name(),
        phaser.sample_id(),
        phaser.het_sites.len()
    );
    let mut ret = ReadFingerprintMap::default();
    let mut reads_with_flanking_indels = HashSet::<String>::default();
    let het_sites = phaser.het_sites.clone();
    for (fingerprint_index, allele_site) in het_sites.iter().enumerate() {
        log::trace!(
            "Updating read fingerprint map for heterozygous-site index {fingerprint_index}"
        );
        update_fingerprint_map(
            phaser,
            &mut ret,
            &mut reads_with_flanking_indels,
            exclude_reads,
            partial_deletion_reads,
            (allele_site, fingerprint_index, min_mapq),
            tid,
            absent_base_per_site,
        )?;
    }

    if check_clip
        && std::cmp::max(
            phaser.clip_3p_positions.len(),
            phaser.clip_5p_positions.len(),
        ) > 0
    {
        handle_clip_step(phaser, &mut ret, min_clip_len, &het_sites, tid, clip_buffer)?;
    }

    Ok(ret)
}

impl Phaser {
    /// Build read fingerprints across selected sites, with optional clip-derived
    /// synthetic markers and absent-base suppression.
    #[allow(clippy::too_many_arguments)]
    pub fn haplotypes_from_reads(
        &mut self,
        exclude_reads: Option<&BTreeSet<String>>,
        kept_sites: &[CandidateSite],
        add_sites: Option<&[CandidateSite]>,
        partial_deletion_reads: Option<&BTreeSet<String>>,
        options: (u8, bool, Option<u32>),
        tid: i32,
        clip_buffer: Option<i32>,
        hom_sites: &[CandidateSite],
    ) -> Result<ReadFingerprintMap, DError> {
        const CLIP_OFFSET: i64 = CLIP_OFFSET_HAPS_FROM_READS;
        const ACGT: &[u8] = b"ACGT";
        let offset = self.try_offset()?;

        let (min_mapq, check_clip, min_clip_len) = options;
        let add_sites = add_sites.unwrap_or(&self.add_sites[..]).to_vec();
        let min_clip_len = min_clip_len.unwrap_or(50);
        log::debug!(
            "Building read fingerprints: min_mapq={min_mapq}, check_clip={check_clip}, min_clip_len={min_clip_len}"
        );
        let mut raw_read_haps = self.haplotypes_from_reads_step(
            exclude_reads,
            min_clip_len,
            check_clip,
            partial_deletion_reads,
            min_mapq,
            tid,
            clip_buffer,
            None,
        )?;
        log::debug!(
            "Initial read fingerprints after first pass: count={}",
            raw_read_haps.len()
        );
        log::debug!(
            "Filtering low-support variant sites from fingerprints: fingerprint_len={}, kept_site_count={}",
            raw_read_haps.values().next().map_or(0, |x| x.len()),
            kept_sites.len()
        );
        let absent_base_per_site = self.remove_var(&raw_read_haps, kept_sites, hom_sites);

        if !self.het_sites.is_empty() {
            self.het_sites.extend_from_slice(&add_sites);
            self.het_sites.sort();
            self.het_sites.dedup();
        }

        let faidx = self.make_local_faidx()?;
        let chr = self.local_chr().ok_or(Exception::new("chr must be set"))?;
        log::debug!("Using local reference contig for synthetic clip-site anchors: {chr}");
        let seq = faidx
            .fetch_seq(&chr, 0, i64::MAX as usize)?
            .to_ascii_uppercase();
        debug_assert!((0..seq.len()).all(|idx| seq[idx]
            == faidx
                .fetch_seq(&chr, idx, idx + 1)
                .ok()
                .and_then(|x| x.first().copied())
                .map_or(seq[idx], |b| b.to_ascii_uppercase())));
        drop(faidx);

        let num_clip = self.clip_5p_positions.len();
        for i in 0..num_clip {
            let clip_pos = self.clip_5p_positions[i];
            let has_var_before_clip = if i == 0 {
                self.het_sites.iter().any(|x| x.pos < clip_pos)
            } else {
                self.het_sites
                    .iter()
                    .any(|x| x.pos > self.clip_5p_positions[i - 1] && x.pos < clip_pos)
            };
            if !has_var_before_clip {
                let var_pos = clip_pos - CLIP_OFFSET;
                let ref_base = seq[(var_pos - offset) as usize];
                let var_base =
                    ACGT.iter()
                        .copied()
                        .find(|&x| x != ref_base)
                        .ok_or(Exception::new(format!(
                        "Expected a variant base at 5' {clip_pos} position with var_pos = {var_pos}"
                    )))?;
                let new_var = CandidateSite::new(var_pos, vec![ref_base], vec![var_base]);
                self.het_sites.push(new_var);
            }
        }

        let num_clip = self.clip_3p_positions.len();
        for i in (0..num_clip).rev() {
            let clip_pos = self.clip_3p_positions[i];
            let has_var_after_clip = if i == num_clip - 1 {
                self.het_sites.iter().any(|x| x.pos > clip_pos)
            } else {
                self.het_sites
                    .iter()
                    .any(|x| x.pos < self.clip_3p_positions[i + 1] && x.pos > clip_pos)
            };
            if !has_var_after_clip {
                let var_pos = clip_pos + CLIP_OFFSET;
                let ref_base = seq[(var_pos - offset) as usize];
                let var_base =
                    ACGT.iter()
                        .copied()
                        .find(|&x| x != ref_base)
                        .ok_or(Exception::new(format!(
                        "Expected a variant base at 3' {clip_pos} position with var_pos = {var_pos}"
                    )))?;
                let new_var = CandidateSite::new(var_pos, vec![ref_base], vec![var_base]);
                self.het_sites.push(new_var);
            }
        }

        self.het_sites.sort();
        self.het_sites.dedup();

        raw_read_haps = self.haplotypes_from_reads_step(
            exclude_reads,
            min_clip_len,
            check_clip,
            partial_deletion_reads,
            min_mapq,
            tid,
            clip_buffer,
            Some(&absent_base_per_site),
        )?;
        Ok(raw_read_haps)
    }

    /// Remove low-support variant sites and return per-site absent-base hints used
    /// in the second fingerprinting pass.
    pub fn remove_var(
        &mut self,
        raw_read_haps: &ReadFingerprintMap,
        kept_sites: &[CandidateSite],
        hom_sites: &[CandidateSite],
    ) -> BTreeMap<CandidateSite, u8> {
        remove_var_impl(self, raw_read_haps, kept_sites, hom_sites)
    }

    /// Inject deletion markers into read fingerprints for configured/observed large deletions.
    pub fn update_for_deletions(
        &mut self,
        hap_map: &mut ReadFingerprintMap,
    ) -> Result<BTreeMap<char, String>, DError> {
        update_for_deletions_impl(self, hap_map)
    }

    /// Create a unique identifier for the read alignment.
    pub(crate) fn get_read_name_free(record: &bam::Record, use_supplementary: bool) -> String {
        get_read_name_free(record, use_supplementary)
    }

    /// Create a unique identifier for the read alignment.
    pub(crate) fn get_read_name(&self, record: &bam::Record) -> String {
        get_read_name(self, record)
    }

    /// Yields names for this read, including supplementary aliasing when configured.
    pub fn get_read_names(
        &self,
        record: &bam::Record,
        partial_deletion_reads: Option<&BTreeSet<String>>,
    ) -> Vec<String> {
        get_read_names(self, record, partial_deletion_reads)
    }

    /// Internal fingerprinting pass over current `het_sites`.
    ///
    /// Optionally applies clip-position updates and absent-base exclusions.
    #[allow(clippy::too_many_arguments)]
    pub fn haplotypes_from_reads_step(
        &mut self,
        exclude_reads: Option<&BTreeSet<String>>,
        min_clip_len: u32,
        check_clip: bool,
        partial_deletion_reads: Option<&BTreeSet<String>>,
        min_mapq: u8,
        tid: i32,
        clip_buffer: Option<i32>,
        absent_base_per_site: Option<&BTreeMap<CandidateSite, u8>>,
    ) -> Result<ReadFingerprintMap, DError> {
        haplotypes_from_reads_step_impl(
            self,
            exclude_reads,
            min_clip_len,
            check_clip,
            partial_deletion_reads,
            min_mapq,
            tid,
            clip_buffer,
            absent_base_per_site,
        )
    }

    /// Create a unique identifier for an alignment.
    #[must_use]
    pub fn labeled_read_name(record: &bam::Record) -> String {
        labeled_read_name(record)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config;
    use crate::phaser;
    use crate::phaser::Phaser;
    use crate::toolkit::deletion::Datum as DeletionDatum;
    use crate::toolkit::site_selection::CandidateSite;
    use crate::toolkit::util;
    use std::str::FromStr;
    use vstr::VString;

    fn build_smn1_phaser() -> Option<Phaser> {
        let outdir = tempfile::TempDir::new().expect("tempdir should build");
        let Some(genome_path) = std::env::var("HG38")
            .ok()
            .map(|x| x.trim_end_matches(".mmi").to_string())
        else {
            log::warn!("Skipping read_abstraction test because the HG38 environment variable is not configured.");
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

    fn stringify_sites(sites: &[CandidateSite]) -> Vec<String> {
        sites.iter().map(std::string::ToString::to_string).collect()
    }

    #[test]
    fn update_for_deletions_matches_python_cases() {
        let Some(mut phaser) = build_smn1_phaser() else {
            return;
        };
        let base_sites = vec![
            CandidateSite::from_str("70917101_A_C").unwrap(),
            CandidateSite::from_str("70917111_A_C").unwrap(),
            CandidateSite::from_str("70917150_A_C").unwrap(),
            CandidateSite::from_str("70917200_A_C").unwrap(),
            CandidateSite::from_str("70917300_A_C").unwrap(),
            CandidateSite::from_str("70917400_A_C").unwrap(),
        ];

        phaser.het_sites = base_sites.clone();
        let mut mid_del = DeletionDatum::new(
            Range64::new(70917198, 70917200),
            None,
            Some(70917198),
            Some(70917199),
            Some(70917200),
            Some(70917200),
        );
        mid_del.del_reads_partial.insert(String::from("r1"));
        mid_del.del_reads_partial.insert(String::from("r2"));
        phaser.del_data = vec![mid_del];
        let mut hap_map =
            ReadFingerprintMap::from([(ReadAlignmentId::from_name("r1"), VString::from("111111"))]);
        phaser
            .update_for_deletions(&mut hap_map)
            .expect("middle deletion update should succeed");
        assert_eq!(
            hap_map[&ReadAlignmentId::from_name("r1")].to_string(),
            "111311"
        );
        assert_eq!(
            stringify_sites(&phaser.het_sites),
            vec![
                String::from("70917101_A_C"),
                String::from("70917111_A_C"),
                String::from("70917150_A_C"),
                String::from("70917200_A_C"),
                String::from("70917300_A_C"),
                String::from("70917400_A_C"),
            ]
        );

        let Some(mut phaser) = build_smn1_phaser() else {
            return;
        };
        phaser.het_sites = base_sites.clone();
        let mut start_del = DeletionDatum::new(
            Range64::new(70917099, 70917099),
            None,
            Some(70917099),
            Some(70917099),
            Some(70917099),
            Some(70917099),
        );
        start_del.del_reads_partial.insert(String::from("r1"));
        start_del.del_reads_partial.insert(String::from("r2"));
        phaser.del_data = vec![start_del];
        let mut hap_map =
            ReadFingerprintMap::from([(ReadAlignmentId::from_name("r1"), VString::from("111111"))]);
        phaser
            .update_for_deletions(&mut hap_map)
            .expect("leading deletion update should succeed");
        assert_eq!(
            hap_map[&ReadAlignmentId::from_name("r1")].to_string(),
            "3111111"
        );
        assert_eq!(
            stringify_sites(&phaser.het_sites),
            vec![
                String::from("70917100_del_0"),
                String::from("70917101_A_C"),
                String::from("70917111_A_C"),
                String::from("70917150_A_C"),
                String::from("70917200_A_C"),
                String::from("70917300_A_C"),
                String::from("70917400_A_C"),
            ]
        );

        let Some(mut phaser) = build_smn1_phaser() else {
            return;
        };
        phaser.het_sites = base_sites;
        let mut end_del = DeletionDatum::new(
            Range64::new(70917499, 70917599),
            None,
            Some(70917499),
            Some(70917520),
            Some(70917580),
            Some(70917599),
        );
        end_del.del_reads_partial.insert(String::from("r1"));
        end_del.del_reads_partial.insert(String::from("r2"));
        phaser.del_data = vec![end_del];
        let mut hap_map =
            ReadFingerprintMap::from([(ReadAlignmentId::from_name("r1"), VString::from("111111"))]);
        phaser
            .update_for_deletions(&mut hap_map)
            .expect("trailing deletion update should succeed");
        assert_eq!(
            hap_map[&ReadAlignmentId::from_name("r1")].to_string(),
            "1111113"
        );
        assert_eq!(
            stringify_sites(&phaser.het_sites),
            vec![
                String::from("70917101_A_C"),
                String::from("70917111_A_C"),
                String::from("70917150_A_C"),
                String::from("70917200_A_C"),
                String::from("70917300_A_C"),
                String::from("70917400_A_C"),
                String::from("70917500_del_100"),
            ]
        );
    }
}
