use crate::phaser::{Exception, Phaser};
use crate::toolkit::range::I64 as Range64;
use crate::toolkit::site_selection::{
    pileup as site_pileup, CandidateSite, FilteredSites, RawVariantCounts,
};
use crate::toolkit::util::DError;
use vstr::VString;

/// Takes a slice and returns a slice for the remainder of the string following all numeric
/// characters.
#[must_use]
fn skip_digits(x: &[u8]) -> &[u8] {
    let seq_start = x
        .iter()
        .position(|&x| !x.is_ascii_digit())
        .unwrap_or(x.len());
    &x[seq_start..]
}

/// Takes a deletion length + count field, extracts the digit-only portion, and then converts it
/// to an integer.
/// "1N"
fn extract_del_length(x: &[u8]) -> Result<i32, DError> {
    if x.is_empty() {
        return Err(Exception::new("extract_del_length got empty input").into());
    }
    let non_digit_idx = x
        .iter()
        .position(|&x| !x.is_ascii_digit())
        .unwrap_or(x.len());
    if non_digit_idx == 0 {
        return Err(Exception::new(format!(
            "extract_del_length expected leading digits, found: {}",
            vstr::VStr::from(x)
        ))
        .into());
    }
    let mut out = 0i32;
    for c in &x[..non_digit_idx] {
        out = out
            .checked_mul(10)
            .and_then(|acc| acc.checked_add(i32::from(*c - b'0')))
            .ok_or_else(|| Exception::new("extract_del_length overflowed i32"))?;
    }
    Ok(out)
}

impl Phaser {
    /// Add homozygous sites to fingerprinting.
    #[must_use]
    pub(crate) fn add_hom_sites(
        &mut self,
        min_no_var_region_size: Option<i64>,
        max_hom_var_to_add: Option<usize>,
        ref_seq: &[u8],
    ) -> Vec<CandidateSite> {
        let offset = match self.try_offset() {
            Ok(offset) => offset,
            Err(e) => {
                log::warn!(
                    "Failed to compute coordinate offset while adding homozygous anchor sites: {e}; skipping anchor-site augmentation."
                );
                return Vec::new();
            }
        };
        let min_no_var_region_size = min_no_var_region_size.unwrap_or(10000);
        let max_hom_var_to_add = max_hom_var_to_add.unwrap_or(10);
        let mut ret = Vec::<CandidateSite>::new();
        self.het_sites.sort_by_key(|a| a.pos);
        self.hom_sites.sort_by_key(|a| a.pos);
        let het_pos = self
            .het_sites
            .iter()
            .filter(|x| x.pos < self.right_boundary_0based() && x.pos > self.left_boundary_0based())
            .map(|x| x.pos)
            .collect::<Vec<_>>();
        if het_pos.is_empty() {
            self.hom_sites
                .iter()
                .filter(|site| site.is_unit_length())
                .for_each(|x| ret.push(x.clone()));
            if self.hom_sites.is_empty() {
                let full_range = self.right_boundary_0based() - self.left_boundary_0based();
                let interval_size = full_range / 4;
                let positions = vec![
                    self.left_boundary_0based() + interval_size,
                    self.left_boundary_0based() + interval_size * 2,
                    self.left_boundary_0based() + interval_size * 3,
                ];
                for var_pos in positions {
                    let pos_on_ref = var_pos - offset;
                    let ref_base_u8 = ref_seq[pos_on_ref as usize];
                    let non_ref_bases = [b'A', b'C', b'G', b'T']
                        .iter()
                        .filter(|x| **x != ref_base_u8)
                        .copied()
                        .collect::<Vec<_>>();
                    let var_base_u8 = non_ref_bases.first().copied().unwrap_or(b'N');
                    let ref_base = char::from(ref_base_u8).to_string();
                    let var_base = char::from(var_base_u8).to_string();
                    let new_variant = CandidateSite::new(var_pos, ref_base, var_base);
                    ret.push(new_variant);
                }
            }
        } else {
            let (Some(min_pos), Some(max_pos)) = (het_pos.iter().min(), het_pos.iter().max())
            else {
                return ret;
            };
            if *min_pos - self.left_boundary_0based() > min_no_var_region_size {
                self.hom_sites
                    .iter()
                    .filter(|site| site.pos < *min_pos && site.is_unit_length())
                    .for_each(|x| ret.push(x.clone()));
            }
            if self.right_boundary_0based() - *max_pos > min_no_var_region_size {
                self.hom_sites
                    .iter()
                    .filter(|site| site.pos > *max_pos && site.is_unit_length())
                    .for_each(|x| ret.push(x.clone()));
            }
            let het_sites_no_del = self
                .het_sites
                .iter()
                .filter(|x| x.is_unit_length())
                .collect::<Vec<_>>();
            let het_site_num = het_sites_no_del.len();
            if het_site_num > 1 {
                for i in 0..(het_site_num - 1) {
                    let interval_start = het_sites_no_del[i].pos;
                    let interval_end = het_sites_no_del[i + 1].pos;
                    if interval_end - interval_start > min_no_var_region_size {
                        self.hom_sites
                            .iter()
                            .filter(|site| {
                                site.pos < interval_end
                                    && site.pos > interval_start
                                    && site.is_unit_length()
                            })
                            .for_each(|x| ret.push(x.clone()));
                    }
                }
            }
        }
        if !ret.is_empty() {
            ret.sort();
            let num_sites = ret.len();
            for hom_site in ret.iter().step_by(num_sites.div_ceil(max_hom_var_to_add)) {
                log::debug!("Adding homozygous anchor site into phasing set: {hom_site:?}");
                self.het_sites.push((*hom_site).clone());
            }
            self.het_sites.sort();
        }
        ret
    }

    /// Remove variant sites within any regions marked as `noisy_regions`.
    pub(crate) fn remove_noisy_sites(&mut self) {
        self.het_sites.retain(|site| {
            !self
                .noisy_regions
                .iter()
                .any(|region| region.contains(&site.pos))
        });
    }

    /// Generate candidate sites from pileup across the locus interval.
    ///
    /// Applies configured site-selection settings plus optional `min_vaf` override.
    pub fn get_candidate_pos(
        &mut self,
        regions_to_check: &[Range64],
        seq: &[u8],
        min_vaf: Option<f64>,
    ) -> Result<(FilteredSites, RawVariantCounts), DError> {
        let mut bam_handle = self.try_realigned_bam()?;

        let mut this_setting = self.settings.site_selection_settings.clone();
        if let Some(min_vaf_value) = min_vaf {
            this_setting.min_vaf = min_vaf_value;
        }
        if let Some(user_min_vaf) = self.settings.min_variant_frequency {
            this_setting.min_vaf = user_min_vaf;
        }
        this_setting.targeted = self.settings.targeted;

        let (filtered_variants, raw_variants) = site_pileup(
            &mut bam_handle,
            &Range64::new(self.left_boundary_0based(), self.right_boundary_0based()),
            seq,
            self,
            Some(&this_setting),
            regions_to_check,
            &self.local_chr().ok_or_else(|| {
                Exception::new(format!(
                    "Missing chromosome in realign region '{}'",
                    self.realign_region
                ))
            })?,
        )?;

        Ok((filtered_variants, raw_variants))
    }

    /// Whether deletion-support markers are allowed at this position.
    ///
    /// This is true when the position falls within a partially supported
    /// configured deletion interval.
    #[must_use]
    pub fn allow_del_bases(&self, pos: i64) -> bool {
        self.del_data.iter().any(|deletion| {
            !deletion.del_reads_partial.is_empty()
                && deletion.threep().start <= pos
                && pos <= deletion.fivep().end
        })
    }

    /// Compute normalized REF/ALT strings and size for one pileup indel token.
    pub fn process_indel(
        &self,
        pos: i64,
        ref_seq: &[u8],
        var_seq: &[u8],
        cached_faidx: &[u8],
    ) -> Result<(VString, VString, i32), DError> {
        Self::free_process_indel(pos, ref_seq, var_seq, cached_faidx)
    }

    /// Stateless helper for indel token normalization.
    ///
    /// Supports pileup insertion forms (`A+3CCC`) and deletion forms (`A-3NNN`).
    pub fn free_process_indel(
        pos: i64,
        ref_seq: &[u8],
        var_seq: &[u8],
        cached_faidx: &[u8],
    ) -> Result<(VString, VString, i32), DError> {
        let indel_size: i32;
        let mut var_ret = VString::default();
        let mut ref_ret = VString::default();
        if let Some(plus_index) = var_seq.iter().position(|&x| x == b'+') {
            let insertion = var_seq.split_at(plus_index + 1).1;
            let insertion = skip_digits(insertion);
            indel_size = i32::try_from(insertion.len())?;
            var_ret.extend_from_slice(ref_seq);
            var_ret.extend_from_slice(insertion);
            ref_ret.extend_from_slice(ref_seq);
        } else {
            var_ret.extend_from_slice(ref_seq);
            let minus_index = if let Some(index) = var_seq.iter().position(|&x| x == b'-') {
                Result::<usize, simple_error::SimpleError>::Ok(index)
            } else {
                simple_error::bail!(
                    "- and + not found in indel seqs. var seq: {}",
                    vstr::VStr::from(var_seq)
                )
            }?;
            let deletion_len = extract_del_length(var_seq.split_at(minus_index + 1).1)?;
            indel_size = deletion_len;
            let deletion_len = deletion_len as usize;
            let pos = pos as usize;
            let end = pos
                .checked_add(deletion_len + 1)
                .ok_or_else(|| Exception::new("process_indel deletion length overflow"))?;
            if end > cached_faidx.len() {
                return Err(Exception::new(format!(
                    "process_indel deletion range {pos}..{end} exceeds reference length {}",
                    cached_faidx.len()
                ))
                .into());
            }
            let cached_seq = &cached_faidx[pos..end];
            debug_assert_eq!(cached_seq, cached_seq.to_ascii_uppercase());
            debug_assert_eq!(cached_seq.len(), deletion_len + 1);
            ref_ret.extend_from_slice(cached_seq);
        }
        ref_ret.make_ascii_uppercase();
        var_ret.make_ascii_uppercase();
        Ok((ref_ret, var_ret, indel_size))
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn skip_digits_ok() {
        use super::skip_digits;
        assert_eq!(skip_digits(b"ACGT"), b"ACGT");
        assert_eq!(skip_digits(b"123ACGT"), b"ACGT");
        assert_eq!(skip_digits(b"123"), &[] as &[u8]);
    }

    #[test]
    fn extract_del_len_ok() {
        use super::extract_del_length;
        assert_eq!(extract_del_length(b"1N").unwrap(), 1);
        assert_eq!(extract_del_length(b"100N").unwrap(), 100);
        assert_eq!(extract_del_length(b"37FD34").unwrap(), 37);
        assert!(extract_del_length(b"").is_err());
        assert!(extract_del_length(b"N10").is_err());
    }

    #[test]
    fn process_indel_ok() {
        use crate::phaser::Phaser;
        use crate::toolkit::util;
        let tmp = tempfile::tempdir().unwrap();
        let local_fa = tmp.path().join("smn1_ref.fa");
        std::fs::copy(util::test_file("ref/smn1_ref.fa"), &local_fa).unwrap();
        crate::phaser::build_faidx(&local_fa).unwrap();
        let faidx = rust_htslib::faidx::Reader::from_path(&local_fa).unwrap();
        let all_ref_seq = util::load_all_seqs(&faidx);
        let offset = 70_890_000;
        // both offset and pos are 1-based, which is okay.

        let (ref_seq, var_seq, indel_len) = Phaser::free_process_indel(
            70_940_935 - offset,
            b"A",
            b"A+3CCC",
            //70_889_999,88
            &all_ref_seq[0],
        )
        .unwrap();
        assert_eq!(*ref_seq, &b"A"[..]);
        assert_eq!(*var_seq, &b"ACCC"[..]);
        assert_eq!(indel_len, 3, "Insertion has wrong length.");

        let (ref_seq, var_seq, indel_len) = Phaser::free_process_indel(
            70_940_935 - offset,
            b"A",
            b"A-2NN",
            //70_889_999,
            &all_ref_seq[0],
        )
        .unwrap();
        assert_eq!(*ref_seq, &b"ACT"[..]);
        assert_eq!(*var_seq, &b"A"[..]);
        assert_eq!(indel_len, 2, "Deletion has wrong length");
    }
}
