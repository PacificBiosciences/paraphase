use crate::io::json::GeneCall;
use crate::phaser::Phaser;
use crate::realign::reference_length;
use crate::toolkit::hapcmp::HapCompare;
use crate::toolkit::util::{output_bam_header, DError, DResult};

use vstr::VStr;

use rand::{Rng, SeedableRng};
use rust_htslib::bam::{self, Read};

use std::path::{Path, PathBuf};

pub mod colors {
    pub const READ: &str = "166,206,227";
    pub const READ_ALLELE1: &str = "178,223,138";
    pub const READ_ALLELE2: &str = "177,156,217";
}

pub struct BamWriter<'a> {
    phaser: &'a Phaser,
    call: &'a GeneCall,
}

pub struct IOTuple(pub PathBuf, pub PathBuf, pub String, pub bool);

impl IOTuple {
    #[must_use]
    pub fn source_bam(&self) -> &Path {
        &self.0
    }
    #[must_use]
    pub fn dest_bam(&self) -> &Path {
        &self.1
    }
    #[must_use]
    pub fn chromosome_name(&self) -> &str {
        &self.2
    }
    #[must_use]
    pub fn is_gene2(&self) -> bool {
        self.3
    }
}

impl std::fmt::Debug for IOTuple {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "IOTuple{{Source: {:?}. Dest: {:?}. Chrom: {:?}. Primary or secondary gene: {}",
            self.source_bam(),
            self.dest_bam(),
            self.chromosome_name(),
            if self.is_gene2() {
                "secondary"
            } else {
                "primary"
            }
        )
    }
}

impl<'a> BamWriter<'a> {
    #[must_use]
    /// Build a BAM writer bound to one phasing run and its gene call outputs.
    pub fn new(phaser: &'a Phaser, call: &'a GeneCall) -> Self {
        Self { phaser, call }
    }

    /// Attach phasing tags (`RN`, `HP`, optional color tags) to one read record.
    ///
    /// This prefers unique-support assignments, then falls back to heuristic
    /// matching or seeded random assignment for non-unique supporting reads.
    fn add_tag_to_read(
        &self,
        record: &mut bam::Record,
        use_supp: bool,
        is_gene2: bool,
        rng: Option<&mut rand::rngs::SmallRng>,
    ) -> DResult {
        record.push_aux(b"RN", bam::record::Aux::String(self.phaser.gene_name()))?;
        let haps = &self.call.final_haplotypes;
        if haps.is_empty() {
            record.push_aux(b"HP", bam::record::Aux::String("Unassigned"))?;
            return Ok(());
        }
        let nonunique_reads = &self.call.nonunique_supporting_reads;

        let mut alleles = Vec::new();
        let alleles_in_call = self.call.region_specific_info.get("raw_alleles");
        if let Some(raw_alleles) = alleles_in_call.and_then(|x| x.as_array()) {
            for allele in raw_alleles {
                if let Some(allele_parsed) = allele.as_array() {
                    let allele_string = allele_parsed
                        .iter()
                        .filter_map(|x| x.as_str().map(str::to_string))
                        .collect::<Vec<_>>();
                    if !allele_string.is_empty() {
                        alleles.push(allele_string);
                    }
                }
            }
        }

        let read_details = &self.call.read_details;
        let qname = VStr::from(record.qname()).to_string();
        let qname = if use_supp && record.is_supplementary() && !is_gene2 {
            let ref_length = reference_length(record);
            let ref_start = record.pos();
            format!("{qname}_sup_{ref_start}_{ref_length}")
        } else {
            qname
        };
        for (hap, hap_name) in &self.call.final_haplotypes {
            if let Some(reads) = self.call.unique_supporting_reads.get(hap) {
                if reads.contains(&qname) {
                    record.push_aux(b"HP", bam::record::Aux::String(hap_name))?;
                    let color = bam::record::Aux::String(if alleles.is_empty() {
                        colors::READ
                    } else if alleles[0].contains(hap_name) {
                        colors::READ_ALLELE1
                    } else if alleles.len() > 1 && alleles[1].contains(hap_name) {
                        colors::READ_ALLELE2
                    } else {
                        colors::READ
                    });
                    record.push_aux(b"YC", color)?;
                    return Ok(());
                }
            }
        }
        if let Some(fingerprint) = read_details.get(&qname) {
            let mut mismatches = Vec::with_capacity(self.call.final_haplotypes.len());
            for (result, hap) in self
                .call
                .final_haplotypes
                .keys()
                .map(|hap| (HapCompare::from_haps(fingerprint, hap), hap))
            {
                let result = result?;
                mismatches.push((hap, result.mismatches));
            }
            mismatches.sort_by(|a, b| a.1.cmp(&b.1));
            if mismatches.len() > 1
                && (1..=2).contains(&mismatches[0].1)
                && mismatches[1].1 >= mismatches[0].1 + 2
            {
                let hp_match =
                    self.call
                        .final_haplotypes
                        .get(mismatches[0].0)
                        .ok_or_else(|| {
                            crate::phaser::Exception::new(format!(
                                "Failed to resolve haplotype match key '{}' from final_haplotypes",
                                mismatches[0].0
                            ))
                        })?;
                record.push_aux(b"HP", bam::record::Aux::String(hp_match))?;
                let color = bam::record::Aux::String(if alleles.is_empty() {
                    colors::READ
                } else if alleles[0].contains(hp_match) {
                    colors::READ_ALLELE1
                } else if alleles.len() > 1 && alleles[1].contains(hp_match) {
                    colors::READ_ALLELE2
                } else {
                    colors::READ
                });
                record.push_aux(b"YC", color)?;
                return Ok(());
            }
        }
        /*
        let assignment = if nonunique_reads.contains_key(&qname) {
            let mut possible_assignments = nonunique_reads.get(&qname).unwrap().clone();
            possible_assignments.sort();
            let first_assignment = possible_assignments.first().unwrap();
            &first_assignment.to_string()
        } else {
            "Unassigned"
        };
        */
        let assignment = rng
            .and_then(|rng| {
                nonunique_reads.get(&qname).map(|possible| {
                    let mut possible_assignments = possible.clone();
                    possible_assignments.sort();
                    if possible_assignments.is_empty() {
                        "Unassigned".to_string()
                    } else {
                        let random_index = rng.gen_range(0..possible_assignments.len());
                        possible_assignments[random_index].clone()
                    }
                })
            })
            .unwrap_or_else(|| String::from("Unassigned"));

        if let Some(assignment1) = self.call.final_haplotypes.get(&assignment) {
            record.push_aux(b"HP", bam::record::Aux::String(assignment1))?;
        } else {
            record.push_aux(b"HP", bam::record::Aux::String(&assignment))?;
        }
        Ok(())
    }

    /// Write tagged BAM outputs for gene1 and, when configured, gene2.
    ///
    /// Returns the list of written BAM paths and builds corresponding BAI files.
    pub fn write_bams(&self) -> Result<Vec<PathBuf>, DError> {
        let mut written_bams = Vec::new();
        let gene1_input_bam = self.phaser.realigned_bam_path();
        let gene1_output_bam = self.phaser.realigned_tagged_bam_path();
        let chr = self.phaser.chr().ok_or_else(|| {
            crate::phaser::Exception::new(format!(
                "Failed to determine chromosome for primary region of gene '{}'",
                self.phaser.gene_name()
            ))
        })?;
        let gene1_inputs = IOTuple(gene1_input_bam, gene1_output_bam, chr.to_string(), false);
        // Assign non-uniquely supporting reads randomly to one haplotype
        self.write_bam(gene1_inputs, Some(0))?;
        written_bams.push(self.phaser.realigned_tagged_bam_path());
        bam::index::build(
            self.phaser.realigned_tagged_bam_path(),
            None,
            bam::index::Type::Bai,
            1,
        )?;

        let region = self
            .phaser
            .locus_config()
            .gene2_region(self.phaser.settings.genome == "37");
        // write to gene2
        if let Some(gene2_region) = region {
            let gene2_chr = gene2_region
                .split_terminator(':')
                .next()
                .ok_or_else(|| {
                    crate::phaser::Exception::new(format!(
                        "Malformed gene2 region '{}': missing chromosome",
                        gene2_region
                    ))
                })?
                .to_string();
            let gene2_output_bam = self.phaser.realigned_tagged_gene2_bam_path();
            let gene2_input_bam = self.phaser.realigned_gene2_bam_path();
            let gene2_inputs = IOTuple(
                gene2_input_bam,
                gene2_output_bam,
                gene2_chr,
                /* is_gene_2= */ true,
            );
            self.write_bam(gene2_inputs, Some(0))?;
            written_bams.push(self.phaser.realigned_tagged_gene2_bam_path());
            bam::index::build(
                self.phaser.realigned_tagged_gene2_bam_path(),
                None,
                bam::index::Type::Bai,
                1,
            )?;
        }
        Ok(written_bams)
    }

    /// Attempts to write a bam
    /// # Inputs
    /// 1. `IOTuple` - inbam, outbam, isgene2
    /// 2. Seed - `Option<usize>`. If None, random assign is disabled. Otherwise, seeds a RNG.
    /// # Return
    /// 1. Path to output bam.
    /// 2. bool - whether or not written. If no haplotypes found, no bam is written.
    pub fn write_bam(&self, tuple: IOTuple, seed: Option<u64>) -> DResult {
        let mut rng = seed.map(rand::rngs::SmallRng::seed_from_u64);
        let use_supp = self.phaser.use_supplementary();
        log::debug!("Resolved IO tuple: {tuple:?}");
        let mut reader = bam::IndexedReader::from_path(tuple.source_bam())?;
        let mut tmp_bam_writer = bam::Writer::from_path(
            tuple.dest_bam(),
            &output_bam_header(reader.header()),
            bam::Format::Bam,
        )?;

        reader.fetch(tuple.chromosome_name())?;
        let mut record = bam::Record::new();
        while let Some(rc) = reader.read(&mut record) {
            if rc.is_err() {
                continue;
            }
            if record.is_secondary() {
                continue;
            }
            self.add_tag_to_read(&mut record, use_supp, tuple.is_gene2(), rng.as_mut())?;
            tmp_bam_writer.write(&record)?;
        }
        Ok(())
    }
}
