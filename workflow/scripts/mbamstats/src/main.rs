use core::str;
use std::{
    collections::HashMap, io::{stdout, BufWriter, Write}, path::Path
};

use itertools::Itertools;
use noodles::{
    bam,
    sam::alignment::record::data::field::{value::Array, Value},
};
use rayon::iter::{ParallelBridge, ParallelIterator};

struct ReadModProb {
    read_name: String,
    // Total positions associated with modification.
    // ex. A => AT
    seq_content: usize,
    /// Modification
    mtype: String,
    /// Probabilities
    probs: Vec<f32>,
}

impl ReadModProb {
    pub fn valid_proportion(&self) -> f32 {
        self.probs.len() as f32 / self.seq_content as f32
    }
}

#[inline(always)]
pub fn rc_nt(nt: char) -> char {
    match nt {
        'A' => 'T',
        'T' => 'A',
        'G' => 'C',
        'C' => 'G',
        _ => unreachable!("Invalid base"),
    }
}

fn collect_read_mod_stats(bamfile: impl AsRef<Path>, mtypes: &HashMap<&str, f32>) -> Vec<ReadModProb> {
    let mut bam_reader = bam::io::reader::Builder::default()
        .build_from_path(bamfile)
        .unwrap();
    let _header = bam_reader.read_header().unwrap();

    bam_reader
        .records()
        .par_bridge()
        .flatten()
        .flat_map(|rec| {
            let mut indices: HashMap<&str, (usize, usize)> = HashMap::new();

            let tags = rec.data();
            let (Some(Ok(Value::String(mm))), Some(Ok(Value::Array(Array::UInt8(ml))))) =
                (tags.get(b"MM"), tags.get(b"ML"))
            else {
                eprintln!("ML not array of u8 or MM not a string.");
                return None;
            };

            let mut offset = 0;
            for j in str::from_utf8(mm)
                .unwrap()
                .split(';')
                .filter(|mtype| !mtype.is_empty())
            {
                // mtype,positions
                let Some((mtype, positions)) = j.split_once(',') else {
                    continue;
                };
                if !mtypes.contains_key(mtype) {
                    continue;
                }
                // TODO: Map positions back to query sequence.
                let n_positions = positions.split(',').count();
                indices.insert(mtype, (offset, offset + n_positions));
                offset += n_positions;
            }

            let seq = rec.sequence();
            let seq_content: HashMap<char, usize> = seq.iter().fold(HashMap::new(), |mut acc, nt| {
                acc.entry(nt as char)
                    .and_modify(|cnt| *cnt += 1)
                    .or_insert(1);
                acc
            });

            let read_name = rec.name().unwrap_or_default();

            Some(
                indices
                    .into_iter()
                    .flat_map(|(mtype, (st_idx, end_idx))| {
                        let Some(nt) = mtype
                            .split('+')
                            .collect_tuple()
                            .and_then(|(nt, _)| nt.chars().next())
                        else {
                            eprintln!("{mtype} can't split.");
                            return None;
                        };
                        let seq_content = seq_content[&nt]
                            + seq_content.get(&rc_nt(nt)).cloned().unwrap_or_default();
                        let thr_mtype = mtypes.get(mtype).cloned().unwrap_or_default();
                        let mtype = mtype.to_string();
                        let probs = ml
                            .iter()
                            .get(st_idx..end_idx)
                            .flatten()
                            .filter_map(|prob| {
                                let new_prob = prob as f32 / 256.0;
                                (new_prob > thr_mtype).then_some(new_prob)
                            })
                            .collect_vec();

                        Some(ReadModProb {
                            seq_content,
                            mtype,
                            probs,
                            read_name: read_name.to_string(),
                        })
                    })
                    .collect_vec(),
            )
        })
        .flatten()
        .collect()
}

fn calculate_summary_stats(probs: &[f32]) -> (f32, f32, f32, f32) {
    let n =  probs.len() as f32;
    let mean = probs.iter().sum::<f32>() / n;
    let var = probs.iter().map(|x| (x - mean).powi(2)).sum::<f32>() / n;
    let stdev = var.sqrt();
    (mean, var, stdev, n)
}

fn main() {
    let args = std::env::args().collect_vec();

    // /project/logsdon_shared/projects/HGSVC3/NA20355_Dimelo/data/original_run/CENPA/20250610_DMLseq_GM20355_ULK114.bam
    let infile = args.get(1).expect("Need bam file.");

    let thr_mtypes: HashMap<&'static str, f32> = HashMap::from_iter([("A+a.", 0.8)]);
    let records = collect_read_mod_stats(infile, &thr_mtypes);

    let mut writer = BufWriter::new(stdout());

    // Calculate proportion of valid A+a. sites
    let all_probs = records.iter().map(|rec| rec.valid_proportion()).collect_vec();

    let (global_mean_prob, _global_var_prob, global_stdev_prob, _global_n) = calculate_summary_stats(all_probs.as_slice());

    for rec in records
    {
        // Calculate 1 stdev, 2 stdev, 3 stdev, mean.
        let valid_prop = rec.valid_proportion();
        let zscore = (valid_prop - global_mean_prob) / global_stdev_prob;
        writeln!(&mut writer, "{}\t{}\t{valid_prop}\t{zscore}", rec.read_name, rec.mtype).unwrap();
    }
}
