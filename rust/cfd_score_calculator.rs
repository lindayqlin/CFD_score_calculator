//! # Cutting frequency determination (CFD) score calculator
//! Written by Linda Lin 3/23/2025 (adapted from python)
//!
//! Command line usage:
//!     ./cfd_score_calculator <20nt_aligned_spacer_sequence> <20nt_aligned_protospacer_sequence>
//!     <last_2nt_of_PAM>
//! Example:
//!     ./cfd_score_calculator CTAACAGTTGCTTTTATCAC tT-ACAGcTGCaTTTATCAC GG
//!
//! Module usage:
//!     use cfd_score_calculator::calculate_cfd;
//!
//! Input: mismatch_scores.txt & pam_scores.txt
//! Output: CFD score

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};

/// Calculate CFD score
pub fn calculate_cfd(spacer: &str, protospacer: &str, pam: &str) -> f32 {
    // Check for expected input lengths
    if spacer.len() != 20 || protospacer.len() != 20 || pam.len() != 2 {
        eprintln!("Incorrect input sequence length\n\
                   Usage: cfd_score_calculator <20nt_aligned_spacer_sequence> \
                   <20nt_aligned_protospacer_sequence> <last_2nt_of_PAM>");
        std::process::exit(1);
    }

    // Get empirical scoring matrices and pre-process sequences
    let mm_scores = parse_scoring_matrix("mismatch_scores.txt");
    let pam_scores = parse_scoring_matrix("pam_scores.txt");
    let spacer_list: Vec<char> = spacer.to_uppercase().replace('T', "U").chars().collect();
    let protospacer_list: Vec<char> = protospacer.to_uppercase().replace('T', "U")
                                                 .chars().collect();

    // Calculate CFD score for alignment by nucleotide
    let mut score = 1.0;
    for (i, &nt) in protospacer_list.iter().enumerate() {
        if spacer_list[i] == nt {
            // No penalty for perfect match
            score *= 1.0;
        } else if i == 0 && (spacer_list[i] == '-' || nt == '-'){
            // No penalty for gap at most PAM-distal nucleotide
            // Conservative choice to maximize sensitivity given lack of empirical data
            score *= 1.0;
        }
        else {
            // Incorporate score for given RNA-DNA basepair at this position
            let key = format!("r{}:d{},{}",
                              spacer_list[i], reverse_complement_nt(&nt.to_string()), i + 1);
            score *= mm_scores.get(&key).expect("Invalid basepair");
        }
    }
    // Incorporate PAM score
    score *= pam_scores.get(&pam.to_uppercase()).expect("Invalid PAM");
    score
}

/// Parse scoring matrix from space-delimited file
fn parse_scoring_matrix(file_path: &str) -> HashMap<String, f32> {
    // Open file
    let file = match File::open(file_path) {
        Ok(file) => file,
        Err(_) => {
            eprintln!("Requires {} in directory", file_path);
            std::process::exit(1);
        }
    };
    // Read file
    let reader = BufReader::new(file);
    let mut matrix = HashMap::new();
    for line in reader.lines() {
        let line = line.unwrap();
        let parts: Vec<&str> = line.split_whitespace().collect();
        matrix.insert(parts[0].to_string(), parts[1].parse::<f32>().unwrap());
    }
    matrix
}

/// Get reverse complement of nucleotide (supports bulges)
fn reverse_complement_nt(seq: &str) -> String {
    seq.chars().map(|x| match x {
        'A' => 'T',
        'C' => 'G',
        'T' | 'U' => 'A',
        'G' => 'C',
        '-' => '-',
        _ => x,
    }).collect::<String>()
}

fn main() {
    // Get args
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 4 {
        eprintln!("Too few arguments\nUsage: cfd_score_calculator <20nt_aligned_spacer_sequence> <20nt_aligned_protospacer_sequence> <last_2nt_of_PAM>");
        std::process::exit(1);
    }
    let spacer = &args[1];
    let protospacer = &args[2];
    let pam = &args[3];

    // Calculate CFD score
    println!("{}", calculate_cfd(spacer, protospacer, pam));
}