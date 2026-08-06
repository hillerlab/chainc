// Copyright (c) 2026 The Hiller Lab at the Senckenberg Gesellschaft für Naturforschung
// Distributed under the terms of the Apache License, Version 2.0.

use std::path::PathBuf;

use clap::Parser;

#[derive(Debug, Clone, Parser)]
#[command(
    name = "chainc",
    version,
    about = "Remove chain-breaking alignments identified by a hierarchical NET"
)]
pub struct CleanerConfig {
    /// Input UCSC chain file.
    #[arg(long, short = 'c', value_name = "CHAIN")]
    pub chains: PathBuf,
    /// Hierarchical NET. `None` nets the input chains in memory.
    #[arg(long, short = 'n', value_name = "NET")]
    pub net: Option<PathBuf>,
    /// Ordered reference sizes, used when no NET is given.
    ///
    /// Defaults to the reference sequence file's own order. Generated NET
    /// sections follow this order, exactly as `chainNet` follows its tSizes.
    #[arg(long, short = 'R', value_name = "CHROM_SIZES")]
    pub reference_sizes: Option<PathBuf>,
    /// Ordered query sizes, used when no NET is given.
    #[arg(long, short = 'Q', value_name = "CHROM_SIZES")]
    pub query_sizes: Option<PathBuf>,
    /// Write the NET used for discovery, mainly to inspect a generated one.
    #[arg(long, short = 'N', value_name = "NET", hide = true)]
    pub dump_net: Option<PathBuf>,
    /// Reference genome in .2bit or FASTA format.
    #[arg(long, short = 'r', value_name = "SEQUENCE")]
    pub reference: PathBuf,
    /// Query genome in .2bit or FASTA format.
    #[arg(long, short = 'q', value_name = "SEQUENCE")]
    pub query: PathBuf,
    /// Sorted output chain file.
    #[arg(long, short = 'o', value_name = "CHAIN")]
    pub output: PathBuf,
    /// BED9 output describing removed suspects.
    #[arg(long, short = 'b', value_name = "BED")]
    pub removed_bed: PathBuf,
    /// Linear gap model: loose, medium, or a Kent-style gap file.
    #[arg(long, short = 'g', value_name = "MODEL")]
    pub linear_gap: String,
    /// Blastz/LASTZ-format substitution score matrix.
    #[arg(long, short = 's', value_name = "FILE")]
    pub score_scheme: Option<PathBuf>,
    #[arg(long, short = 'F', default_value_t = 0.0)]
    pub fold_threshold: f64,
    #[arg(long, short = 'L', default_value_t = 2.5)]
    pub lr_fold_threshold: f64,
    #[arg(long, short = 'P', default_value_t = 10.0)]
    pub lr_fold_threshold_pairs: f64,
    #[arg(long, short = 'B', default_value_t = i32::MAX as f64)]
    pub max_suspect_bases: f64,
    #[arg(long, short = 'M', default_value_t = 100_000.0)]
    pub max_suspect_score: f64,
    #[arg(long, short = 'S', default_value_t = 50_000.0)]
    pub min_broken_chain_score: f64,
    #[arg(long, short = 'G', default_value_t = 0)]
    pub min_lr_gap_size: u32,
    /// Also test adjacent pairs of suspects.
    #[arg(long, short = 'p')]
    pub do_pairs: bool,
    #[arg(long, short = 'D', default_value_t = 10_000)]
    pub max_pair_distance: u32,
    /// Write NEW_CHAIN_ID<TAB>ORIGINAL_BREAKING_CHAIN_ID.
    #[arg(long, short = 'i', value_name = "FILE")]
    pub new_chain_id_dict: Option<PathBuf>,
    /// Dry-run table of every single suspect; disables removal and pairs.
    #[arg(long, short = 'd', value_name = "FILE")]
    pub suspect_data_file: Option<PathBuf>,
    /// Number of worker threads used by parallel parsers.
    #[arg(long, short = 't')]
    pub threads: Option<usize>,
    /// Reject ranking anomalies and zero-local-score candidates.
    #[arg(long, short = 'e')]
    pub strict: bool,
    /// Write the compatibility-ordered candidate trace.
    #[arg(long, short = 'x', value_name = "FILE", hide = true)]
    pub dump_candidate_order: Option<PathBuf>,
}
