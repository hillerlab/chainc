// Copyright (c) 2026 The Hiller Lab at the Senckenberg Gesellschaft für Naturforschung
// Distributed under the terms of the Apache License, Version 2.0.

use std::fs;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

use chainc::config::CleanerConfig;

static NEXT_ID: AtomicU64 = AtomicU64::new(0);

fn temporary_directory() -> PathBuf {
    let id = NEXT_ID.fetch_add(1, Ordering::Relaxed);
    let path = std::env::temp_dir().join(format!("chainc-test-{}-{id}", std::process::id()));
    fs::create_dir_all(&path).unwrap();
    path
}

fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures")
        .join(name)
}

/// Permissive configuration: every threshold is disabled so a fixture's
/// geometry alone decides what happens.
fn config(stem: &str, directory: &Path) -> CleanerConfig {
    CleanerConfig {
        chains: fixture(&format!("{stem}.chain")),
        net: Some(fixture(&format!("{stem}.net"))),
        reference_sizes: None,
        query_sizes: None,
        dump_net: None,
        threads: None,
        reference: fixture("reference.fa"),
        query: fixture("query.fa"),
        output: directory.join("clean.chain"),
        removed_bed: directory.join("removed.bed"),
        linear_gap: "loose".into(),
        score_scheme: None,
        fold_threshold: 0.0,
        lr_fold_threshold: 0.0,
        lr_fold_threshold_pairs: 0.0,
        max_suspect_bases: f64::MAX,
        max_suspect_score: f64::MAX,
        min_broken_chain_score: 0.0,
        min_lr_gap_size: 0,
        do_pairs: false,
        max_pair_distance: 10_000,
        new_chain_id_dict: None,
        suspect_data_file: None,
        dump_candidate_order: None,
        strict: false,
    }
}

#[test]
fn removes_one_suspect_and_rescores_output() {
    let directory = temporary_directory();
    let mut config = config("single", &directory);
    config.new_chain_id_dict = Some(directory.join("ids.tsv"));
    config.dump_candidate_order = Some(directory.join("candidates.tsv"));
    let output = config.output.clone();
    let bed = config.removed_bed.clone();

    chainc::run(config).unwrap();

    let output_text = fs::read_to_string(output).unwrap();
    assert!(output_text.starts_with("#fixture=single-removal\n"));
    assert!(output_text.contains(" 3\n10\n\n"));
    assert!(output_text.contains(" 1\n10\t80\t80\n10\n\n"));
    let bed_text = fs::read_to_string(bed).unwrap();
    assert!(bed_text.starts_with("chr1\t45\t55\tbreakingChainID_1_"));
    assert!(bed_text.ends_with("\t0,0,153\n"));
    assert_eq!(
        fs::read_to_string(directory.join("ids.tsv")).unwrap(),
        "3\t1\n"
    );
    assert_eq!(
        fs::read_to_string(directory.join("candidates.tsv"))
            .unwrap()
            .lines()
            .count(),
        1
    );

    fs::remove_dir_all(directory).unwrap();
}

#[test]
fn suspect_data_mode_is_a_dry_run() {
    let directory = temporary_directory();
    let mut config = config("single", &directory);
    config.do_pairs = true;
    config.suspect_data_file = Some(directory.join("suspects.tsv"));
    let output = config.output.clone();
    let bed = config.removed_bed.clone();

    chainc::run(config).unwrap();
    assert_eq!(fs::read_to_string(bed).unwrap(), "");
    assert_eq!(
        fs::read_to_string(directory.join("suspects.tsv"))
            .unwrap()
            .lines()
            .count(),
        1
    );
    let output_text = fs::read_to_string(output).unwrap();
    assert!(!output_text.contains(" 3\n"));
    assert!(output_text.contains(" 1\n10\t35\t35\n10\t35\t35\n10\n\n"));

    fs::remove_dir_all(directory).unwrap();
}

#[test]
fn removes_adjacent_suspects_as_a_pair() {
    let directory = temporary_directory();
    let mut config = config("pair", &directory);
    config.reference = fixture("pair_reference.fa");
    config.query = fixture("pair_query.fa");
    config.lr_fold_threshold = 3.0;
    config.lr_fold_threshold_pairs = 2.5;
    config.do_pairs = true;
    config.max_pair_distance = 30;
    let output = config.output.clone();
    let bed = config.removed_bed.clone();

    chainc::run(config).unwrap();
    let bed_text = fs::read_to_string(bed).unwrap();
    assert!(bed_text.starts_with("chr1\t45\t95\t"));
    assert!(bed_text.ends_with("\t0,100,255\n"));
    let output_text = fs::read_to_string(output).unwrap();
    assert!(output_text.contains(" 3\n10\t30\t30\n10\n\n"));
    assert!(output_text.contains(" 1\n10\t120\t120\n10\n\n"));

    fs::remove_dir_all(directory).unwrap();
}

#[test]
fn output_is_identical_across_thread_counts() {
    let root = temporary_directory();
    let run_with_threads = |threads: usize, label: &str| {
        let directory = root.join(label);
        fs::create_dir_all(&directory).unwrap();
        let mut config = config("single", &directory);
        config.new_chain_id_dict = Some(directory.join("ids.tsv"));
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap()
            .install(|| chainc::run(config))
            .unwrap();
        directory
    };

    let one = run_with_threads(1, "one");
    let four = run_with_threads(4, "four");
    for file in ["clean.chain", "removed.bed", "ids.tsv"] {
        assert_eq!(
            fs::read(one.join(file)).unwrap(),
            fs::read(four.join(file)).unwrap(),
            "{file} differs by thread count"
        );
    }

    fs::remove_dir_all(root).unwrap();
}

/// Without `--net`, the chains are netted in memory by the ported
/// `chainNet -minScore=0` plus `NetFilterNonNested.perl -minScore1 3000`.
#[test]
fn nets_the_input_chains_when_no_net_is_given() {
    let directory = temporary_directory();
    let mut config = config("nonet", &directory);
    config.net = None;
    config.reference = fixture("nonet_reference.fa");
    config.query = fixture("nonet_query.fa");
    config.lr_fold_threshold = 2.5;
    config.lr_fold_threshold_pairs = 10.0;
    config.min_broken_chain_score = 10_000.0;
    config.max_suspect_score = 100_000.0;
    config.dump_net = Some(directory.join("generated.net"));
    config.new_chain_id_dict = Some(directory.join("ids.tsv"));
    let output = config.output.clone();
    let bed = config.removed_bed.clone();

    chainc::run(config).unwrap();

    let net_text = fs::read_to_string(directory.join("generated.net")).unwrap();
    assert_eq!(
        net_text.lines().next().unwrap(),
        "net chr1 400",
        "generated NET:\n{net_text}"
    );
    assert_eq!(
        net_text
            .lines()
            .filter(|line| line.trim_start().starts_with("fill"))
            .count(),
        3,
        "generated NET:\n{net_text}"
    );

    // The suspect is the breaking chain's middle block between the two gaps
    // that hold the broken chain's pieces.
    let bed_text = fs::read_to_string(bed).unwrap();
    assert!(
        bed_text.starts_with("chr1\t200\t220\tbreakingChainID_1_"),
        "{bed_text}"
    );
    assert!(bed_text.ends_with("\t0,0,153\n"), "{bed_text}");
    assert_eq!(
        fs::read_to_string(directory.join("ids.tsv")).unwrap(),
        "3\t1\n"
    );

    // The breaking chain keeps its flanks joined by one gap, and the removed
    // suspect becomes a standalone 20 bp chain.
    let output_text = fs::read_to_string(output).unwrap();
    assert!(
        output_text.contains(" 1\n100\t220\t220\n80\n\n"),
        "{output_text}"
    );
    assert!(output_text.contains(" 3\n20\n\n"), "{output_text}");

    fs::remove_dir_all(directory).unwrap();
}

/// The generated NET must lead to exactly the same result as feeding that same
/// NET back in through `--net`.
#[test]
fn generated_net_matches_the_same_net_supplied_as_a_file() {
    let root = temporary_directory();
    let generated_directory = root.join("generated");
    let file_directory = root.join("file");
    fs::create_dir_all(&generated_directory).unwrap();
    fs::create_dir_all(&file_directory).unwrap();
    let net = generated_directory.join("generated.net");

    let mut generated = config("nonet", &generated_directory);
    generated.net = None;
    generated.reference = fixture("nonet_reference.fa");
    generated.query = fixture("nonet_query.fa");
    generated.lr_fold_threshold = 2.5;
    generated.dump_net = Some(net.clone());
    chainc::run(generated).unwrap();

    let mut from_file = config("nonet", &file_directory);
    from_file.net = Some(net);
    from_file.reference = fixture("nonet_reference.fa");
    from_file.query = fixture("nonet_query.fa");
    from_file.lr_fold_threshold = 2.5;
    chainc::run(from_file).unwrap();

    for file in ["clean.chain", "removed.bed"] {
        assert_eq!(
            fs::read(generated_directory.join(file)).unwrap(),
            fs::read(file_directory.join(file)).unwrap(),
            "{file} differs between the generated and file NET"
        );
    }

    fs::remove_dir_all(root).unwrap();
}
