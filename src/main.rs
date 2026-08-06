// Copyright (c) 2026 The Hiller Lab at the Senckenberg Gesellschaft für Naturforschung
// Distributed under the terms of the Apache License, Version 2.0.

use clap::Parser;

use chainc::config::CleanerConfig;
use chainc::error::Error;

fn main() {
    let config = CleanerConfig::parse();

    let result = match config.threads {
        Some(0) => Err(Error::Consistency(
            "--threads must be greater than zero".into(),
        )),
        Some(threads) => rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .map_err(|error| Error::Consistency(format!("cannot create thread pool: {error}")))
            .and_then(|pool| pool.install(|| chainc::run(config))),
        None => chainc::run(config),
    };

    if let Err(error) = result {
        eprintln!("chainc: {error}");
        std::process::exit(1);
    }
}
