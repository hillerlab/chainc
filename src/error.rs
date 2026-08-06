// Copyright (c) 2026 The Hiller Lab at the Senckenberg Gesellschaft für Naturforschung
// Distributed under the terms of the Apache License, Version 2.0.

use std::path::PathBuf;

#[derive(Debug, thiserror::Error)]
pub enum Error {
    #[error("I/O error for {path}: {source}")]
    Io {
        path: PathBuf,
        #[source]
        source: std::io::Error,
    },

    #[error("chain input error: {0}")]
    Chain(String),

    #[error("NET input error: {0}")]
    Net(String),

    #[error("sequence input error: {0}")]
    Sequence(String),

    #[error("score configuration error: {0}")]
    Score(String),

    #[error("inconsistent input: {0}")]
    Consistency(String),

    #[error("coordinate arithmetic overflow: {0}")]
    Overflow(String),
}

pub type Result<T> = std::result::Result<T, Error>;

pub(crate) fn io(path: impl Into<PathBuf>, source: std::io::Error) -> Error {
    Error::Io {
        path: path.into(),
        source,
    }
}
