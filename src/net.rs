// Copyright (c) 2026 The Hiller Lab at the Senckenberg Gesellschaft für Naturforschung
// Distributed under the terms of the Apache License, Version 2.0.

//! NET acquisition: either a user-supplied file or an in-memory NET generated
//! from the input chains.
//!
//! The original program shells out to `chainNet -minScore=0` followed by
//! `NetFilterNonNested.perl -minScore1 3000` and reads the resulting temporary
//! file back. `netools` 0.0.2 ports both steps, so the same NET is built in
//! memory and consumed directly: no temporary file, no shell command, and the
//! raw (unfiltered) NET is never materialized.

use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};

use netools::chainnet::{ChainNetBuilder, ChainNetOptions, SequenceSizes, SizeSource};
use netools::{Net, NetRef, OwnedNet, Reader, Severity, ValidationMode, ValidationReport};

use crate::error::{Error, Result};

/// Where the reference-side NET sections come from.
pub enum NetSource {
    /// Sections parsed from a `--net` file.
    Parsed(Reader<Net>),
    /// Sections built from the input chains by the ported `chainNet` +
    /// `NetFilterNonNested` pipeline.
    Generated(Vec<OwnedNet>),
}

impl NetSource {
    /// Parse a NET file with mmap and parallel section parsing.
    pub fn from_file(path: &Path) -> Result<Self> {
        let reader = Reader::<Net>::from_path_parallel(path)
            .map_err(|error| Error::Net(format!("{path}: {error}", path = path.display())))?;
        Ok(Self::Parsed(reader))
    }

    /// Build the cleaner-compatible NET from already parsed chains.
    ///
    /// This is `chainNet -minScore=0` restricted to the reference side, fused
    /// with `NetFilterNonNested.perl -minScore1 3000`.
    pub fn generate(
        chains: &chaintools::Reader<chaintools::Chain>,
        reference_sizes: SizeSource,
        query_sizes: SizeSource,
    ) -> Result<Self> {
        // Threads are left unset so section construction runs on the pool this
        // process already installed instead of nesting a second one.
        let options = ChainNetOptions::chaincleaner_compatible();
        let generated = ChainNetBuilder::new(options)
            .reference_sizes(reference_sizes)
            .query_sizes(query_sizes)
            .build(chains)
            .map_err(|error| Error::Net(format!("cannot net the input chains: {error}")))?;
        Ok(Self::Generated(generated.into_reference()))
    }

    /// Reference-side sections in NET order.
    pub fn sections(&self) -> Box<dyn Iterator<Item = NetRef<'_>> + '_> {
        match self {
            Self::Parsed(reader) => Box::new(reader.nets()),
            Self::Generated(nets) => Box::new(nets.iter().map(OwnedNet::as_ref)),
        }
    }

    /// Reject structurally invalid input.
    ///
    /// Only file input is validated: generated sections come from a builder
    /// that maintains the invariants by construction, and a second full pass
    /// over millions of records is pure overhead.
    pub fn validate(&self) -> Result<()> {
        let Self::Parsed(reader) = self else {
            return Ok(());
        };
        let report = reader.validate(ValidationMode::Strict);
        first_error(&report)
    }

    /// Write the sections back out, mainly to inspect a generated NET.
    pub fn write(&self, path: &Path) -> Result<()> {
        let mut writer = netools::Writer::from_path(path)
            .map_err(|error| Error::Net(format!("{}: {error}", path.display())))?;
        for net in self.sections() {
            writer
                .write_net(net)
                .map_err(|error| Error::Net(format!("{}: {error}", path.display())))?;
        }
        writer
            .finish()
            .map_err(|error| Error::Net(format!("{}: {error}", path.display())))
    }
}

fn first_error(report: &ValidationReport) -> Result<()> {
    let Some(issue) = report
        .issues()
        .iter()
        .find(|issue| issue.severity == Severity::Error)
    else {
        return Ok(());
    };
    Err(Error::Net(format!(
        "{} in NET section {}{}: {}",
        issue.code.as_str(),
        issue.net.get(),
        issue
            .node
            .map(|node| format!(", node {}", node.get()))
            .unwrap_or_default(),
        issue.message
    )))
}

/// Resolve the ordered sequence sizes used for NET construction.
///
/// Section order follows this source, exactly as `chainNet` follows the order
/// of its `tSizes` file. An explicit two-column file wins; otherwise the sizes
/// come from the sequence file itself, which for `.2bit` means the index only.
pub fn size_source(explicit: Option<&PathBuf>, sequence: &Path) -> Result<SizeSource> {
    if let Some(path) = explicit {
        return Ok(SizeSource::ChromSizes(path.clone()));
    }
    if is_two_bit(sequence) {
        return Ok(SizeSource::TwoBit(sequence.to_path_buf()));
    }
    Ok(SizeSource::Provided(fasta_sizes(sequence)?))
}

fn is_two_bit(path: &Path) -> bool {
    path.extension()
        .is_some_and(|extension| extension.eq_ignore_ascii_case("2bit"))
}

/// Collect FASTA record names and lengths in file order without keeping any
/// sequence, so a NET can be generated before sequences are loaded.
fn fasta_sizes(path: &Path) -> Result<SequenceSizes> {
    let file = std::fs::File::open(path).map_err(|error| crate::error::io(path, error))?;
    let reader: Box<dyn BufRead> = if path
        .extension()
        .is_some_and(|extension| extension.eq_ignore_ascii_case("gz"))
    {
        Box::new(BufReader::with_capacity(
            1 << 20,
            flate2::read::MultiGzDecoder::new(file),
        ))
    } else {
        Box::new(BufReader::with_capacity(1 << 20, file))
    };

    let mut entries: Vec<(Vec<u8>, u32)> = Vec::new();
    let mut length = 0u64;
    for line in reader.split(b'\n') {
        let line = line.map_err(|error| crate::error::io(path, error))?;
        let line = line.strip_suffix(b"\r").unwrap_or(&line);
        if let Some(header) = line.strip_prefix(b">") {
            if let Some(last) = entries.last_mut() {
                last.1 = finite_length(length, &last.0)?;
            }
            length = 0;
            let name = header
                .split(|byte| byte.is_ascii_whitespace())
                .next()
                .unwrap_or_default()
                .to_vec();
            if name.is_empty() {
                return Err(Error::Sequence(format!(
                    "{}: FASTA record without a name",
                    path.display()
                )));
            }
            entries.push((name, 0));
        } else {
            length += line
                .iter()
                .filter(|byte| !byte.is_ascii_whitespace())
                .count() as u64;
        }
    }
    if let Some(last) = entries.last_mut() {
        last.1 = finite_length(length, &last.0)?;
    }
    if entries.is_empty() {
        return Err(Error::Sequence(format!(
            "{}: no FASTA records",
            path.display()
        )));
    }
    SequenceSizes::new(entries).map_err(|error| Error::Sequence(format!("{error}")))
}

fn finite_length(length: u64, name: &[u8]) -> Result<u32> {
    u32::try_from(length).map_err(|_| {
        Error::Overflow(format!(
            "sequence {} is longer than 2^32-1 bases",
            String::from_utf8_lossy(name)
        ))
    })
}
