// Copyright (c) 2026 The Hiller Lab at the Senckenberg Gesellschaft für Naturforschung
// Distributed under the terms of the Apache License, Version 2.0.

use std::cmp::Ordering;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use crate::chain::{EditableChain, OriginalChain};
use crate::clean::RemovalEvent;
use crate::config::CleanerConfig;
use crate::error::{Error, Result, io};

enum OutputRecord<'a> {
    Original(&'a OriginalChain),
    Editable(&'a EditableChain),
}

impl OutputRecord<'_> {
    fn score(&self) -> i64 {
        match self {
            Self::Original(chain) => chain.score(),
            Self::Editable(chain) => chain.effective_score(),
        }
    }

    fn ordinal(&self) -> u64 {
        match self {
            Self::Original(chain) => chain.input_ordinal,
            Self::Editable(chain) => chain.input_ordinal,
        }
    }

    fn id(&self) -> u64 {
        match self {
            Self::Original(chain) => chain.id(),
            Self::Editable(chain) => chain.id,
        }
    }
}

#[allow(clippy::too_many_arguments)]
pub fn write_outputs(
    config: &CleanerConfig,
    chains: &HashMap<u64, EditableChain>,
    originals: &[OriginalChain],
    events: &[RemovalEvent],
    removed_bed: &[String],
    id_dictionary: &[String],
    suspect_data: &[String],
    metadata: &[Vec<u8>],
) -> Result<()> {
    let mut records: Vec<_> = originals.iter().map(OutputRecord::Original).collect();
    records.extend(chains.values().map(OutputRecord::Editable));
    records.extend(
        events
            .iter()
            .map(|event| OutputRecord::Editable(&event.removed_chain)),
    );
    records.sort_unstable_by(output_order);

    let mut writer = buffered_file(&config.output)?;
    chaintools::write_metadata_lines(&mut writer, metadata)
        .map_err(|error| Error::Chain(format!("{error:?}")))?;
    for record in records {
        let owned = match record {
            OutputRecord::Original(chain) => chain.to_owned_chain(),
            OutputRecord::Editable(chain) => chain.to_owned_chain()?,
        };
        chaintools::write_chain_dense(&mut writer, &owned)
            .map_err(|error| Error::Chain(format!("{error:?}")))?;
    }
    writer.flush().map_err(|error| io(&config.output, error))?;

    write_lines(&config.removed_bed, removed_bed)?;
    if let Some(path) = &config.new_chain_id_dict {
        write_lines(path, id_dictionary)?;
    }
    if let Some(path) = &config.suspect_data_file {
        write_lines(path, suspect_data)?;
    }
    Ok(())
}

fn buffered_file(path: &Path) -> Result<BufWriter<File>> {
    let file = File::create(path).map_err(|error| io(path, error))?;
    Ok(BufWriter::with_capacity(1 << 20, file))
}

fn write_lines(path: &Path, lines: &[String]) -> Result<()> {
    let mut writer = buffered_file(path)?;
    for line in lines {
        writer
            .write_all(line.as_bytes())
            .and_then(|_| writer.write_all(b"\n"))
            .map_err(|error| io(path, error))?;
    }
    writer.flush().map_err(|error| io(path, error))
}

/// Descending score, then a total tie-break the original program lacks: its
/// `chainSort` compares only the score through C `qsort`.
fn output_order(left: &OutputRecord<'_>, right: &OutputRecord<'_>) -> Ordering {
    right
        .score()
        .cmp(&left.score())
        .then_with(|| left.ordinal().cmp(&right.ordinal()))
        .then_with(|| left.id().cmp(&right.id()))
}
