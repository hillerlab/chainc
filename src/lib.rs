// Copyright (c) 2026 The Hiller Lab at the Senckenberg Gesellschaft für Naturforschung
// Distributed under the terms of the Apache License, Version 2.0.

pub mod candidate;
pub mod chain;
pub mod clean;
pub mod config;
pub mod error;
pub mod net;
pub mod output;
pub mod score;

use std::collections::{HashMap, HashSet};
use std::fs;

use candidate::discover_candidates;
use chain::{EditableChain, OriginalChain, validate_parsed_chain};
use clean::Cleaner;
use config::CleanerConfig;
use error::{Error, Result, io};
use net::NetSource;
use score::{SequenceSet, load_gap_model, load_score_matrix};

#[derive(Debug)]
struct ChainHeader {
    id: u64,
    score: i64,
    reference_name: Vec<u8>,
    reference_size: u32,
}

pub fn run(config: CleanerConfig) -> Result<()> {
    let chain_reader = chaintools::Reader::<chaintools::Chain>::from_path_parallel(&config.chains)
        .map_err(|e| Error::Chain(format!("{e:?}")))?;
    let metadata: Vec<Vec<u8>> = chain_reader.metadata_lines().map(<[u8]>::to_vec).collect();

    let nets = match &config.net {
        Some(path) => NetSource::from_file(path)?,
        None => NetSource::generate(
            &chain_reader,
            net::size_source(config.reference_sizes.as_ref(), &config.reference)?,
            net::size_source(config.query_sizes.as_ref(), &config.query)?,
        )?,
    };
    nets.validate()?;
    if let Some(path) = &config.dump_net {
        nets.write(path)?;
    }
    let discovery = discover_candidates(nets.sections())?;

    if let Some(path) = &config.dump_candidate_order {
        let trace = discovery.trace();
        fs::write(path, trace).map_err(|e| io(path, e))?;
    }

    let mut chains = HashMap::with_capacity(discovery.interest_chain_ids.len());
    let mut originals = Vec::with_capacity(
        chain_reader
            .len()
            .saturating_sub(discovery.interest_chain_ids.len()),
    );
    let mut headers = HashMap::with_capacity(chain_reader.len());
    let mut maximum_chain_id = None;
    let mut reference_names = HashSet::new();
    let mut query_names = HashSet::new();
    let mut reverse_query_names = HashSet::new();

    for (ordinal, parsed) in chain_reader.chains().enumerate() {
        validate_parsed_chain(parsed)?;
        if headers.contains_key(&parsed.id) {
            return Err(Error::Consistency(format!(
                "duplicate chain ID {}",
                parsed.id
            )));
        }
        headers.insert(
            parsed.id,
            ChainHeader {
                id: parsed.id,
                score: parsed.score,
                reference_name: parsed.reference_name.as_bytes().to_vec(),
                reference_size: parsed.reference_size,
            },
        );
        maximum_chain_id = Some(maximum_chain_id.map_or(parsed.id, |old: u64| old.max(parsed.id)));
        if discovery.interest_chain_ids.contains(&parsed.id) {
            let editable = EditableChain::from_parsed(parsed, ordinal as u64);
            reference_names.insert(editable.reference_name.clone());
            if editable.query_strand == chaintools::Strand::Minus {
                reverse_query_names.insert(editable.query_name.clone());
            }
            query_names.insert(editable.query_name.clone());
            chains.insert(parsed.id, editable);
        } else {
            originals.push(OriginalChain::new(parsed, ordinal as u64));
        }
    }

    for chain_id in &discovery.net_chain_ids {
        if !headers.contains_key(chain_id) {
            return Err(Error::Consistency(format!(
                "NET references missing chain ID {chain_id}"
            )));
        }
    }

    validate_net_sizes(&nets, &headers)?;
    if config.strict {
        validate_id_ranking(&headers)?;
    }

    let maximum_chain_id =
        maximum_chain_id.ok_or_else(|| Error::Consistency("chain input is empty".into()))?;
    let matrix = load_score_matrix(config.score_scheme.as_deref())?;
    let gaps = load_gap_model(&config.linear_gap)?;
    let sequences = SequenceSet::load(
        &config.reference,
        &config.query,
        &reference_names,
        &query_names,
        &reverse_query_names,
    )?;
    Cleaner::new(
        config,
        discovery,
        chains,
        originals,
        chain_reader.len() as u64,
        maximum_chain_id,
        matrix,
        gaps,
        sequences,
        metadata,
    )
    .run()
}

fn validate_net_sizes(nets: &NetSource, chains: &HashMap<u64, ChainHeader>) -> Result<()> {
    for net in nets.sections() {
        for fill in net.fills() {
            if let Some(id) = fill.chain_id()
                && let Some(chain) = chains.get(&id)
            {
                if chain.reference_name.as_slice() != net.reference_name_bytes() {
                    return Err(Error::Consistency(format!(
                        "NET section {} contains chain {} from reference {}",
                        String::from_utf8_lossy(net.reference_name_bytes()),
                        id,
                        String::from_utf8_lossy(&chain.reference_name)
                    )));
                }
                if chain.reference_size != net.reference_size() {
                    return Err(Error::Consistency(format!(
                        "reference size mismatch for chain {id}: chain {}, NET {}",
                        chain.reference_size,
                        net.reference_size()
                    )));
                }
            }
        }
    }
    Ok(())
}

fn validate_id_ranking(chains: &HashMap<u64, ChainHeader>) -> Result<()> {
    let mut ordered: Vec<_> = chains.values().collect();
    ordered.sort_by_key(|chain| chain.id);
    for pair in ordered.windows(2) {
        if pair[0].score < pair[1].score {
            return Err(Error::Consistency(format!(
                "strict mode: chain ID ranking is inconsistent (ID {} score {} precedes ID {} score {})",
                pair[0].id, pair[0].score, pair[1].id, pair[1].score
            )));
        }
    }
    Ok(())
}
