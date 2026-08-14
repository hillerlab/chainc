// Copyright (c) 2026 The Hiller Lab at the Senckenberg Gesellschaft für Naturforschung
// Distributed under the terms of the Apache License, Version 2.0.

use std::collections::{HashMap, HashSet};
use std::path::Path;

use chaintools::seq::revcomp::reverse_complement_in_place;
use chaintools::seq::sequence::SequenceResolver;
use chaintools::{CompactMatrix, GapCalc, ScoreMatrix, score_absolute_block};

use crate::chain::{EditableChain, SubchainView};
use crate::error::{Error, Result};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SubchainMetrics {
    pub global_score: i64,
    pub local_score: i64,
    pub aligned_bases: u64,
}

#[derive(Debug, Clone)]
pub struct SequenceSet {
    reference: SequenceResolver,
    query: SequenceResolver,
    reverse_queries: HashMap<Vec<u8>, Vec<u8>>,
}

impl SequenceSet {
    /// Load only the sequences the chains of interest need.
    ///
    /// `reverse_query_names` holds the query sequences used by at least one
    /// minus-strand chain; only those are reverse complemented, mirroring the
    /// original program's on-demand minus-strand cache.
    pub fn load(
        reference_path: &Path,
        query_path: &Path,
        reference_names: &HashSet<Vec<u8>>,
        query_names: &HashSet<Vec<u8>>,
        reverse_query_names: &HashSet<Vec<u8>>,
    ) -> Result<Self> {
        let reference = SequenceResolver::preload(reference_path, Some(reference_names))
            .map_err(|error| Error::Sequence(format!("{error:?}")))?;
        let query = SequenceResolver::preload(query_path, Some(query_names))
            .map_err(|error| Error::Sequence(format!("{error:?}")))?;
        let mut reverse_queries = HashMap::with_capacity(reverse_query_names.len());
        for name in reverse_query_names {
            let mut sequence = query
                .sequence(name)
                .map_err(|error| Error::Sequence(format!("{error:?}")))?
                .to_vec();
            reverse_complement_in_place(&mut sequence);
            reverse_queries.insert(name.clone(), sequence);
        }
        Ok(Self {
            reference,
            query,
            reverse_queries,
        })
    }

    fn sequences_for<'a>(&'a self, chain: &EditableChain) -> Result<(&'a [u8], &'a [u8])> {
        let reference = self
            .reference
            .sequence(&chain.reference_name)
            .map_err(|error| Error::Sequence(format!("{error:?}")))?;
        let query = match chain.query_strand {
            chaintools::Strand::Plus => self
                .query
                .sequence(&chain.query_name)
                .map_err(|error| Error::Sequence(format!("{error:?}")))?,
            chaintools::Strand::Minus => self
                .reverse_queries
                .get(&chain.query_name)
                .map(Vec::as_slice)
                .ok_or_else(|| {
                    Error::Sequence(format!(
                        "missing reverse-complemented query {}",
                        String::from_utf8_lossy(&chain.query_name)
                    ))
                })?,
        };
        if reference.len() != chain.reference_size as usize {
            return Err(Error::Consistency(format!(
                "reference sequence {} has length {}, chain {} declares {}",
                String::from_utf8_lossy(&chain.reference_name),
                reference.len(),
                chain.id,
                chain.reference_size
            )));
        }
        if query.len() != chain.query_size as usize {
            return Err(Error::Consistency(format!(
                "query sequence {} has length {}, chain {} declares {}",
                String::from_utf8_lossy(&chain.query_name),
                query.len(),
                chain.id,
                chain.query_size
            )));
        }
        Ok((reference, query))
    }
}

#[derive(Debug, Clone)]
pub struct ParityScorer {
    matrix: CompactMatrix,
    gap: GapCalc,
    sequences: SequenceSet,
}

impl ParityScorer {
    pub fn new(matrix: ScoreMatrix, gap: GapCalc, sequences: SequenceSet) -> Self {
        Self {
            matrix: matrix.compact(),
            gap,
            sequences,
        }
    }

    pub fn score(&self, chain: &EditableChain, view: &SubchainView) -> Result<SubchainMetrics> {
        let (reference, query) = self.sequences.sequences_for(chain)?;
        let mut global_score = 0i64;
        let mut running_score = 0i64;
        let mut local_score = 0i64;
        let mut aligned_bases = 0u64;

        for (index, block) in view.blocks.iter().copied().enumerate() {
            let block_score = score_absolute_block(block, query, reference, &self.matrix)
                .map_err(|error| Error::Score(format!("{error:?}")))?;
            global_score = global_score
                .checked_add(block_score)
                .ok_or_else(|| Error::Overflow(format!("global score for chain {}", chain.id)))?;
            running_score = running_score
                .checked_add(block_score)
                .ok_or_else(|| Error::Overflow(format!("local score for chain {}", chain.id)))?;
            local_score = local_score.max(running_score);
            aligned_bases = aligned_bases
                .checked_add(u64::from(block.reference_len()))
                .ok_or_else(|| Error::Overflow(format!("aligned bases for chain {}", chain.id)))?;

            if let Some(next) = view.blocks.get(index + 1) {
                let reference_gap = next
                    .reference_start
                    .checked_sub(block.reference_end)
                    .ok_or_else(|| {
                        Error::Consistency(format!(
                            "overlapping reference blocks in chain {}",
                            chain.id
                        ))
                    })?;
                let query_gap = next
                    .query_start
                    .checked_sub(block.query_end)
                    .ok_or_else(|| {
                        Error::Consistency(format!(
                            "overlapping query blocks in chain {}",
                            chain.id
                        ))
                    })?;
                let reference_gap = i32::try_from(reference_gap)
                    .map_err(|_| Error::Overflow(format!("reference gap in chain {}", chain.id)))?;
                let query_gap = i32::try_from(query_gap)
                    .map_err(|_| Error::Overflow(format!("query gap in chain {}", chain.id)))?;
                let cost = i64::from(self.gap.cost(query_gap, reference_gap));
                global_score = global_score.checked_sub(cost).ok_or_else(|| {
                    Error::Overflow(format!("global gap score for chain {}", chain.id))
                })?;
                running_score = running_score.checked_sub(cost).ok_or_else(|| {
                    Error::Overflow(format!("local gap score for chain {}", chain.id))
                })?;
                if running_score < 0 {
                    running_score = 0;
                }
            }
        }
        Ok(SubchainMetrics {
            global_score,
            local_score,
            aligned_bases,
        })
    }
}

pub fn load_score_matrix(path: Option<&Path>) -> Result<ScoreMatrix> {
    path.map_or_else(
        || Ok(ScoreMatrix::default_dna()),
        |path| {
            ScoreMatrix::from_blastz_file(path).map_err(|error| Error::Score(format!("{error:?}")))
        },
    )
}

pub fn load_gap_model(value: &str) -> Result<GapCalc> {
    GapCalc::from_linear_gap(value).map_err(|error| Error::Score(format!("{error:?}")))
}

#[cfg(test)]
mod tests {
    use chaintools::{AbsoluteBlock, GapCalc, ScoreMatrix};

    #[test]
    fn default_components_match_known_kent_values() {
        let matrix = ScoreMatrix::default_dna();
        assert_eq!(matrix.pair(b'A', b'A'), 91);
        assert_eq!(matrix.pair(b'C', b'C'), 100);
        assert_eq!(GapCalc::default_costs().cost(2, 0), 360);
        let block = AbsoluteBlock {
            reference_start: 0,
            reference_end: 4,
            query_start: 0,
            query_end: 4,
        };
        assert_eq!(
            chaintools::score_absolute_block(block, b"ACGT", b"ACGT", &matrix.compact()).unwrap(),
            382
        );
    }
}
