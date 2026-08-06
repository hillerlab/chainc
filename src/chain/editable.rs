// Copyright (c) 2026 The Hiller Lab at the Senckenberg Gesellschaft für Naturforschung
// Distributed under the terms of the Apache License, Version 2.0.

use chaintools::{AbsoluteBlock, Chain, OwnedChain, Strand, absolute_to_dense_blocks};

use crate::error::{Error, Result};

/// Validate dense coordinates without materializing absolute blocks.
pub fn validate_parsed_chain(chain: &Chain) -> Result<()> {
    if chain.reference_strand != Strand::Plus {
        return Err(Error::Consistency(format!(
            "chain {} has unsupported reference strand '-'",
            chain.id
        )));
    }
    if chain.reference_start > chain.reference_end
        || chain.reference_end > chain.reference_size
        || chain.query_start > chain.query_end
        || chain.query_end > chain.query_size
    {
        return Err(Error::Consistency(format!(
            "chain {} has invalid header coordinates",
            chain.id
        )));
    }
    let dense = chain.blocks.as_slice();
    if dense.is_empty() {
        return Err(Error::Consistency(format!(
            "chain {} contains no blocks",
            chain.id
        )));
    }
    let mut reference_position = chain.reference_start;
    let mut query_position = chain.query_start;
    for (index, block) in dense.iter().enumerate() {
        if block.size == 0 {
            return Err(Error::Consistency(format!(
                "chain {} block {} has zero size",
                chain.id, index
            )));
        }
        reference_position = reference_position.checked_add(block.size).ok_or_else(|| {
            Error::Overflow(format!("chain {} reference block {}", chain.id, index))
        })?;
        query_position = query_position
            .checked_add(block.size)
            .ok_or_else(|| Error::Overflow(format!("chain {} query block {}", chain.id, index)))?;
        if index + 1 < dense.len() {
            reference_position = reference_position
                .checked_add(block.gap_reference)
                .ok_or_else(|| {
                    Error::Overflow(format!("chain {} reference gap {}", chain.id, index))
                })?;
            query_position = query_position.checked_add(block.gap_query).ok_or_else(|| {
                Error::Overflow(format!("chain {} query gap {}", chain.id, index))
            })?;
        } else if block.gap_reference != 0 || block.gap_query != 0 {
            return Err(Error::Consistency(format!(
                "chain {} final block has non-zero gap",
                chain.id
            )));
        }
    }
    if reference_position != chain.reference_end || query_position != chain.query_end {
        return Err(Error::Consistency(format!(
            "chain {} block coordinates end at reference/query {}/{} but header ends at {}/{}",
            chain.id, reference_position, query_position, chain.reference_end, chain.query_end
        )));
    }
    Ok(())
}

/// Cheap mmap-backed handle retained for chains that cannot be modified.
#[derive(Debug, Clone)]
pub struct OriginalChain {
    parsed: Chain,
    pub input_ordinal: u64,
}

impl OriginalChain {
    pub fn new(parsed: &Chain, input_ordinal: u64) -> Self {
        Self {
            parsed: parsed.clone(),
            input_ordinal,
        }
    }

    pub fn id(&self) -> u64 {
        self.parsed.id
    }

    pub fn score(&self) -> i64 {
        self.parsed.score
    }

    pub fn to_owned_chain(&self) -> OwnedChain {
        OwnedChain {
            score: self.parsed.score,
            reference_name: self.parsed.reference_name.as_bytes().to_vec(),
            reference_size: self.parsed.reference_size,
            reference_strand: self.parsed.reference_strand,
            reference_start: self.parsed.reference_start,
            reference_end: self.parsed.reference_end,
            query_name: self.parsed.query_name.as_bytes().to_vec(),
            query_size: self.parsed.query_size,
            query_strand: self.parsed.query_strand,
            query_start: self.parsed.query_start,
            query_end: self.parsed.query_end,
            id: self.parsed.id,
            blocks: self.parsed.blocks.as_slice().to_vec(),
        }
    }
}

#[derive(Debug, Clone)]
pub struct EditableChain {
    pub id: u64,
    pub reference_name: Vec<u8>,
    pub reference_size: u32,
    pub reference_strand: Strand,
    pub reference_start: u32,
    pub reference_end: u32,
    pub query_name: Vec<u8>,
    pub query_size: u32,
    pub query_strand: Strand,
    pub query_start: u32,
    pub query_end: u32,
    pub header_score: i64,
    pub final_score: Option<i64>,
    pub blocks: Vec<AbsoluteBlock>,
    pub revision: u32,
    pub dirty: bool,
    pub input_ordinal: u64,
}

impl EditableChain {
    /// Materialize a chain of interest as absolute blocks.
    ///
    /// [`validate_parsed_chain`] has already proved that every block size is
    /// non-zero, that no cumulative coordinate overflows, and that the blocks
    /// end exactly at the header bounds, so the expansion below only walks.
    pub fn from_parsed(chain: &Chain, input_ordinal: u64) -> Self {
        let mut reference_position = chain.reference_start;
        let mut query_position = chain.query_start;
        let dense = chain.blocks.as_slice();
        let mut blocks = Vec::with_capacity(dense.len());
        for block in dense {
            let reference_end = reference_position + block.size;
            let query_end = query_position + block.size;
            blocks.push(AbsoluteBlock {
                reference_start: reference_position,
                reference_end,
                query_start: query_position,
                query_end,
            });
            reference_position = reference_end + block.gap_reference;
            query_position = query_end + block.gap_query;
        }

        Self {
            id: chain.id,
            reference_name: chain.reference_name.as_bytes().to_vec(),
            reference_size: chain.reference_size,
            reference_strand: chain.reference_strand,
            reference_start: chain.reference_start,
            reference_end: chain.reference_end,
            query_name: chain.query_name.as_bytes().to_vec(),
            query_size: chain.query_size,
            query_strand: chain.query_strand,
            query_start: chain.query_start,
            query_end: chain.query_end,
            header_score: chain.score,
            final_score: None,
            blocks,
            revision: 0,
            dirty: false,
            input_ordinal,
        }
    }

    pub fn effective_score(&self) -> i64 {
        self.final_score.unwrap_or(self.header_score)
    }

    pub fn to_owned_chain(&self) -> Result<OwnedChain> {
        let blocks = absolute_to_dense_blocks(&self.blocks)
            .map_err(|error| Error::Chain(format!("{error:?}")))?;
        Ok(OwnedChain {
            score: self.effective_score(),
            reference_name: self.reference_name.clone(),
            reference_size: self.reference_size,
            reference_strand: self.reference_strand,
            reference_start: self.reference_start,
            reference_end: self.reference_end,
            query_name: self.query_name.clone(),
            query_size: self.query_size,
            query_strand: self.query_strand,
            query_start: self.query_start,
            query_end: self.query_end,
            id: self.id,
            blocks,
        })
    }

    pub(crate) fn owned_from_blocks(
        &self,
        id: u64,
        score: i64,
        blocks: Vec<AbsoluteBlock>,
        input_ordinal: u64,
    ) -> Result<EditableChain> {
        let first = blocks.first().ok_or_else(|| {
            Error::Consistency(format!("cannot create empty subset of chain {}", self.id))
        })?;
        let last = blocks.last().expect("checked non-empty");
        Ok(EditableChain {
            id,
            reference_name: self.reference_name.clone(),
            reference_size: self.reference_size,
            reference_strand: self.reference_strand,
            reference_start: first.reference_start,
            reference_end: last.reference_end,
            query_name: self.query_name.clone(),
            query_size: self.query_size,
            query_strand: self.query_strand,
            query_start: first.query_start,
            query_end: last.query_end,
            header_score: score,
            final_score: Some(score),
            blocks,
            revision: 0,
            dirty: false,
            input_ordinal,
        })
    }
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use chaintools::{Block, BlockSlice, ByteSlice, Chain};

    use super::*;

    fn bytes(value: &[u8]) -> ByteSlice {
        let storage = chaintools::io::storage::SharedBytes::from_owned(value.to_vec());
        ByteSlice::new(storage, 0..value.len())
    }

    #[test]
    fn expands_dense_blocks() {
        let blocks = Arc::new(vec![
            Block {
                size: 10,
                gap_reference: 5,
                gap_query: 7,
            },
            Block {
                size: 4,
                gap_reference: 0,
                gap_query: 0,
            },
        ]);
        let parsed = Chain {
            score: 42,
            reference_name: bytes(b"chr1"),
            reference_size: 100,
            reference_strand: Strand::Plus,
            reference_start: 3,
            reference_end: 22,
            query_name: bytes(b"chrQ"),
            query_size: 100,
            query_strand: Strand::Plus,
            query_start: 8,
            query_end: 29,
            id: 1,
            blocks: BlockSlice::new(blocks, 0..2),
        };
        let chain = EditableChain::from_parsed(&parsed, 0);
        assert_eq!(
            chain.blocks,
            vec![
                AbsoluteBlock {
                    reference_start: 3,
                    reference_end: 13,
                    query_start: 8,
                    query_end: 18,
                },
                AbsoluteBlock {
                    reference_start: 18,
                    reference_end: 22,
                    query_start: 25,
                    query_end: 29,
                }
            ]
        );
    }
}
