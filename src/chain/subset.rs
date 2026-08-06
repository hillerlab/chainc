// Copyright (c) 2026 The Hiller Lab at the Senckenberg Gesellschaft für Naturforschung
// Distributed under the terms of the Apache License, Version 2.0.

use std::ops::Range;

use chaintools::AbsoluteBlock;

use super::EditableChain;

#[derive(Debug, Clone)]
pub struct SubchainView {
    pub chain_id: u64,
    pub revision: u32,
    pub requested: Range<u32>,
    pub blocks: Vec<AbsoluteBlock>,
}

impl EditableChain {
    /// Kent `chainSubsetOnT` semantics on the reference axis.
    pub fn subset_on_reference(&self, start: u32, end: u32) -> Option<SubchainView> {
        if start >= end {
            return None;
        }
        let first = self
            .blocks
            .partition_point(|block| block.reference_end <= start);
        let mut blocks = Vec::new();
        for block in &self.blocks[first..] {
            if block.reference_start >= end {
                break;
            }
            let mut clipped = *block;
            if clipped.reference_start < start {
                let amount = start - clipped.reference_start;
                clipped.reference_start = start;
                clipped.query_start += amount;
            }
            if clipped.reference_end > end {
                let amount = clipped.reference_end - end;
                clipped.reference_end = end;
                clipped.query_end -= amount;
            }
            if clipped.reference_start < clipped.reference_end {
                blocks.push(clipped);
            }
        }
        (!blocks.is_empty()).then_some(SubchainView {
            chain_id: self.id,
            revision: self.revision,
            requested: start..end,
            blocks,
        })
    }
}

#[cfg(test)]
mod tests {
    use chaintools::{AbsoluteBlock, Strand};

    use super::*;

    fn chain() -> EditableChain {
        EditableChain {
            id: 7,
            reference_name: b"chr1".to_vec(),
            reference_size: 100,
            reference_strand: Strand::Plus,
            reference_start: 10,
            reference_end: 50,
            query_name: b"q".to_vec(),
            query_size: 100,
            query_strand: Strand::Plus,
            query_start: 20,
            query_end: 60,
            header_score: 1,
            final_score: None,
            blocks: vec![
                AbsoluteBlock {
                    reference_start: 10,
                    reference_end: 20,
                    query_start: 20,
                    query_end: 30,
                },
                AbsoluteBlock {
                    reference_start: 30,
                    reference_end: 40,
                    query_start: 40,
                    query_end: 50,
                },
                AbsoluteBlock {
                    reference_start: 45,
                    reference_end: 50,
                    query_start: 55,
                    query_end: 60,
                },
            ],
            revision: 0,
            dirty: false,
            input_ordinal: 0,
        }
    }

    #[test]
    fn clips_first_and_last_blocks() {
        let view = chain().subset_on_reference(15, 47).unwrap();
        assert_eq!(
            view.blocks,
            vec![
                AbsoluteBlock {
                    reference_start: 15,
                    reference_end: 20,
                    query_start: 25,
                    query_end: 30,
                },
                AbsoluteBlock {
                    reference_start: 30,
                    reference_end: 40,
                    query_start: 40,
                    query_end: 50,
                },
                AbsoluteBlock {
                    reference_start: 45,
                    reference_end: 47,
                    query_start: 55,
                    query_end: 57,
                },
            ]
        );
    }

    #[test]
    fn no_overlap_is_none() {
        assert!(chain().subset_on_reference(20, 30).is_none());
    }
}
