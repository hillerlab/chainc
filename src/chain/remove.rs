// Copyright (c) 2026 The Hiller Lab at the Senckenberg Gesellschaft für Naturforschung
// Distributed under the terms of the Apache License, Version 2.0.

use std::ops::Range;

use crate::error::{Error, Result};

use super::EditableChain;

impl EditableChain {
    /// Remove complete internal blocks following Kent `chainRemoveBlocks`.
    pub fn remove_internal_blocks(&mut self, suspect: Range<u32>) -> Result<Range<usize>> {
        let downstream_of_start = self
            .blocks
            .partition_point(|block| block.reference_start < suspect.start);
        if downstream_of_start == 0 {
            return Err(Error::Consistency(format!(
                "removing {}-{} would remove the first block of chain {}",
                suspect.start, suspect.end, self.id
            )));
        }
        let upstream = downstream_of_start - 1;
        let downstream = self.blocks[upstream + 1..]
            .iter()
            .position(|block| block.reference_start >= suspect.end)
            .map(|offset| upstream + 1 + offset)
            .ok_or_else(|| {
                Error::Consistency(format!(
                    "removing {}-{} would remove the final block of chain {}",
                    suspect.start, suspect.end, self.id
                ))
            })?;
        if downstream == upstream + 1 {
            return Err(Error::Consistency(format!(
                "suspect {}-{} contains no complete internal block in chain {}",
                suspect.start, suspect.end, self.id
            )));
        }
        let removed = upstream + 1..downstream;
        self.blocks.drain(removed.clone());
        self.revision = self
            .revision
            .checked_add(1)
            .ok_or_else(|| Error::Overflow(format!("chain {} revision", self.id)))?;
        self.dirty = true;
        Ok(removed)
    }
}

#[cfg(test)]
mod tests {
    use chaintools::{AbsoluteBlock, Strand};

    use super::*;

    fn chain() -> EditableChain {
        EditableChain {
            id: 1,
            reference_name: b"chr1".to_vec(),
            reference_size: 100,
            reference_strand: Strand::Plus,
            reference_start: 0,
            reference_end: 50,
            query_name: b"q".to_vec(),
            query_size: 100,
            query_strand: Strand::Plus,
            query_start: 0,
            query_end: 50,
            header_score: 1,
            final_score: None,
            blocks: (0..5)
                .map(|i| AbsoluteBlock {
                    reference_start: i * 10,
                    reference_end: i * 10 + 5,
                    query_start: i * 10,
                    query_end: i * 10 + 5,
                })
                .collect(),
            revision: 0,
            dirty: false,
            input_ordinal: 0,
        }
    }

    #[test]
    fn removes_only_internal_blocks() {
        let mut chain = chain();
        chain.remove_internal_blocks(10..30).unwrap();
        let starts: Vec<_> = chain
            .blocks
            .iter()
            .map(|block| block.reference_start)
            .collect();
        assert_eq!(starts, vec![0, 30, 40]);
        assert!(chain.dirty);
        assert_eq!(chain.revision, 1);
    }
}
