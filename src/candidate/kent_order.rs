// Copyright (c) 2026 The Hiller Lab at the Senckenberg Gesellschaft für Naturforschung
// Distributed under the terms of the Apache License, Version 2.0.

use std::collections::HashSet;

const INITIAL_POWER: u32 = 12;

fn kent_hash(value: u64) -> u32 {
    value.to_string().bytes().fold(0u32, |result, byte| {
        result
            .wrapping_add(result.wrapping_shl(3))
            .wrapping_add(u32::from(byte))
    })
}

#[derive(Debug)]
struct KentHash {
    buckets: Vec<Vec<u64>>,
    members: HashSet<u64>,
}

impl KentHash {
    fn new() -> Self {
        Self {
            buckets: vec![Vec::new(); 1usize << INITIAL_POWER],
            members: HashSet::new(),
        }
    }

    fn insert_unique(&mut self, value: u64) {
        if !self.members.insert(value) {
            return;
        }
        let bucket = kent_hash(value) as usize & (self.buckets.len() - 1);
        self.buckets[bucket].insert(0, value);
        if self.members.len() > self.buckets.len() {
            self.resize();
        }
    }

    fn resize(&mut self) {
        let mut resized = vec![Vec::new(); self.buckets.len() * 2];
        for value in self.buckets.iter().flatten().copied() {
            let bucket = kent_hash(value) as usize & (resized.len() - 1);
            resized[bucket].insert(0, value);
        }
        for bucket in &mut resized {
            bucket.reverse();
        }
        self.buckets = resized;
    }

    fn traverse(&self) -> Vec<u64> {
        self.buckets.iter().flatten().copied().collect()
    }

    fn element_list(&self) -> Vec<u64> {
        let mut values = self.traverse();
        values.reverse();
        values
    }
}

/// Order used by `hashTraverseEls`.
pub fn traversal_order(insertion_order: impl IntoIterator<Item = u64>) -> Vec<u64> {
    let mut hash = KentHash::new();
    for value in insertion_order {
        hash.insert_unique(value);
    }
    hash.traverse()
}

/// Order used by `hashElListHash`.
pub fn element_list_order(insertion_order: impl IntoIterator<Item = u64>) -> Vec<u64> {
    let mut hash = KentHash::new();
    for value in insertion_order {
        hash.insert_unique(value);
    }
    hash.element_list()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn element_list_reverses_bucket_traversal() {
        let inserted = [1, 2, 3, 4, 5];
        let mut traversed = traversal_order(inserted);
        traversed.reverse();
        assert_eq!(element_list_order(inserted), traversed);
    }
}
