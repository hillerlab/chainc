// Copyright (c) 2026 The Hiller Lab at the Senckenberg Gesellschaft für Naturforschung
// Distributed under the terms of the Apache License, Version 2.0.

use netools::{AlignmentSpan, NetRange};

#[derive(Debug, Clone)]
pub struct SpanIndex {
    spans: Vec<AlignmentSpan>,
    prefix_max_end: Vec<u32>,
}

impl SpanIndex {
    pub fn new(mut spans: Vec<AlignmentSpan>) -> Self {
        spans.sort_by(|left, right| {
            left.reference_range
                .start
                .cmp(&right.reference_range.start)
                .then(left.reference_range.end.cmp(&right.reference_range.end))
                .then(left.chain_id.cmp(&right.chain_id))
                .then(left.fill.cmp(&right.fill))
        });
        let mut maximum = 0;
        let prefix_max_end = spans
            .iter()
            .map(|span| {
                maximum = maximum.max(span.reference_range.end);
                maximum
            })
            .collect();
        Self {
            spans,
            prefix_max_end,
        }
    }

    pub fn any_overlap_where(
        &self,
        range: NetRange,
        predicate: impl Fn(&AlignmentSpan) -> bool,
    ) -> bool {
        if range.is_empty() {
            return false;
        }
        let first = self
            .prefix_max_end
            .partition_point(|maximum_end| *maximum_end <= range.start);
        self.spans[first..]
            .iter()
            .take_while(|span| span.reference_range.start < range.end)
            .any(|span| span.reference_range.overlaps(range) && predicate(span))
    }
}
