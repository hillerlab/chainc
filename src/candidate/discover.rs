// Copyright (c) 2026 The Hiller Lab at the Senckenberg Gesellschaft für Naturforschung
// Distributed under the terms of the Apache License, Version 2.0.

use std::collections::{HashMap, HashSet};
use std::ops::Range;
use std::sync::Arc;

use netools::{NestedFillGapContext, NetRange, NetRef};

use crate::error::{Error, Result};

use super::interval_index::SpanIndex;
use super::kent_order::{element_list_order, traversal_order};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BreakCandidate {
    pub ordinal: u64,
    pub net_id: u32,
    pub depth: u16,
    pub reference_name: Arc<[u8]>,
    pub broken_chain_id: u64,
    pub breaking_chain_id: u64,
    pub left_fill: Range<u32>,
    pub right_fill: Range<u32>,
    pub left_gap: Range<u32>,
    pub right_gap: Range<u32>,
    pub suspect: Range<u32>,
}

#[derive(Debug, Clone)]
pub struct CandidateGroup {
    pub compatibility_ordinal: u64,
    pub breaking_chain_id: u64,
    pub candidates: Vec<BreakCandidate>,
}

#[derive(Debug)]
pub struct CandidateDiscovery {
    pub groups: Vec<CandidateGroup>,
    pub interest_chain_ids: HashSet<u64>,
    pub net_chain_ids: HashSet<u64>,
}

impl CandidateDiscovery {
    pub fn trace(&self) -> Vec<u8> {
        let mut output = String::new();
        for group in &self.groups {
            for candidate in &group.candidates {
                use std::fmt::Write;
                let _ = writeln!(
                    output,
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    group.compatibility_ordinal,
                    candidate.ordinal,
                    candidate.breaking_chain_id,
                    candidate.broken_chain_id,
                    String::from_utf8_lossy(candidate.reference_name.as_ref()),
                    candidate.depth,
                    candidate.left_fill.start,
                    candidate.left_fill.end,
                    candidate.right_fill.start,
                    candidate.right_fill.end
                );
            }
        }
        output.into_bytes()
    }
}

#[derive(Clone, Copy)]
struct OrderedContext {
    context: NestedFillGapContext,
    net_id: u32,
}

/// Discover break candidates from reference NET sections in NET order.
///
/// Sections are numbered by position rather than by their own `NetId` so that
/// generated single-section stores, which all report id 0, stay distinct.
pub fn discover_candidates<'a, I>(sections: I) -> Result<CandidateDiscovery>
where
    I: IntoIterator<Item = NetRef<'a>>,
{
    let mut contexts_by_chain: HashMap<u64, Vec<OrderedContext>> = HashMap::new();
    let mut nested_chain_insertion_order = Vec::new();
    let mut indexes = HashMap::new();
    let mut reference_names: HashMap<u32, Arc<[u8]>> = HashMap::new();
    let mut net_chain_ids = HashSet::new();

    for (ordinal, net) in sections.into_iter().enumerate() {
        let net_id = u32::try_from(ordinal)
            .map_err(|_| Error::Overflow("more than 2^32-1 NET sections".into()))?;
        reference_names.insert(net_id, Arc::from(net.reference_name_bytes()));
        net_chain_ids.extend(net.used_chain_ids());
        let mut spans = Vec::new();
        for fill in net.fills() {
            spans.extend(fill.uninterrupted_reference_spans());
        }
        indexes.insert(net_id, SpanIndex::new(spans));
        for context in net.nested_fill_gap_contexts() {
            nested_chain_insertion_order.push(context.fill_chain_id);
            contexts_by_chain
                .entry(context.fill_chain_id)
                .or_default()
                .push(OrderedContext { context, net_id });
        }
    }

    let broken_chain_order = traversal_order(nested_chain_insertion_order);
    let mut candidates_by_breaker: HashMap<u64, Vec<BreakCandidate>> = HashMap::new();
    let mut breaker_insertion_order = Vec::new();
    let mut interest_chain_ids = HashSet::new();
    let mut next_ordinal = 0u64;

    for broken_chain_id in broken_chain_order {
        let Some(contexts) = contexts_by_chain.get(&broken_chain_id) else {
            continue;
        };
        if contexts.len() < 2 {
            continue;
        }
        for pair in contexts.windows(2) {
            let left = pair[0];
            let right = pair[1];
            if left.net_id != right.net_id
                || left.context.depth != right.context.depth
                || left.context.parent_chain_id != right.context.parent_chain_id
            {
                continue;
            }
            if left.context.parent_gap_range == right.context.parent_gap_range {
                continue;
            }
            let between =
                NetRange::new(left.context.fill_range.end, right.context.fill_range.start);
            let index = indexes.get(&left.net_id).ok_or_else(|| {
                Error::Consistency(format!("missing span index for NET {}", left.net_id))
            })?;
            let breaking_chain_id = left.context.parent_chain_id;
            if index.any_overlap_where(between, |span| {
                span.chain_id < broken_chain_id && span.chain_id != breaking_chain_id
            }) {
                continue;
            }

            let suspect_start = left.context.parent_gap_range.end;
            let suspect_end = right.context.parent_gap_range.start;
            if suspect_start >= suspect_end
                || left.context.fill_range.start >= suspect_start
                || left.context.fill_range.end > suspect_start
                || right.context.fill_range.start < suspect_end
                || right.context.fill_range.end <= suspect_end
            {
                return Err(Error::Consistency(format!(
                    "invalid NET break geometry for broken chain {} in {}",
                    broken_chain_id,
                    String::from_utf8_lossy(
                        reference_names
                            .get(&left.net_id)
                            .map(Arc::as_ref)
                            .unwrap_or(b"?")
                    )
                )));
            }

            let candidate = BreakCandidate {
                ordinal: next_ordinal,
                net_id: left.net_id,
                depth: left.context.depth,
                reference_name: reference_names
                    .get(&left.net_id)
                    .cloned()
                    .ok_or_else(|| Error::Consistency("missing NET reference name".into()))?,
                broken_chain_id,
                breaking_chain_id,
                left_fill: left.context.fill_range.start..left.context.fill_range.end,
                right_fill: right.context.fill_range.start..right.context.fill_range.end,
                left_gap: left.context.parent_gap_range.start..left.context.parent_gap_range.end,
                right_gap: right.context.parent_gap_range.start..right.context.parent_gap_range.end,
                suspect: suspect_start..suspect_end,
            };
            next_ordinal = next_ordinal
                .checked_add(1)
                .ok_or_else(|| Error::Overflow("candidate ordinal".into()))?;

            if !candidates_by_breaker.contains_key(&breaking_chain_id) {
                breaker_insertion_order.push(breaking_chain_id);
            }
            candidates_by_breaker
                .entry(breaking_chain_id)
                .or_default()
                .push(candidate);
            interest_chain_ids.insert(broken_chain_id);
            interest_chain_ids.insert(breaking_chain_id);
        }
    }

    let group_order = element_list_order(breaker_insertion_order);
    let groups = group_order
        .into_iter()
        .enumerate()
        .filter_map(|(ordinal, breaking_chain_id)| {
            candidates_by_breaker
                .remove(&breaking_chain_id)
                .map(|candidates| CandidateGroup {
                    compatibility_ordinal: ordinal as u64,
                    breaking_chain_id,
                    candidates,
                })
        })
        .collect();

    Ok(CandidateDiscovery {
        groups,
        interest_chain_ids,
        net_chain_ids,
    })
}
