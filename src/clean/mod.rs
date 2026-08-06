// Copyright (c) 2026 The Hiller Lab at the Senckenberg Gesellschaft für Naturforschung
// Distributed under the terms of the Apache License, Version 2.0.

use std::collections::HashMap;
use std::ops::Range;

use chaintools::{GapCalc, ScoreMatrix};

use crate::candidate::{BreakCandidate, CandidateDiscovery, CandidateGroup};
use crate::chain::{EditableChain, OriginalChain, SubchainView};
use crate::config::CleanerConfig;
use crate::error::{Error, Result};
use crate::output;
use crate::score::{MetricsKey, ParityScorer, SequenceSet, SubchainMetrics};

#[derive(Debug)]
pub struct RemovalEvent {
    pub compatibility_ordinal: u64,
    pub removed_chain: EditableChain,
}

#[derive(Debug)]
struct CandidateEvaluation {
    suspect_metrics: SubchainMetrics,
    ratio_left: f64,
    ratio_right: f64,
    suspect_view: SubchainView,
    breaking_header_score: i64,
    broken_header_score: i64,
    should_remove: bool,
}

pub struct Cleaner {
    config: CleanerConfig,
    discovery: CandidateDiscovery,
    chains: HashMap<u64, EditableChain>,
    originals: Vec<OriginalChain>,
    maximum_chain_id: u64,
    scorer: ParityScorer,
    metrics_cache: HashMap<MetricsKey, SubchainMetrics>,
    events: Vec<RemovalEvent>,
    removed_bed: Vec<String>,
    id_dictionary: Vec<String>,
    suspect_data: Vec<String>,
    suspect_id: u64,
    next_output_ordinal: u64,
    metadata: Vec<Vec<u8>>,
}

impl Cleaner {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        config: CleanerConfig,
        discovery: CandidateDiscovery,
        chains: HashMap<u64, EditableChain>,
        originals: Vec<OriginalChain>,
        input_chain_count: u64,
        maximum_chain_id: u64,
        matrix: ScoreMatrix,
        gaps: GapCalc,
        sequences: SequenceSet,
        metadata: Vec<Vec<u8>>,
    ) -> Self {
        let next_output_ordinal = input_chain_count;
        Self {
            config,
            discovery,
            chains,
            originals,
            maximum_chain_id,
            scorer: ParityScorer::new(matrix, gaps, sequences),
            metrics_cache: HashMap::new(),
            events: Vec::new(),
            removed_bed: Vec::new(),
            id_dictionary: Vec::new(),
            suspect_data: Vec::new(),
            suspect_id: 0,
            next_output_ordinal,
            metadata,
        }
    }

    pub fn run(mut self) -> Result<()> {
        let mut groups = std::mem::take(&mut self.discovery.groups);
        for group in &mut groups {
            self.process_group(group)?;
        }
        self.rescore_dirty_chains()?;
        output::write_outputs(
            &self.config,
            &self.chains,
            &self.originals,
            &self.events,
            &self.removed_bed,
            &self.id_dictionary,
            &self.suspect_data,
            &self.metadata,
        )
    }

    fn process_group(&mut self, group: &mut CandidateGroup) -> Result<()> {
        let mut list = CandidateList::new(group.candidates.len());
        loop {
            loop {
                let mut any_neighbor_updated = false;
                let mut cursor = list.head;
                while let Some(index) = cursor {
                    let following = list.next(index);
                    let previous = list.previous(index);
                    let candidate = group.candidates[index].clone();
                    if let Some(evaluation) = self.evaluate(&candidate, false)?
                        && evaluation.should_remove
                    {
                        self.commit_removal(
                            group.compatibility_ordinal,
                            &candidate,
                            &evaluation,
                            false,
                        )?;
                        if update_neighbors(&mut group.candidates, &candidate, previous, following)
                        {
                            any_neighbor_updated = true;
                        }
                        list.remove(index);
                    }
                    cursor = following;
                }
                if !any_neighbor_updated || list.is_empty() {
                    break;
                }
            }

            if !self.config.do_pairs || self.config.suspect_data_file.is_some() {
                break;
            }

            let mut any_pair_neighbor_updated = false;
            let mut cursor = list.head;
            while let Some(upstream) = cursor {
                let Some(downstream) = list.next(upstream) else {
                    break;
                };
                if !valid_pair(
                    &group.candidates[upstream],
                    &group.candidates[downstream],
                    self.config.max_pair_distance,
                ) {
                    cursor = Some(downstream);
                    continue;
                }
                let pair =
                    pair_candidate(&group.candidates[upstream], &group.candidates[downstream]);
                let before = list.previous(upstream);
                let after = list.next(downstream);
                let evaluation = self.evaluate(&pair, true)?;
                if let Some(evaluation) = evaluation
                    && evaluation.should_remove
                {
                    self.commit_removal(group.compatibility_ordinal, &pair, &evaluation, true)?;
                    if update_neighbors(&mut group.candidates, &pair, before, after) {
                        any_pair_neighbor_updated = true;
                    }
                    list.remove(upstream);
                    list.remove(downstream);
                    cursor = after;
                } else {
                    cursor = Some(downstream);
                }
            }
            if !any_pair_neighbor_updated || list.is_empty() {
                break;
            }
        }
        Ok(())
    }

    fn evaluate(
        &mut self,
        candidate: &BreakCandidate,
        is_pair: bool,
    ) -> Result<Option<CandidateEvaluation>> {
        // Both header scores are read before any subset is scored: a subset
        // that spans a whole chain is that chain, so scoring it overwrites the
        // chain's score in the original program (see `metrics`).
        let breaking_header_score = self
            .chains
            .get(&candidate.breaking_chain_id)
            .expect("validated chain ID")
            .header_score;
        let broken_header_score = self
            .chains
            .get(&candidate.broken_chain_id)
            .expect("validated chain ID")
            .header_score;

        let suspect = self.metrics(candidate.breaking_chain_id, candidate.suspect.clone(), true)?;
        let Some((suspect_view, suspect_metrics)) = suspect else {
            // The original program treats this as a suspect already removed by
            // an earlier overlapping candidate.
            return Ok(None);
        };
        let (_, fill) = self
            .metrics(
                candidate.broken_chain_id,
                candidate.left_fill.start..candidate.right_fill.end,
                false,
            )?
            .expect("required subset errors rather than returning None");
        let (_, left_fill) = self
            .metrics(
                candidate.broken_chain_id,
                candidate.left_fill.start..candidate.suspect.end,
                false,
            )?
            .expect("required subset errors rather than returning None");
        let (_, right_fill) = self
            .metrics(
                candidate.broken_chain_id,
                candidate.suspect.start..candidate.right_fill.end,
                false,
            )?
            .expect("required subset errors rather than returning None");

        let denominator = suspect_metrics.local_score as f64;
        let ratio = fill.global_score as f64 / denominator;
        let ratio_left = left_fill.global_score as f64 / denominator;
        let ratio_right = right_fill.global_score as f64 / denominator;

        if self.config.suspect_data_file.is_some() {
            self.suspect_id = self
                .suspect_id
                .checked_add(1)
                .ok_or_else(|| Error::Overflow("suspect data ID".into()))?;
            self.suspect_data.push(format!(
                "{}\t{}\t{}\t{},{},{},{},{},{},{},{},{},{},{},{},{},{}",
                String::from_utf8_lossy(candidate.reference_name.as_ref()),
                candidate.suspect.start,
                candidate.suspect.end,
                self.suspect_id,
                candidate.breaking_chain_id,
                reported_score(breaking_header_score),
                candidate.broken_chain_id,
                reported_score(broken_header_score),
                reported_score(suspect_metrics.local_score),
                reported_score(fill.global_score),
                reported_score(left_fill.global_score),
                reported_score(right_fill.global_score),
                suspect_metrics.aligned_bases,
                candidate.left_gap.end - candidate.left_gap.start,
                candidate.right_gap.end - candidate.right_gap.start,
                reported_score(left_fill.local_score),
                reported_score(right_fill.local_score)
            ));
        }

        let lr_threshold = if is_pair {
            self.config.lr_fold_threshold_pairs
        } else {
            self.config.lr_fold_threshold
        };
        let zero_rejected = self.config.strict && suspect_metrics.local_score == 0;
        let should_remove = self.config.suspect_data_file.is_none()
            && !zero_rejected
            && ratio_left >= lr_threshold
            && ratio_right >= lr_threshold
            && ratio >= self.config.fold_threshold
            && suspect_metrics.local_score as f64 <= self.config.max_suspect_score
            && suspect_metrics.aligned_bases as f64 <= self.config.max_suspect_bases
            && broken_header_score as f64 >= self.config.min_broken_chain_score
            && candidate.left_gap.end - candidate.left_gap.start >= self.config.min_lr_gap_size
            && candidate.right_gap.end - candidate.right_gap.start >= self.config.min_lr_gap_size;

        Ok(Some(CandidateEvaluation {
            suspect_metrics,
            ratio_left,
            ratio_right,
            suspect_view,
            breaking_header_score,
            broken_header_score,
            should_remove,
        }))
    }

    fn metrics(
        &mut self,
        chain_id: u64,
        range: Range<u32>,
        optional: bool,
    ) -> Result<Option<(SubchainView, SubchainMetrics)>> {
        let chain = self
            .chains
            .get(&chain_id)
            .ok_or_else(|| Error::Consistency(format!("missing chain ID {chain_id}")))?;
        let key = MetricsKey {
            chain_id,
            revision: chain.revision,
            start: range.start,
            end: range.end,
        };
        let Some(view) = chain.subset_on_reference(range.start, range.end) else {
            if optional {
                return Ok(None);
            }
            return Err(Error::Consistency(format!(
                "chain {chain_id} has no blocks in required reference subset {}-{}",
                range.start, range.end
            )));
        };
        if let Some(metrics) = self.metrics_cache.get(&key).copied() {
            return Ok(Some((view, metrics)));
        }
        let spans_whole_chain =
            range.start <= chain.reference_start && range.end >= chain.reference_end;
        let metrics = self.scorer.score(chain, &view)?;
        self.metrics_cache.insert(key, metrics);
        if spans_whole_chain {
            // `chainSubsetOnT` returns the original chain when the requested
            // interval covers it, and `getChainScore` then writes the freshly
            // computed global score into that chain's header. Later threshold
            // tests, BED records, and the written chain all see the new value.
            self.chains
                .get_mut(&chain_id)
                .expect("chain exists")
                .header_score = metrics.global_score;
        }
        Ok(Some((view, metrics)))
    }

    fn commit_removal(
        &mut self,
        compatibility_ordinal: u64,
        candidate: &BreakCandidate,
        evaluation: &CandidateEvaluation,
        is_pair: bool,
    ) -> Result<()> {
        let new_chain_id = self
            .maximum_chain_id
            .checked_add(1)
            .ok_or_else(|| Error::Overflow("new suspect chain ID".into()))?;
        let removed_chain = self
            .chains
            .get(&candidate.breaking_chain_id)
            .expect("validated chain ID")
            .owned_from_blocks(
                new_chain_id,
                evaluation.suspect_metrics.global_score,
                evaluation.suspect_view.blocks.clone(),
                self.next_output_ordinal,
            )?;

        self.chains
            .get_mut(&candidate.breaking_chain_id)
            .expect("validated chain ID")
            .remove_internal_blocks(candidate.suspect.clone())?;
        self.maximum_chain_id = new_chain_id;
        self.next_output_ordinal = self
            .next_output_ordinal
            .checked_add(1)
            .ok_or_else(|| Error::Overflow("output ordinal".into()))?;

        self.removed_bed.push(format!(
            "{}\t{}\t{}\tbreakingChainID_{}_Score_{}_brokenChainID_{}_Score_{}_suspectLocalScore_{}_RatioL_{:.2}_RatioR_{:.2}\t1000\t+\t{}\t{}\t{}",
            String::from_utf8_lossy(candidate.reference_name.as_ref()),
            candidate.suspect.start,
            candidate.suspect.end,
            candidate.breaking_chain_id,
            reported_score(evaluation.breaking_header_score),
            candidate.broken_chain_id,
            reported_score(evaluation.broken_header_score),
            reported_score(evaluation.suspect_metrics.local_score),
            evaluation.ratio_left,
            evaluation.ratio_right,
            candidate.suspect.start,
            candidate.suspect.end,
            if is_pair { "0,100,255" } else { "0,0,153" }
        ));
        self.id_dictionary
            .push(format!("{}\t{}", new_chain_id, candidate.breaking_chain_id));
        self.events.push(RemovalEvent {
            compatibility_ordinal,
            removed_chain,
        });
        Ok(())
    }

    fn rescore_dirty_chains(&mut self) -> Result<()> {
        let dirty: Vec<(u64, Range<u32>)> = self
            .chains
            .values()
            .filter(|chain| chain.dirty)
            .map(|chain| (chain.id, chain.reference_start..chain.reference_end))
            .collect();
        for (id, range) in dirty {
            let (_, metrics) = self
                .metrics(id, range, false)?
                .expect("required full chain subset");
            self.chains.get_mut(&id).expect("chain exists").final_score =
                Some(metrics.global_score);
        }
        Ok(())
    }
}

/// Format a score the way the original program's `%d` of a `(int)`-cast double
/// does.
///
/// Chain scores are doubles in C and are printed after a cast to `int`, so any
/// score outside the 32-bit range is converted out of range. On the reference
/// x86-64 build that conversion yields `INT_MIN`, which whole-genome chains with
/// scores above 2^31-1 do reach. BED and suspect-data records reproduce it;
/// chain records themselves are written from the full value, as `chainWrite`
/// prints the double directly.
fn reported_score(score: i64) -> i32 {
    i32::try_from(score).unwrap_or(i32::MIN)
}

/// Index-based doubly linked list over one group's candidates.
///
/// The original program keeps the remaining breaks of a breaking chain in a
/// doubly linked list and passes `prev`/`next` of the current break into the
/// removal test. This arena reproduces that structure and its neighbor
/// semantics without per-candidate allocation or repeated scans.
struct CandidateList {
    head: Option<usize>,
    previous: Vec<Option<usize>>,
    next: Vec<Option<usize>>,
    live: usize,
}

impl CandidateList {
    fn new(len: usize) -> Self {
        Self {
            head: (len > 0).then_some(0),
            previous: (0..len).map(|index| index.checked_sub(1)).collect(),
            next: (0..len)
                .map(|index| (index + 1 < len).then_some(index + 1))
                .collect(),
            live: len,
        }
    }

    fn previous(&self, index: usize) -> Option<usize> {
        self.previous[index]
    }

    fn next(&self, index: usize) -> Option<usize> {
        self.next[index]
    }

    fn is_empty(&self) -> bool {
        self.live == 0
    }

    fn remove(&mut self, index: usize) {
        let previous = self.previous[index];
        let next = self.next[index];
        match previous {
            Some(previous) => self.next[previous] = next,
            None => self.head = next,
        }
        if let Some(next) = next {
            self.previous[next] = previous;
        }
        self.live -= 1;
    }
}

fn update_neighbors(
    candidates: &mut [BreakCandidate],
    current: &BreakCandidate,
    previous: Option<usize>,
    next: Option<usize>,
) -> bool {
    let mut updated = false;
    if let Some(index) = previous {
        let upstream = &mut candidates[index];
        if upstream.broken_chain_id == current.broken_chain_id
            && upstream.breaking_chain_id == current.breaking_chain_id
            && upstream.right_fill == current.left_fill
        {
            upstream.right_fill.end = current.right_fill.end;
            upstream.right_gap.end = current.right_gap.end;
            updated = true;
        }
    }
    if let Some(index) = next {
        let downstream = &mut candidates[index];
        if downstream.broken_chain_id == current.broken_chain_id
            && downstream.breaking_chain_id == current.breaking_chain_id
            && downstream.left_fill == current.right_fill
        {
            downstream.left_fill.start = current.left_fill.start;
            downstream.left_gap.start = current.left_gap.start;
            updated = true;
        }
    }
    updated
}

fn valid_pair(upstream: &BreakCandidate, downstream: &BreakCandidate, max_distance: u32) -> bool {
    upstream.broken_chain_id == downstream.broken_chain_id
        && upstream.breaking_chain_id == downstream.breaking_chain_id
        && upstream.depth == downstream.depth
        && i64::from(downstream.suspect.start) - i64::from(upstream.suspect.end)
            <= i64::from(max_distance)
        && upstream.right_gap == downstream.left_gap
}

fn pair_candidate(upstream: &BreakCandidate, downstream: &BreakCandidate) -> BreakCandidate {
    BreakCandidate {
        ordinal: upstream.ordinal,
        net_id: upstream.net_id,
        depth: upstream.depth,
        reference_name: upstream.reference_name.clone(),
        broken_chain_id: upstream.broken_chain_id,
        breaking_chain_id: upstream.breaking_chain_id,
        left_fill: upstream.left_fill.clone(),
        right_fill: downstream.right_fill.clone(),
        left_gap: upstream.left_gap.clone(),
        right_gap: downstream.right_gap.clone(),
        suspect: upstream.left_gap.end..downstream.right_gap.start,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn candidate(left: u32, right: u32) -> BreakCandidate {
        BreakCandidate {
            ordinal: 0,
            net_id: 0,
            depth: 2,
            reference_name: std::sync::Arc::from(&b"chr1"[..]),
            broken_chain_id: 2,
            breaking_chain_id: 1,
            left_fill: left - 20..left - 10,
            right_fill: right + 10..right + 20,
            left_gap: left - 10..left,
            right_gap: right..right + 10,
            suspect: left..right,
        }
    }

    #[test]
    fn scores_outside_the_32_bit_range_are_reported_like_c() {
        assert_eq!(reported_score(100_000), 100_000);
        assert_eq!(reported_score(i64::from(i32::MAX)), i32::MAX);
        assert_eq!(reported_score(i64::from(i32::MAX) + 1), i32::MIN);
        assert_eq!(reported_score(3_592_426_336), i32::MIN);
    }

    #[test]
    fn pair_requires_shared_middle_gap() {
        let first = candidate(20, 30);
        let mut second = candidate(40, 50);
        second.left_gap = first.right_gap.clone();
        assert!(valid_pair(&first, &second, 10));
        second.left_gap.start += 1;
        assert!(!valid_pair(&first, &second, 10));
    }
}
