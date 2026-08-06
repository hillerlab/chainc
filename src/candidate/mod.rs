// Copyright (c) 2026 The Hiller Lab at the Senckenberg Gesellschaft für Naturforschung
// Distributed under the terms of the Apache License, Version 2.0.

mod discover;
mod interval_index;
mod kent_order;

pub use discover::{BreakCandidate, CandidateDiscovery, CandidateGroup, discover_candidates};
