//! strain2bscan — prototype 2bRAD-tag strain-level profiler.
//!
//! Architecture (decided 2026-06): **Fast2bRAD-M core + ported StrainScan Layer-2**.
//!
//! Pipeline:
//!   genomes ──digest──▶ markers ──build──▶ StrainDb (sparse, unique-marker aware)
//!   reads   ──digest──▶ marker counts ──profile (Layer-2)──▶ strain calls + abundance
//!
//! The key accuracy idea ported from real StrainScan and absent from `strainscan-rust`:
//! score candidate strains on **unique / not-yet-covered** markers with a per-strain
//! support threshold, then estimate abundance with a **non-negative** Elastic Net and
//! filter by coverage + relative abundance. 2bRAD tags are a sparse, taxonomy-specific
//! marker set, so they supply exactly the low-redundancy markers this algorithm needs.

pub mod bench;
pub mod cst;
pub mod db;
pub mod enzymes;
pub mod identify;
pub mod iibdb;
pub mod markers;
pub mod parallel;

pub use markers::Marker;
