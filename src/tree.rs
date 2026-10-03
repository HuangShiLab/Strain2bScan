//! Layer-1: Cluster Search Tree descent (StrainScan port).

use crate::cst::Cst;
use crate::db::StrainDb;
use crate::depth::{detectable_fraction, marker_panel_evidence, support_count, PanelStats};
use crate::identify::{Layer1, Params};
use crate::markers::{Marker, MarkerCounts};

/// Minimum markers for a node's own set to be worth testing — and, crucially, to be trusted to
/// rule its whole subtree OUT.
pub const MIN_NODE_MARKERS: usize = 30;

/// Widest clade the internal-node fallback may report.
pub const MAX_FALLBACK_CLADE: usize = 8;

/// What the tree stored in a database can actually do for it.
#[derive(Debug, Clone, Copy)]
pub struct TreeUtility {
    /// Clusters whose own unique markers fall below the support floor.
    pub rescuable: usize,
    /// Internal nodes carrying enough markers to be tested.
    pub informative_nodes: usize,
    /// Internal nodes where both children are themselves informative.
    pub fallback_nodes: usize,
}

impl TreeUtility {
    /// Whether the tree is worth descending.
    pub fn worth_descending(&self) -> bool {
        (self.rescuable > 0 && self.informative_nodes > 0) || self.fallback_nodes > 0
    }
}

/// Measure what the stored tree can do. `None` when the database carries no tree.
pub fn tree_utility(db: &StrainDb, support_floor: usize) -> Option<TreeUtility> {
    let tree = db.tree.as_ref()?;
    let rescuable = (0..db.n_strains())
        .filter(|&j| db.unique_marker_count(j) < support_floor)
        .count();
    let informative = |v: usize| {
        tree.node_markers[v]
            .iter()
            .filter(|&&m| db.is_quantifiable(m))
            .count()
            >= MIN_NODE_MARKERS
    };
    let informative_nodes = (0..tree.n_nodes())
        .filter(|&v| !tree.is_leaf(v) && informative(v))
        .count();
    let fallback_nodes = (0..tree.n_nodes())
        .filter(|&v| {
            !tree.is_leaf(v)
                && informative(v)
                && tree.children[v].is_some_and(|(a, b)| informative(a) && informative(b))
        })
        .count();
    Some(TreeUtility {
        rescuable,
        informative_nodes,
        fallback_nodes,
    })
}

/// Resolve [`Layer1::Auto`] against the database. Explicit choices pass through untouched.
pub fn resolve_layer1(db: &StrainDb, p: &Params) -> Layer1 {
    match p.layer1 {
        Layer1::Auto => match tree_utility(db, p.min_support_markers) {
            Some(u) if u.worth_descending() => Layer1::Cst,
            _ => Layer1::Unique,
        },
        explicit => explicit,
    }
}

/// One leaf accepted by the tree descent.
#[derive(Debug, Clone)]
pub struct TreeCall {
    /// The node the descent resolved to.
    pub node: usize,
    /// Leaf ids beneath `node`.
    pub desc_leaves: Vec<usize>,
    /// Markers pooled along the unique path.
    pub panel: usize,
    pub detected: usize,
    pub coverage: f64,
    pub depth: f64,
    /// Nodes whose markers were pooled.
    pub path: Vec<usize>,
}

/// Evidence for one marker set.
pub fn set_evidence(markers: &[Marker], counts: &MarkerCounts) -> PanelStats {
    marker_panel_evidence(markers, counts)
}

/// Descend the Cluster Search Tree with no cross-species restriction.
pub fn descend_tree(cst: &Cst, counts: &MarkerCounts, p: &Params) -> Vec<TreeCall> {
    let node_markers: Vec<Vec<Marker>> = cst
        .node_markers
        .iter()
        .map(|set| set.iter().copied().collect())
        .collect();
    descend_tree_inner(cst, counts, p, &node_markers)
}

/// The tree descent, honouring `multi-profile`'s cross-species marker restriction.
pub fn descend_tree_masked(
    cst: &Cst,
    counts: &MarkerCounts,
    p: &Params,
    masked_nodes: Option<&[Vec<Marker>]>,
) -> Vec<TreeCall> {
    let owned: Vec<Vec<Marker>>;
    let node_markers: &[Vec<Marker>] = match masked_nodes {
        Some(nodes) => nodes,
        None => {
            owned = cst
                .node_markers
                .iter()
                .map(|set| set.iter().copied().collect())
                .collect();
            &owned
        }
    };
    descend_tree_inner(cst, counts, p, node_markers)
}

pub fn descend_tree_inner(
    cst: &Cst,
    counts: &MarkerCounts,
    p: &Params,
    node_markers: &[Vec<Marker>],
) -> Vec<TreeCall> {
    if cst.n_leaves() == 0 {
        return Vec::new();
    }
    if cst.n_leaves() == 1 {
        let ms = &node_markers[0];
        let ev = set_evidence(ms, counts);
        if ev.panel > 0 {
            let support = support_count(&ev, p.adaptive_singleton);
            let coverage = ev.detected1 as f64 / ev.panel as f64;
            if support >= p.min_support_markers && coverage >= p.min_coverage {
                return vec![TreeCall {
                    node: 0,
                    desc_leaves: vec![0],
                    panel: ev.panel,
                    detected: support,
                    coverage,
                    depth: ev.depth,
                    path: vec![0],
                }];
            }
        }
        return Vec::new();
    }

    let fires = |v: usize| -> bool {
        let ms = &node_markers[v];
        if ms.len() < MIN_NODE_MARKERS {
            return false;
        }
        let ev = set_evidence(ms, counts);
        let support = support_count(&ev, p.adaptive_singleton);
        let coverage = ev.detected1 as f64 / ev.panel as f64;
        support >= p.min_support_markers && coverage >= p.min_coverage
    };
    let informative = |v: usize| node_markers[v].len() >= MIN_NODE_MARKERS;

    let mut entered: Vec<bool> = vec![false; cst.n_nodes()];
    let mut reached: Vec<usize> = Vec::new();
    let mut stack = vec![cst.root];
    entered[cst.root] = true;
    while let Some(v) = stack.pop() {
        if cst.is_leaf(v) {
            reached.push(v);
            continue;
        }
        let (a, b) = cst.children[v].expect("internal node has children");
        let mut descended = false;
        for c in [a, b] {
            if !informative(c) || fires(c) {
                entered[c] = true;
                stack.push(c);
                descended = true;
            }
        }
        if !descended {
            reached.push(v);
        }
    }

    let mut out: Vec<TreeCall> = Vec::new();
    let mut rejected: Vec<usize> = Vec::new();
    for &node in &reached {
        let mut path = vec![node];
        let mut pooled: Vec<Marker> = node_markers[node].clone();
        let mut v = node;
        while let Some(par) = cst.parent[v] {
            match cst.sibling(v) {
                Some(s) if !entered[s] => {
                    pooled.extend(node_markers[par].iter().copied());
                    path.push(par);
                    v = par;
                }
                _ => break,
            }
        }
        let ev = set_evidence(&pooled, counts);
        if ev.panel == 0 {
            continue;
        }
        let support = support_count(&ev, p.adaptive_singleton);
        if support < p.min_support_markers {
            rejected.push(node);
            continue;
        }
        let coverage = ev.detected1 as f64 / ev.panel as f64;
        if coverage < p.min_coverage {
            rejected.push(node);
            continue;
        }
        let expected = detectable_fraction(ev.depth);
        if expected > 0.0 && coverage / expected < p.min_consistency {
            rejected.push(node);
            continue;
        }
        out.push(TreeCall {
            node,
            desc_leaves: cst.desc_leaves[node].clone(),
            panel: ev.panel,
            detected: support,
            coverage,
            depth: ev.depth,
            path,
        });
    }

    if !rejected.is_empty() {
        let accepted_under = |v: usize, out: &[TreeCall]| -> bool {
            out.iter()
                .any(|c| cst.desc_leaves[v].contains(&cst.desc_leaves[c.node][0]))
        };
        let mut added: Vec<usize> = Vec::new();
        for &r in &rejected {
            let mut v = r;
            while let Some(par) = cst.parent[v] {
                if accepted_under(par, &out) || added.contains(&par) {
                    break;
                }
                if cst.desc_leaves[par].len() > MAX_FALLBACK_CLADE {
                    break;
                }
                let ms = &node_markers[par];
                if ms.len() >= MIN_NODE_MARKERS {
                    let ev = set_evidence(ms, counts);
                    let support = support_count(&ev, p.adaptive_singleton);
                    let coverage = ev.detected1 as f64 / ev.panel as f64;
                    let expected = detectable_fraction(ev.depth);
                    let consistent = expected <= 0.0 || coverage / expected >= p.min_consistency;
                    if support >= p.min_support_markers && coverage >= p.min_coverage && consistent
                    {
                        out.push(TreeCall {
                            node: par,
                            desc_leaves: cst.desc_leaves[par].clone(),
                            panel: ev.panel,
                            detected: support,
                            coverage,
                            depth: ev.depth,
                            path: vec![par],
                        });
                        added.push(par);
                        break;
                    }
                }
                v = par;
            }
        }
    }
    out
}
