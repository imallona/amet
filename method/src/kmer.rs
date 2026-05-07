//! Per-cell L-mer counting at fixed lag.
//!
//! For a feature, the cell's calls are placed onto the reference's CpG positions, leaving
//! gaps where the cell has no call. A pair (X_i, X_{i+k}) is counted only when both
//! reference positions i and i+k inside the feature are observed in this cell.

use crate::MethCall;
use crate::features::Feature;
use crate::reference::CpgReference;

/// One feature's per-cell binary observation: one entry per reference CpG in the feature,
/// either Some(0/1) if observed in this cell or None if missing.
pub struct CellWindow<'a> {
    pub feature: &'a Feature,
    pub calls: Vec<Option<u8>>,
}

impl<'a> CellWindow<'a> {
    pub fn n_observed(&self) -> usize {
        self.calls.iter().filter(|c| c.is_some()).count()
    }

    pub fn mean_meth(&self) -> Option<f64> {
        let mut sum = 0u32;
        let mut n = 0u32;
        for v in self.calls.iter().flatten() {
            sum += *v as u32;
            n += 1;
        }
        if n == 0 {
            None
        } else {
            Some(sum as f64 / n as f64)
        }
    }
}

/// Counts of (X_i, X_{i+k}) pairs across all valid positions in a feature for one cell.
/// Index: 2*x_i + x_{i+k}.
#[derive(Debug, Clone, Copy, Default)]
pub struct PairCounts {
    pub counts: [u32; 4],
}

impl PairCounts {
    pub fn total(&self) -> u32 {
        self.counts.iter().sum()
    }

    pub fn add(&mut self, other: &PairCounts) {
        for i in 0..4 {
            self.counts[i] += other.counts[i];
        }
    }
}

/// Marginal counts of single positions used for I_k. Stored as [count_of_0, count_of_1].
#[derive(Debug, Clone, Copy, Default)]
pub struct MarginalCounts {
    pub counts: [u32; 2],
}

impl MarginalCounts {
    pub fn total(&self) -> u32 {
        self.counts[0] + self.counts[1]
    }

    pub fn add(&mut self, other: &MarginalCounts) {
        self.counts[0] += other.counts[0];
        self.counts[1] += other.counts[1];
    }
}

/// Build the per-feature observation vector for one cell.
///
/// `calls` must be sorted by (chrom_id, pos). The function aligns them to the reference
/// positions in the feature's `cpg_start_idx..cpg_end_idx` range.
pub fn build_window<'a>(
    feature: &'a Feature,
    reference: &CpgReference,
    calls: &[MethCall],
    threshold: f64,
    min_reads: u32,
) -> CellWindow<'a> {
    let positions = &reference.positions[feature.chrom_id as usize];
    let n = feature.cpg_end_idx - feature.cpg_start_idx;
    let mut window = vec![None; n];

    if n == 0 {
        return CellWindow {
            feature,
            calls: window,
        };
    }

    let feature_start_pos = positions[feature.cpg_start_idx];
    let feature_end_pos = positions[feature.cpg_end_idx - 1];

    let lo = calls.partition_point(|c| (c.chrom_id, c.pos) < (feature.chrom_id, feature_start_pos));
    let hi = calls.partition_point(|c| (c.chrom_id, c.pos) <= (feature.chrom_id, feature_end_pos));

    let mut ref_idx = 0;
    for call in &calls[lo..hi] {
        if call.chrom_id != feature.chrom_id {
            continue;
        }
        while ref_idx < n && positions[feature.cpg_start_idx + ref_idx] < call.pos {
            ref_idx += 1;
        }
        if ref_idx >= n {
            break;
        }
        if positions[feature.cpg_start_idx + ref_idx] == call.pos {
            if let Some(b) = call.binarize(threshold, min_reads) {
                window[ref_idx] = Some(b);
            }
            ref_idx += 1;
        }
    }

    CellWindow {
        feature,
        calls: window,
    }
}

/// Count (X_i, X_{i+k}) pairs in a cell window for a single lag k.
pub fn pair_counts(window: &CellWindow, lag: usize) -> PairCounts {
    let mut pc = PairCounts::default();
    let n = window.calls.len();
    if lag == 0 || lag >= n {
        return pc;
    }
    for i in 0..(n - lag) {
        if let (Some(a), Some(b)) = (window.calls[i], window.calls[i + lag]) {
            pc.counts[(a as usize) * 2 + b as usize] += 1;
        }
    }
    pc
}

/// Marginal counts over the whole window.
pub fn marginal_counts(window: &CellWindow) -> MarginalCounts {
    let mut mc = MarginalCounts::default();
    for v in window.calls.iter().flatten() {
        mc.counts[*v as usize] += 1;
    }
    mc
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ref3() -> CpgReference {
        CpgReference {
            chrom_names: vec!["chr1".into()],
            chrom_id_of: [("chr1".into(), 0u32)].into_iter().collect(),
            positions: vec![vec![10, 20, 30, 40, 50]],
        }
    }

    fn feat_full() -> Feature {
        Feature {
            feature_id: "f".into(),
            chrom_id: 0,
            start: 0,
            end: 100,
            cpg_start_idx: 0,
            cpg_end_idx: 5,
        }
    }

    #[test]
    fn build_window_aligns_calls_to_reference_positions() {
        let r = ref3();
        let f = feat_full();
        // Cell observes CpGs 10, 30, 50 only; 20 and 40 are gaps.
        let calls = vec![
            MethCall {
                chrom_id: 0,
                pos: 10,
                m: 1,
                t: 1,
            },
            MethCall {
                chrom_id: 0,
                pos: 30,
                m: 0,
                t: 1,
            },
            MethCall {
                chrom_id: 0,
                pos: 50,
                m: 1,
                t: 1,
            },
        ];
        let w = build_window(&f, &r, &calls, 0.0, 1);
        assert_eq!(w.calls, vec![Some(1), None, Some(0), None, Some(1)]);
    }

    #[test]
    fn pair_counts_lag1_skips_gaps() {
        let r = ref3();
        let f = feat_full();
        let calls = vec![
            MethCall {
                chrom_id: 0,
                pos: 10,
                m: 1,
                t: 1,
            },
            MethCall {
                chrom_id: 0,
                pos: 20,
                m: 1,
                t: 1,
            },
            // 30 missing
            MethCall {
                chrom_id: 0,
                pos: 40,
                m: 0,
                t: 1,
            },
            MethCall {
                chrom_id: 0,
                pos: 50,
                m: 0,
                t: 1,
            },
        ];
        let w = build_window(&f, &r, &calls, 0.0, 1);
        // calls = [1,1,_,0,0]; lag-1 pairs: (1,1), (0,0). (1,_) and (_,0) and (0,0) excluded.
        let pc = pair_counts(&w, 1);
        // 11 (idx 3) = 1; 00 (idx 0) = 1.
        assert_eq!(pc.counts[3], 1);
        assert_eq!(pc.counts[0], 1);
        assert_eq!(pc.total(), 2);
    }

    #[test]
    fn pair_counts_lag2_jumps_over_gap() {
        let r = ref3();
        let f = feat_full();
        let calls = vec![
            MethCall {
                chrom_id: 0,
                pos: 10,
                m: 1,
                t: 1,
            },
            MethCall {
                chrom_id: 0,
                pos: 20,
                m: 0,
                t: 1,
            },
            MethCall {
                chrom_id: 0,
                pos: 30,
                m: 1,
                t: 1,
            },
            MethCall {
                chrom_id: 0,
                pos: 40,
                m: 0,
                t: 1,
            },
            MethCall {
                chrom_id: 0,
                pos: 50,
                m: 1,
                t: 1,
            },
        ];
        let w = build_window(&f, &r, &calls, 0.0, 1);
        // calls = [1,0,1,0,1]; lag-2 pairs: (1,1), (0,0), (1,1).
        let pc = pair_counts(&w, 2);
        assert_eq!(pc.counts[3], 2); // (1,1)
        assert_eq!(pc.counts[0], 1); // (0,0)
    }

    #[test]
    fn marginal_counts_balanced() {
        let r = ref3();
        let f = feat_full();
        let calls = vec![
            MethCall {
                chrom_id: 0,
                pos: 10,
                m: 1,
                t: 1,
            },
            MethCall {
                chrom_id: 0,
                pos: 20,
                m: 0,
                t: 1,
            },
            MethCall {
                chrom_id: 0,
                pos: 30,
                m: 1,
                t: 1,
            },
        ];
        let w = build_window(&f, &r, &calls, 0.0, 1);
        let mc = marginal_counts(&w);
        assert_eq!(mc.counts, [1, 2]);
    }

    #[test]
    fn n_observed_correct() {
        let r = ref3();
        let f = feat_full();
        let calls = vec![
            MethCall {
                chrom_id: 0,
                pos: 10,
                m: 1,
                t: 1,
            },
            MethCall {
                chrom_id: 0,
                pos: 30,
                m: 0,
                t: 1,
            },
        ];
        let w = build_window(&f, &r, &calls, 0.0, 1);
        assert_eq!(w.n_observed(), 2);
    }

    #[test]
    fn min_reads_filter_drops_low_coverage() {
        let r = ref3();
        let f = feat_full();
        let calls = vec![
            MethCall {
                chrom_id: 0,
                pos: 10,
                m: 1,
                t: 1,
            },
            MethCall {
                chrom_id: 0,
                pos: 20,
                m: 5,
                t: 5,
            },
        ];
        let w = build_window(&f, &r, &calls, 0.0, 2);
        assert_eq!(w.calls, vec![None, Some(1), None, None, None]);
    }

    #[test]
    fn lag_too_large_yields_zero_pairs() {
        let r = ref3();
        let f = feat_full();
        let calls = vec![
            MethCall {
                chrom_id: 0,
                pos: 10,
                m: 1,
                t: 1,
            },
            MethCall {
                chrom_id: 0,
                pos: 20,
                m: 1,
                t: 1,
            },
        ];
        let w = build_window(&f, &r, &calls, 0.0, 1);
        let pc = pair_counts(&w, 10);
        assert_eq!(pc.total(), 0);
    }
}
