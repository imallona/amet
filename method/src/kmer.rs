//! Per-cell pair counting at fixed lag.
//!
//! For a feature, the cell's calls are placed onto the reference's CpG positions. Only
//! observed positions are kept in a compact list, sorted by reference-CpG index. A pair
//! (X_i, X_{i+k}) at lag k is counted only when both reference positions i and i+k inside
//! the feature are observed in this cell.

use crate::MethCall;
use crate::features::Feature;
use crate::reference::CpgReference;

/// One observed methylation call inside a feature, addressed by the reference-CpG index
/// and the genomic position of the C on the + strand.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Observation {
    pub ref_idx: u32,
    pub value: u8,
    pub pos: u64,
}

/// Per-cell observations for a feature, sorted by `ref_idx` ascending.
pub struct CellWindow<'a> {
    pub feature: &'a Feature,
    pub observed: Vec<Observation>,
}

impl<'a> CellWindow<'a> {
    pub fn n_observed(&self) -> usize {
        self.observed.len()
    }

    pub fn mean_meth(&self) -> Option<f64> {
        if self.observed.is_empty() {
            return None;
        }
        let sum: u32 = self.observed.iter().map(|o| o.value as u32).sum();
        Some(sum as f64 / self.observed.len() as f64)
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

/// Build the per-feature observation list for one cell.
///
/// `calls` must be sorted by (chrom_id, pos). The function aligns them to the reference
/// positions in the feature's `cpg_start_idx..cpg_end_idx` range and emits one entry per
/// observed CpG that passes the binarization filter.
pub fn build_window<'a>(
    feature: &'a Feature,
    reference: &'a CpgReference,
    calls: &[MethCall],
    threshold: f64,
    min_reads: u32,
) -> CellWindow<'a> {
    let positions = &reference.positions[feature.chrom_id as usize];
    let feature_positions = &positions[feature.cpg_start_idx..feature.cpg_end_idx];
    let n = feature_positions.len();

    if n == 0 {
        return CellWindow {
            feature,
            observed: Vec::new(),
        };
    }

    let feature_start_pos = feature_positions[0];
    let feature_end_pos = feature_positions[n - 1];

    let lo = calls.partition_point(|c| (c.chrom_id, c.pos) < (feature.chrom_id, feature_start_pos));
    let hi = calls.partition_point(|c| (c.chrom_id, c.pos) <= (feature.chrom_id, feature_end_pos));

    let mut observed = Vec::new();
    let mut ref_idx = 0usize;
    for call in &calls[lo..hi] {
        if call.chrom_id != feature.chrom_id {
            continue;
        }
        while ref_idx < n && feature_positions[ref_idx] < call.pos {
            ref_idx += 1;
        }
        if ref_idx >= n {
            break;
        }
        if feature_positions[ref_idx] == call.pos {
            if let Some(b) = call.binarize(threshold, min_reads) {
                observed.push(Observation {
                    ref_idx: ref_idx as u32,
                    value: b,
                    pos: call.pos,
                });
            }
            ref_idx += 1;
        }
    }

    CellWindow { feature, observed }
}

/// Count (X_i, X_{i+k}) pairs at a single lag k. Pairs whose genomic distance exceeds
/// `max_distance` are skipped; pass 0 to disable the cap.
pub fn pair_counts(window: &CellWindow, lag: usize, max_distance: u64) -> PairCounts {
    let mut pc = PairCounts::default();
    let obs = &window.observed;
    if lag == 0 || obs.is_empty() {
        return pc;
    }
    let lag_u32 = lag as u32;
    let mut j = 0usize;
    for i in 0..obs.len() {
        let target = obs[i].ref_idx + lag_u32;
        while j < obs.len() && obs[j].ref_idx < target {
            j += 1;
        }
        if j >= obs.len() {
            break;
        }
        if obs[j].ref_idx == target
            && (max_distance == 0 || obs[j].pos - obs[i].pos <= max_distance)
        {
            let idx = (obs[i].value as usize) * 2 + obs[j].value as usize;
            pc.counts[idx] += 1;
        }
    }
    pc
}

/// Count pairs at every lag in 1..=k_max in a single sweep over observed pairs.
/// Breaks the inner walk as soon as the lag exceeds k_max or the genomic distance
/// exceeds `max_distance` (0 disables the cap). Returns a `Vec` of length `k_max`,
/// where index `k - 1` holds the table for lag k.
///
/// Cost is O(n_observed * local_neighbours_within_k_max), independent of the feature's
/// total CpG count. This is the hot path for large features such as heterochromatin
/// blocks where most reference CpGs are not observed in a given single cell.
pub fn pair_counts_all_lags(
    window: &CellWindow,
    k_max: usize,
    max_distance: u64,
) -> Vec<PairCounts> {
    let mut out = vec![PairCounts::default(); k_max];
    if k_max == 0 {
        return out;
    }
    let obs = &window.observed;
    let k_max_u32 = k_max as u32;
    for i in 0..obs.len() {
        let oi = obs[i];
        for oj in &obs[i + 1..] {
            let lag = oj.ref_idx - oi.ref_idx;
            if lag > k_max_u32 {
                break;
            }
            if max_distance > 0 && oj.pos - oi.pos > max_distance {
                break;
            }
            let idx = (oi.value as usize) * 2 + oj.value as usize;
            out[(lag as usize) - 1].counts[idx] += 1;
        }
    }
    out
}

/// Marginal counts over all observed positions in the window.
pub fn marginal_counts(window: &CellWindow) -> MarginalCounts {
    let mut mc = MarginalCounts::default();
    for o in &window.observed {
        mc.counts[o.value as usize] += 1;
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
        assert_eq!(
            w.observed,
            vec![
                Observation {
                    ref_idx: 0,
                    value: 1,
                    pos: 10
                },
                Observation {
                    ref_idx: 2,
                    value: 0,
                    pos: 30
                },
                Observation {
                    ref_idx: 4,
                    value: 1,
                    pos: 50
                },
            ]
        );
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
        // observed = [1@0, 1@1, 0@3, 0@4]; lag-1 pairs: (1,1) ref 0->1, (0,0) ref 3->4.
        let pc = pair_counts(&w, 1, 0);
        assert_eq!(pc.counts[3], 1); // (1,1)
        assert_eq!(pc.counts[0], 1); // (0,0)
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
        // observed values [1,0,1,0,1]; lag-2 pairs: (1,1), (0,0), (1,1).
        let pc = pair_counts(&w, 2, 0);
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
        // pos 10 dropped (t < 2); pos 20 kept (value 1, ref_idx 1).
        assert_eq!(w.observed.len(), 1);
        assert_eq!(
            w.observed[0],
            Observation {
                ref_idx: 1,
                value: 1,
                pos: 20
            }
        );
    }

    #[test]
    fn max_distance_drops_far_pairs() {
        let r = CpgReference {
            chrom_names: vec!["chr1".into()],
            chrom_id_of: [("chr1".into(), 0u32)].into_iter().collect(),
            positions: vec![vec![10, 20, 1100]],
        };
        let f = Feature {
            feature_id: "f".into(),
            chrom_id: 0,
            start: 0,
            end: 2000,
            cpg_start_idx: 0,
            cpg_end_idx: 3,
        };
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
                pos: 1100,
                m: 1,
                t: 1,
            },
        ];
        let w = build_window(&f, &r, &calls, 0.0, 1);
        // lag-1: (10,20) distance 10, (20,1100) distance 1080.
        // With cap=1000, only the (10,20) pair survives. Calls are (1,0) so index 2.
        let pc = pair_counts(&w, 1, 1000);
        assert_eq!(pc.counts[2], 1);
        assert_eq!(pc.total(), 1);
        // lag-2: (10,1100) distance 1090, dropped.
        let pc2 = pair_counts(&w, 2, 1000);
        assert_eq!(pc2.total(), 0);
        // Disabled cap keeps all reachable pairs.
        let pc_off = pair_counts(&w, 1, 0);
        assert_eq!(pc_off.total(), 2);
    }

    #[test]
    fn max_distance_boundary_inclusive() {
        // Distance exactly equal to the cap is kept; one over is dropped.
        let r = CpgReference {
            chrom_names: vec!["chr1".into()],
            chrom_id_of: [("chr1".into(), 0u32)].into_iter().collect(),
            positions: vec![vec![0, 100]],
        };
        let f = Feature {
            feature_id: "f".into(),
            chrom_id: 0,
            start: 0,
            end: 200,
            cpg_start_idx: 0,
            cpg_end_idx: 2,
        };
        let calls = vec![
            MethCall {
                chrom_id: 0,
                pos: 0,
                m: 0,
                t: 1,
            },
            MethCall {
                chrom_id: 0,
                pos: 100,
                m: 1,
                t: 1,
            },
        ];
        let w = build_window(&f, &r, &calls, 0.0, 1);
        assert_eq!(pair_counts(&w, 1, 100).total(), 1);
        assert_eq!(pair_counts(&w, 1, 99).total(), 0);
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
        let pc = pair_counts(&w, 10, 0);
        assert_eq!(pc.total(), 0);
    }

    #[test]
    fn pair_counts_all_lags_matches_single_lag() {
        // Five observations spanning ref indices 0..5; cross-check the multi-lag sweep
        // against repeated single-lag calls.
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
        let all = pair_counts_all_lags(&w, 4, 0);
        assert_eq!(all.len(), 4);
        for k in 1..=4 {
            let single = pair_counts(&w, k, 0);
            assert_eq!(
                all[k - 1].counts,
                single.counts,
                "lag {} mismatch: all={:?} single={:?}",
                k,
                all[k - 1].counts,
                single.counts
            );
        }
    }

    /// Brute-force reference implementation against which the optimized sweep is checked.
    /// Mirrors the pre-rewrite dense-window algorithm: walk every reference slot, treat
    /// unobserved slots as gaps, and emit one PairCounts per lag.
    fn brute_force_pair_counts(
        feature_len: usize,
        slot_value: &[Option<u8>],
        slot_pos: &[u64],
        k_max: usize,
        max_distance: u64,
    ) -> Vec<PairCounts> {
        assert_eq!(slot_value.len(), feature_len);
        assert_eq!(slot_pos.len(), feature_len);
        let mut out = vec![PairCounts::default(); k_max];
        for lag in 1..=k_max {
            if lag >= feature_len {
                continue;
            }
            for i in 0..(feature_len - lag) {
                if max_distance > 0 && slot_pos[i + lag] - slot_pos[i] > max_distance {
                    continue;
                }
                if let (Some(a), Some(b)) = (slot_value[i], slot_value[i + lag]) {
                    out[lag - 1].counts[(a as usize) * 2 + b as usize] += 1;
                }
            }
        }
        out
    }

    /// xorshift32 — deterministic, no external rng dependency, good enough for tests.
    fn rng_next(state: &mut u32) -> u32 {
        let mut x = *state;
        x ^= x << 13;
        x ^= x >> 17;
        x ^= x << 5;
        *state = x;
        x
    }

    #[test]
    fn pair_counts_all_lags_matches_brute_force_random() {
        // Stress test: 200 reference CpGs, irregular spacing, ~30% observed, varying
        // values. Run several seeds and several max_distance settings.
        let feature_len = 200usize;
        let mut positions = Vec::with_capacity(feature_len);
        // Irregular but monotonic positions so max_distance has bite.
        let mut p = 100u64;
        let mut spacing_state = 0xdead_beef_u32;
        for _ in 0..feature_len {
            positions.push(p);
            let gap = (rng_next(&mut spacing_state) % 50) as u64 + 1;
            p += gap;
        }

        let chrom_id_of: std::collections::HashMap<String, u32> =
            [("chr1".into(), 0u32)].into_iter().collect();
        let reference = CpgReference {
            chrom_names: vec!["chr1".into()],
            chrom_id_of,
            positions: vec![positions.clone()],
        };
        let feature = Feature {
            feature_id: "f".into(),
            chrom_id: 0,
            start: 0,
            end: positions[feature_len - 1] + 10,
            cpg_start_idx: 0,
            cpg_end_idx: feature_len,
        };

        let k_max = 8usize;
        let max_distances = [0u64, 50, 200, 1000];

        for seed in [1u32, 42, 9999, 0x5a5a5a5a] {
            let mut state = seed;
            // Build synthetic dense slot arrays (Option<u8> per reference slot),
            // and the matching sorted MethCall list for build_window.
            let mut slot_value: Vec<Option<u8>> = vec![None; feature_len];
            let mut calls: Vec<MethCall> = Vec::new();
            for i in 0..feature_len {
                if rng_next(&mut state) % 100 < 30 {
                    let v = (rng_next(&mut state) & 1) as u8;
                    slot_value[i] = Some(v);
                    calls.push(MethCall {
                        chrom_id: 0,
                        pos: positions[i],
                        m: v as u32,
                        t: 1,
                    });
                }
            }
            let window = build_window(&feature, &reference, &calls, 0.0, 1);

            for &md in &max_distances {
                let got = pair_counts_all_lags(&window, k_max, md);
                let want = brute_force_pair_counts(feature_len, &slot_value, &positions, k_max, md);
                for k in 1..=k_max {
                    assert_eq!(
                        got[k - 1].counts,
                        want[k - 1].counts,
                        "seed={seed} max_distance={md} lag={k}: optimized={:?} brute={:?}",
                        got[k - 1].counts,
                        want[k - 1].counts,
                    );
                }
                // Also check the single-lag pair_counts API matches the multi-lag sweep.
                for k in 1..=k_max {
                    let single = pair_counts(&window, k, md);
                    assert_eq!(
                        single.counts,
                        got[k - 1].counts,
                        "single vs sweep mismatch at lag {k}",
                    );
                }
            }
        }
    }

    #[test]
    fn pair_counts_all_lags_respects_max_distance() {
        // Same fixture as max_distance_drops_far_pairs but via the multi-lag sweep.
        let r = CpgReference {
            chrom_names: vec!["chr1".into()],
            chrom_id_of: [("chr1".into(), 0u32)].into_iter().collect(),
            positions: vec![vec![10, 20, 1100]],
        };
        let f = Feature {
            feature_id: "f".into(),
            chrom_id: 0,
            start: 0,
            end: 2000,
            cpg_start_idx: 0,
            cpg_end_idx: 3,
        };
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
                pos: 1100,
                m: 1,
                t: 1,
            },
        ];
        let w = build_window(&f, &r, &calls, 0.0, 1);
        let all = pair_counts_all_lags(&w, 2, 1000);
        assert_eq!(all[0].counts[2], 1); // lag-1 pair (1,0) at distance 10
        assert_eq!(all[0].total(), 1);
        assert_eq!(all[1].total(), 0); // lag-2 pair distance 1090 dropped
    }
}
