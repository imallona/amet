//! Across-cell Jensen-Shannon divergence on per-cell L-mer histograms.
//!
//! Each cell contributes a 4-bin histogram of 2-mer counts at lag 1 within a feature.
//! Cells with zero pair observations are excluded from the divergence.

use super::shannon_entropy;
use crate::kmer::PairCounts;

/// Streaming accumulator for multi-distribution JSD.
///
/// Folds one cell's lag-1 2-mer histogram at a time, so per-cell counts need
/// not be retained for the whole group. Holds a running sum of normalised
/// distributions and of per-cell entropies; both terms of the JSD are means
/// over cells, so they accumulate exactly.
#[derive(Debug, Clone, Default)]
pub struct JsdAccumulator {
    mixture_sum: [f64; 4],
    entropy_sum: f64,
    n: u64,
}

impl JsdAccumulator {
    /// Fold one cell's lag-1 2-mer counts. Cells with no pairs are ignored,
    /// matching the non-empty filter in `multi_jsd`.
    pub fn add(&mut self, cell: &PairCounts) {
        let total = cell.total();
        if total == 0 {
            return;
        }
        let t = total as f64;
        let mut p = [0.0f64; 4];
        for (k, slot) in p.iter_mut().enumerate() {
            *slot = cell.counts[k] as f64 / t;
            self.mixture_sum[k] += *slot;
        }
        self.entropy_sum += entropy_of_distribution(&p);
        self.n += 1;
    }

    /// JSD = H(mean P_i) - mean H(P_i). Returns 0 with fewer than 2 non-empty cells.
    pub fn finish(&self) -> f64 {
        if self.n < 2 {
            return 0.0;
        }
        let n = self.n as f64;
        let mixture = [
            self.mixture_sum[0] / n,
            self.mixture_sum[1] / n,
            self.mixture_sum[2] / n,
            self.mixture_sum[3] / n,
        ];
        let jsd = entropy_of_distribution(&mixture) - self.entropy_sum / n;
        if jsd < 0.0 { 0.0 } else { jsd }
    }
}

/// Multi-distribution generalised JSD (in bits, using log base 2).
///
/// JSD(P_1, ..., P_n) = H(mean P_i) - mean H(P_i).
///
/// Returns 0 when fewer than 2 cells have non-zero histograms.
pub fn multi_jsd(per_cell: &[PairCounts]) -> f64 {
    let mut acc = JsdAccumulator::default();
    for cell in per_cell {
        acc.add(cell);
    }
    acc.finish()
}

fn entropy_of_distribution(p: &[f64]) -> f64 {
    let mut h = 0.0;
    for &q in p {
        if q > 0.0 {
            h -= q * q.log2();
        }
    }
    h
}

/// Convenience: when only entropies are needed (e.g., for testing), compute Shannon
/// entropy on raw counts instead.
pub fn jsd_from_counts(per_cell: &[PairCounts]) -> f64 {
    multi_jsd(per_cell)
}

#[allow(dead_code)]
fn shannon_check(counts: &[u32]) -> f64 {
    shannon_entropy(counts)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn pc(c00: u32, c01: u32, c10: u32, c11: u32) -> PairCounts {
        PairCounts {
            counts: [c00, c01, c10, c11],
        }
    }

    #[test]
    fn identical_cells_give_zero_jsd() {
        let cells = vec![pc(10, 10, 10, 10), pc(20, 20, 20, 20), pc(5, 5, 5, 5)];
        let j = multi_jsd(&cells);
        assert!(j.abs() < 1e-12, "expected 0, got {}", j);
    }

    #[test]
    fn maximally_divergent_cells_give_max_jsd() {
        // Two cells, each entirely on a different bin.
        let cells = vec![pc(10, 0, 0, 0), pc(0, 10, 0, 0)];
        let j = multi_jsd(&cells);
        // H(mix) = H({0.5, 0.5, 0, 0}) = 1; H_avg = 0; JSD = 1 bit.
        assert!((j - 1.0).abs() < 1e-12);
    }

    #[test]
    fn single_cell_yields_zero() {
        let cells = vec![pc(10, 10, 10, 10)];
        assert_eq!(multi_jsd(&cells), 0.0);
    }

    #[test]
    fn empty_cells_filtered() {
        let cells = vec![pc(0, 0, 0, 0), pc(0, 0, 0, 0), pc(10, 0, 0, 0)];
        assert_eq!(multi_jsd(&cells), 0.0);
    }

    #[test]
    fn jsd_bounded_above_by_log_n_cells() {
        // For n cells, JSD <= log2(n).
        let cells = vec![
            pc(10, 0, 0, 0),
            pc(0, 10, 0, 0),
            pc(0, 0, 10, 0),
            pc(0, 0, 0, 10),
        ];
        let j = multi_jsd(&cells);
        assert!(j > 0.0);
        assert!(j <= 2.0 + 1e-12, "JSD {} exceeded log2(4) = 2", j);
    }

    #[test]
    fn jsd_increases_with_divergence() {
        let close = vec![pc(50, 25, 15, 10), pc(48, 27, 16, 9)];
        let far = vec![pc(50, 25, 15, 10), pc(10, 15, 25, 50)];
        assert!(multi_jsd(&far) > multi_jsd(&close));
    }

    #[test]
    fn accumulator_matches_multi_jsd() {
        // The streaming accumulator must yield bit-identical results to the
        // slice-based multi_jsd, since the workflow relies on the streaming path.
        let cases: Vec<Vec<PairCounts>> = vec![
            vec![pc(10, 10, 10, 10), pc(20, 20, 20, 20), pc(5, 5, 5, 5)],
            vec![pc(50, 25, 15, 10), pc(10, 15, 25, 50), pc(0, 0, 0, 0)],
            vec![pc(10, 0, 0, 0), pc(0, 10, 0, 0), pc(0, 0, 10, 0)],
            vec![pc(7, 3, 1, 9)],
            vec![],
        ];
        for case in &cases {
            let mut acc = JsdAccumulator::default();
            for c in case {
                acc.add(c);
            }
            assert_eq!(
                acc.finish(),
                multi_jsd(case),
                "accumulator diverged from multi_jsd on {:?}",
                case
            );
        }
    }
}
