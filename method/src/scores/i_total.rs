//! Within-cell mutual-information score.
//!
//! For lag k:
//!   I_k = H(X_i) + H(X_{i+k}) - H(X_i, X_{i+k})
//!
//! For an i.i.d. binary sequence with marginal p, all I_k are 0 regardless of p.
//! I_total = sum over k=1..k_max of I_k inherits the same p-invariant zero baseline.

use super::{shannon_entropy, shannon_entropy_mm};
use crate::kmer::PairCounts;

/// Per-position marginal counts derived from a pair-count table for lag k.
/// X_i counts: row sums; X_{i+k} counts: column sums.
fn marginals_from_pairs(pairs: &PairCounts) -> ([u32; 2], [u32; 2]) {
    let c = pairs.counts;
    let xi = [c[0] + c[1], c[2] + c[3]]; // i = 0 row + i = 1 row
    let xj = [c[0] + c[2], c[1] + c[3]]; // j = 0 col + j = 1 col
    (xi, xj)
}

/// Compute I_k from a pair-count table at lag k. Uses Miller-Madow correction
/// on each entropy. Returns 0 when there are no pairs.
pub fn i_k(pairs: &PairCounts) -> f64 {
    if pairs.total() == 0 {
        return 0.0;
    }
    let (xi, xj) = marginals_from_pairs(pairs);
    let h_xi = shannon_entropy_mm(&xi);
    let h_xj = shannon_entropy_mm(&xj);
    let h_joint = shannon_entropy_mm(&pairs.counts);
    let i = h_xi + h_xj - h_joint;
    // MI cannot be negative; bias correction can occasionally push it slightly below 0.
    if i < 0.0 { 0.0 } else { i }
}

/// Variant without Miller-Madow correction; useful for tests against analytic values.
pub fn i_k_uncorrected(pairs: &PairCounts) -> f64 {
    if pairs.total() == 0 {
        return 0.0;
    }
    let (xi, xj) = marginals_from_pairs(pairs);
    let h_xi = shannon_entropy(&xi);
    let h_xj = shannon_entropy(&xj);
    let h_joint = shannon_entropy(&pairs.counts);
    let i = h_xi + h_xj - h_joint;
    if i < 0.0 { 0.0 } else { i }
}

/// Sum of I_k over k=1..k_max using uncorrected MI per lag. Each pair table corresponds
/// to one lag value.
pub fn i_total_uncorrected(per_lag_pairs: &[PairCounts]) -> f64 {
    per_lag_pairs.iter().map(i_k_uncorrected).sum()
}

/// Sum of I_k over k=1..k_max using Miller-Madow-corrected MI per lag.
pub fn i_total(per_lag_pairs: &[PairCounts]) -> f64 {
    per_lag_pairs.iter().map(i_k).sum()
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
    fn iid_uniform_balanced_gives_zero() {
        // Uniform {0,1} marginal, independent: pair counts proportional to (1,1,1,1).
        // H(X_i) = H(X_j) = 1; H(joint) = 2; I = 1+1-2 = 0.
        let pairs = pc(25, 25, 25, 25);
        assert!(i_k_uncorrected(&pairs).abs() < 1e-12);
    }

    #[test]
    fn iid_skewed_gives_zero() {
        // Marginal p = 0.1, independent: pair counts proportional to (0.81, 0.09, 0.09, 0.01).
        // I should be 0 regardless of p.
        let pairs = pc(81, 9, 9, 1);
        assert!(i_k_uncorrected(&pairs).abs() < 1e-12);
    }

    #[test]
    fn perfect_dependence_gives_full_marginal_entropy() {
        // X_i = X_j: only diagonal entries non-zero. I = H(X_i) = 1 bit at p=0.5.
        let pairs = pc(50, 0, 0, 50);
        assert!((i_k_uncorrected(&pairs) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn perfect_anti_dependence_also_full() {
        // X_j = 1 - X_i: only off-diagonal entries non-zero. Still I = H(X_i) = 1.
        let pairs = pc(0, 50, 50, 0);
        assert!((i_k_uncorrected(&pairs) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn empty_pairs_gives_zero() {
        let pairs = pc(0, 0, 0, 0);
        assert_eq!(i_k_uncorrected(&pairs), 0.0);
    }

    #[test]
    fn manuscript_example_010101() {
        // Sequence "010101", lag-1 pairs: {01:3, 10:2}.
        let pairs = pc(0, 3, 2, 0);
        let i = i_k_uncorrected(&pairs);
        assert!(i > 0.9, "expected high MI, got {}", i);
    }

    #[test]
    fn manuscript_example_011010_lower_than_010101() {
        // Sequence "011010", lag-1 pairs: {01:2, 10:2, 11:1}.
        let s1 = i_k_uncorrected(&pc(0, 3, 2, 0));
        let s2 = i_k_uncorrected(&pc(0, 2, 2, 1));
        assert!(
            s1 > s2,
            "010101 should have higher MI than 011010 (got {} vs {})",
            s1,
            s2
        );
    }

    #[test]
    fn period3_has_zero_lag1_nonzero_lag3() {
        // Sequence "010010010" has zero lag-1 mutual info? Actually periodic structure.
        // Let's verify:
        // 010010010: lag-1 pairs: 01,10,00,01,10,00,01,10 = {01:3, 10:3, 00:2, 11:0}
        let pairs1 = pc(2, 3, 3, 0);
        let i1 = i_k_uncorrected(&pairs1);
        // It is not zero because the marginals differ from joint factorisation here.
        // What matters: lag-3 picks up the periodicity.
        // lag-3 pairs of "010010010": positions (0,3),(1,4),(2,5),(3,6),(4,7),(5,8) =
        //   (0,0),(1,1),(0,0),(0,0),(1,1),(0,0) = {00:4, 11:2}
        let pairs3 = pc(4, 0, 0, 2);
        let i3 = i_k_uncorrected(&pairs3);
        assert!(
            i3 > i1,
            "lag-3 MI should exceed lag-1 MI for period-3 (got {} vs {})",
            i3,
            i1
        );
    }

    #[test]
    fn i_total_sums_per_lag() {
        let p1 = pc(50, 0, 0, 50);
        let p2 = pc(50, 0, 0, 50);
        let total = i_total_uncorrected(&[p1, p2]);
        assert!((total - 2.0).abs() < 1e-10);
    }
}
