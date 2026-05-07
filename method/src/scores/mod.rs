pub mod i_total;
pub mod jsd;

/// Plug-in Shannon entropy in bits, with optional Miller-Madow bias correction.
/// Counts that sum to 0 return 0.
pub fn shannon_entropy(counts: &[u32]) -> f64 {
    let n: u32 = counts.iter().sum();
    if n == 0 {
        return 0.0;
    }
    let n_f = n as f64;
    let mut h = 0.0;
    for &c in counts {
        if c > 0 {
            let p = c as f64 / n_f;
            h -= p * p.log2();
        }
    }
    h
}

/// Miller-Madow correction: H + (K_observed - 1) / (2N), in bits divided by ln(2).
/// Note: Miller-Madow is derived in nats; we convert to bits.
pub fn miller_madow_correction(counts: &[u32]) -> f64 {
    let n: u32 = counts.iter().sum();
    if n == 0 {
        return 0.0;
    }
    let k = counts.iter().filter(|&&c| c > 0).count() as f64;
    (k - 1.0) / (2.0 * n as f64) / std::f64::consts::LN_2
}

/// Plug-in Shannon entropy with Miller-Madow correction applied.
pub fn shannon_entropy_mm(counts: &[u32]) -> f64 {
    shannon_entropy(counts) + miller_madow_correction(counts)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn entropy_uniform_2() {
        // {1, 1} → log2(2) = 1.
        assert!((shannon_entropy(&[1, 1]) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn entropy_uniform_4() {
        // {1, 1, 1, 1} → log2(4) = 2.
        assert!((shannon_entropy(&[1, 1, 1, 1]) - 2.0).abs() < 1e-12);
    }

    #[test]
    fn entropy_degenerate_zero() {
        assert_eq!(shannon_entropy(&[5, 0, 0, 0]), 0.0);
    }

    #[test]
    fn entropy_empty_zero() {
        assert_eq!(shannon_entropy(&[0, 0]), 0.0);
    }

    #[test]
    fn miller_madow_zero_for_empty() {
        assert_eq!(miller_madow_correction(&[0, 0]), 0.0);
    }

    #[test]
    fn miller_madow_two_observed_bins_n10() {
        // K=2, N=10: correction = 1 / (2*10) / ln(2) ≈ 0.0721 bits.
        let c = miller_madow_correction(&[5, 5]);
        assert!((c - 0.07213475).abs() < 1e-6);
    }
}
