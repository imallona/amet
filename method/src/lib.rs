pub mod cli;
pub mod error;
pub mod features;
pub mod io;
pub mod kmer;
pub mod manifest;
pub mod parsers;
pub mod reference;
pub mod scores;

pub use error::{AmetError, Result};

// pos is the 0-based start of the CpG dinucleotide on the + strand.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MethCall {
    pub chrom_id: u32,
    pub pos: u64,
    pub m: u32,
    pub t: u32,
}

impl MethCall {
    pub fn binarize(&self, threshold: f64, min_reads: u32) -> Option<u8> {
        if self.t < min_reads {
            return None;
        }
        let frac = self.m as f64 / self.t as f64;
        Some(if frac > threshold { 1 } else { 0 })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn binarize_below_min_reads() {
        let c = MethCall {
            chrom_id: 0,
            pos: 100,
            m: 1,
            t: 1,
        };
        assert_eq!(c.binarize(0.0, 2), None);
    }

    #[test]
    fn binarize_default_threshold() {
        let c = MethCall {
            chrom_id: 0,
            pos: 100,
            m: 1,
            t: 5,
        };
        assert_eq!(c.binarize(0.0, 1), Some(1));
        let c = MethCall {
            chrom_id: 0,
            pos: 100,
            m: 0,
            t: 5,
        };
        assert_eq!(c.binarize(0.0, 1), Some(0));
    }

    #[test]
    fn binarize_majority_threshold() {
        let c = MethCall {
            chrom_id: 0,
            pos: 100,
            m: 2,
            t: 5,
        };
        assert_eq!(c.binarize(0.5, 1), Some(0));
        let c = MethCall {
            chrom_id: 0,
            pos: 100,
            m: 3,
            t: 5,
        };
        assert_eq!(c.binarize(0.5, 1), Some(1));
    }

    #[test]
    fn binarize_crc_threshold() {
        // CRC analysis used 0.1 fractional threshold.
        let c = MethCall {
            chrom_id: 0,
            pos: 100,
            m: 1,
            t: 9,
        };
        assert_eq!(c.binarize(0.1, 1), Some(1)); // 1/9 = 0.111 > 0.1
        let c = MethCall {
            chrom_id: 0,
            pos: 100,
            m: 1,
            t: 10,
        };
        assert_eq!(c.binarize(0.1, 1), Some(0)); // 1/10 = 0.1 not > 0.1
    }
}
