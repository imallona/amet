use crate::error::Result;
use flate2::Compression;
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

/// Open a file, transparently decompressing if its extension is .gz or .bgz.
pub fn open_read(path: &Path) -> Result<Box<dyn BufRead>> {
    let file = File::open(path)?;
    let ext = path.extension().and_then(|s| s.to_str()).unwrap_or("");
    if ext == "gz" || ext == "bgz" {
        Ok(Box::new(BufReader::new(MultiGzDecoder::new(file))))
    } else {
        Ok(Box::new(BufReader::new(file)))
    }
}

/// Open a file for writing, gzipping if the path ends with .gz.
pub fn open_write(path: &Path) -> Result<Box<dyn Write>> {
    let file = File::create(path)?;
    let ext = path.extension().and_then(|s| s.to_str()).unwrap_or("");
    if ext == "gz" {
        Ok(Box::new(BufWriter::new(GzEncoder::new(
            file,
            Compression::default(),
        ))))
    } else {
        Ok(Box::new(BufWriter::new(file)))
    }
}
