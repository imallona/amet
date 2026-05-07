use thiserror::Error;

#[derive(Debug, Error)]
pub enum AmetError {
    #[error("io: {0}")]
    Io(#[from] std::io::Error),

    #[error("parse error in {file} line {line}: {msg}")]
    Parse {
        file: String,
        line: usize,
        msg: String,
    },

    #[error("manifest error: {0}")]
    Manifest(String),

    #[error("config error: {0}")]
    Config(String),
}

pub type Result<T> = std::result::Result<T, AmetError>;
