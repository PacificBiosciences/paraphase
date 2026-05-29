use thiserror::Error;

#[derive(Debug, Error)]
pub enum ParaphaseError {
    #[error("{0}")]
    Message(String),

    #[error(transparent)]
    Io(#[from] std::io::Error),

    #[error(transparent)]
    Utf8(#[from] std::str::Utf8Error),

    #[error(transparent)]
    FromUtf8(#[from] std::string::FromUtf8Error),

    #[error(transparent)]
    ParseInt(#[from] std::num::ParseIntError),

    #[error(transparent)]
    ParseFloat(#[from] std::num::ParseFloatError),

    #[error(transparent)]
    TryFromInt(#[from] std::num::TryFromIntError),

    #[error(transparent)]
    SerdeJson(#[from] serde_json::Error),

    #[error(transparent)]
    SerdeYaml(#[from] serde_yaml::Error),

    #[error(transparent)]
    Regex(#[from] regex::Error),

    #[error(transparent)]
    Nul(#[from] std::ffi::NulError),

    #[error(transparent)]
    RustHtslib(#[from] rust_htslib::errors::Error),

    #[error(transparent)]
    Url(#[from] url::ParseError),

    #[error(transparent)]
    Anyhow(#[from] anyhow::Error),

    #[error(transparent)]
    Simple(#[from] simple_error::SimpleError),

    #[error("{0}")]
    Boxed(String),
}

impl From<&str> for ParaphaseError {
    fn from(value: &str) -> Self {
        Self::Message(value.to_string())
    }
}

impl From<String> for ParaphaseError {
    fn from(value: String) -> Self {
        Self::Message(value)
    }
}

impl From<Box<dyn std::error::Error>> for ParaphaseError {
    fn from(value: Box<dyn std::error::Error>) -> Self {
        Self::Boxed(value.to_string())
    }
}

impl From<Box<dyn std::error::Error + Send + Sync>> for ParaphaseError {
    fn from(value: Box<dyn std::error::Error + Send + Sync>) -> Self {
        Self::Boxed(value.to_string())
    }
}

impl From<crate::assembly::variant_graph::GraphError> for ParaphaseError {
    fn from(value: crate::assembly::variant_graph::GraphError) -> Self {
        Self::Message(value.to_string())
    }
}
