pub type Error = simple_error::SimpleError;
type JsonMap = serde_json::Map<String, serde_json::Value>;

mod decode;
mod errors;
mod model;
mod output;
mod source;
mod transform;
mod types;
mod validate;

pub use model::{ParaphaseOutput, ParsedParaphaseOutputJSON};
pub use output::{write_outputs, GeneCall};
pub use types::{PhasingSite, ReadAlignmentId, ReadFingerprintMap};
