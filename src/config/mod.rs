pub mod gene;
pub mod region;
pub(crate) mod schema;

pub use gene::Config as Gene;
pub use region::Config as Region;
pub use region::Locus;

pub const REGION_CONFIG_HG19: &[u8] =
    std::include_bytes!(concat!(env!("CARGO_MANIFEST_DIR"), "/data/19/config.yaml"));
pub const REGION_CONFIG_HG38: &[u8] =
    std::include_bytes!(concat!(env!("CARGO_MANIFEST_DIR"), "/data/38/config.yaml"));
pub const REGION_CONFIG_CHM13: &[u8] = std::include_bytes!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/data/chm13/config.yaml"
));
lazy_static::lazy_static! {
    pub static ref CONFIG: Region = region::Config::load(None); // Default config: hg38
    pub static ref CONFIG_HG19: Region = region::Config::load(Some(REGION_CONFIG_HG19));
    pub static ref CONFIG_HG38: Region = region::Config::load(Some(REGION_CONFIG_HG38));
}
