use std::error::Error;
use vergen::EmitBuilder;

fn main() -> Result<(), Box<dyn Error>> {
    EmitBuilder::builder()
        .fail_on_error()
        .all_git()
        .git_describe(true, false, Some("ThisPatternShouldNotMatchAnythingEver"))
        .emit()?;

    println!("cargo:rerun-if-changed=Cargo.toml");
    println!("cargo:rerun-if-changed=src");

    Ok(())
}
