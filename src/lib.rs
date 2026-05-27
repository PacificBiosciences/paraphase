#![allow(clippy::cast_sign_loss)]
#![allow(clippy::cast_possible_truncation)]
#![allow(clippy::cast_possible_wrap)]

// Core functionality
pub mod assembly;
pub mod phaser;
pub mod toolkit;

// Gene-specific callers
pub mod genes;

// BAM utilities
pub mod depth;
pub mod realign;

// Configuration
pub mod config;

// Command-line interface
pub mod cli;

// Error types
pub mod error;

// Input/output utilities
pub mod io;
