//! Poseidon hashing implemention with variable lenght input setting. This crate
//! also exposes constant parameters for circuit implementations

#![deny(missing_debug_implementations)]
#![deny(missing_docs)]

mod grain;
mod matrix;
mod permutation;
mod poseidon;
mod spec;

pub use crate::poseidon::Poseidon;
pub use crate::spec::{MDSMatrices, MDSMatrix, SparseMDSMatrix, Spec, State};
