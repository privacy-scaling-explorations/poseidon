use crate::{Spec, State};
use halo2curves::FieldExt;

/// Poseidon hasher that hashes constant input length and single output
#[derive(Debug, Clone)]
pub struct Poseidon<F: FieldExt, const LEN: usize, const T: usize, const RATE: usize> {
    spec: Spec<F, T, RATE>,
}

impl<F: FieldExt, const LEN: usize, const T: usize, const RATE: usize> Poseidon<F, LEN, T, RATE> {
    /// Constructs a clear state poseidon instance
    pub fn new(r_f: usize, r_p: usize) -> Self {
        Self {
            spec: Spec::new(r_f, r_p),
        }
    }

    /// Perform hashing
    pub fn hash(&self, elements: &[F; LEN]) -> F {
        let mut state = State::<F, T>::init_constant_length_mode();

        for chunk in elements.chunks(RATE) {
            // If chunk size is less than RATE pad with 0:
            // Zipping with different lengths will automatically do that.

            // Add new chunk of inputs for the next permutation cycle.
            for (input_element, state) in chunk.iter().zip(state.0.iter_mut().skip(1)) {
                state.add_assign(input_element);
            }

            self.spec.permute(&mut state);
        }
        state.result()
    }
}
