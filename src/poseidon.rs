use crate::{Spec, State};
use halo2curves::group::ff::{FromUniformBytes, PrimeField};

/// Poseidon hasher that maintains state and inputs and yields single element
/// output when desired
#[derive(Debug, Clone)]
pub struct Poseidon<F: PrimeField, const T: usize, const RATE: usize> {
    state: State<F, T>,
    spec: Spec<F, T, RATE>,
    absorbing: Vec<F>,
}

impl<F: FromUniformBytes<64>, const T: usize, const RATE: usize> Poseidon<F, T, RATE> {
    /// Constructs a clear state poseidon instance
    pub fn new(r_f: usize, r_p: usize) -> Self {
        Self {
            spec: Spec::new(r_f, r_p),
            state: State::default(),
            absorbing: Vec::new(),
        }
    }

    /// Appends elements to the absorption line updates state while `RATE` is
    /// full
    pub fn update(&mut self, elements: &[F]) {
        let mut input_elements = self.absorbing.clone();
        input_elements.extend_from_slice(elements);

        for chunk in input_elements.chunks(RATE) {
            if chunk.len() < RATE {
                // Must be the last iteration of this update. Feed unpermutaed inputs to the
                // absorbation line
                self.absorbing = chunk.to_vec();
            } else {
                // Add new chunk of inputs for the next permutation cycle.
                for (input_element, state) in chunk.iter().zip(self.state.0.iter_mut().skip(1)) {
                    state.add_assign(input_element);
                }
                // Perform intermediate permutation
                self.spec.permute(&mut self.state);
                // Flush the absorption line
                self.absorbing.clear();
            }
        }
    }

    /// Results a single element by absorbing already added inputs
    pub fn squeeze(&mut self) -> F {
        let mut last_chunk = self.absorbing.clone();
        {
            // Expect padding offset to be in [0, RATE)
            debug_assert!(last_chunk.len() < RATE);
        }
        // Add the finishing sign of the variable length hashing. Note that this mut
        // also apply when absorbing line is empty
        last_chunk.push(F::ONE);
        // Add the last chunk of inputs to the state for the final permutation cycle

        for (input_element, state) in last_chunk.iter().zip(self.state.0.iter_mut().skip(1)) {
            state.add_assign(input_element);
        }

        // Perform final permutation
        self.spec.permute(&mut self.state);
        // Flush the absorption line
        self.absorbing.clear();
        // Returns the challenge while preserving internal state
        self.state.result()
    }
}

#[cfg(test)]
mod tests {
    use crate::{Poseidon, State};
    use halo2curves::bn256::Fr;
    use halo2curves::group::ff::Field;
    use paste::paste;
    use rand_core::OsRng;

    const R_F: usize = 8;
    const R_P: usize = 57;
    const T: usize = 5;
    const RATE: usize = 4;

    fn gen_random_vec(len: usize) -> Vec<Fr> {
        (0..len).map(|_| Fr::random(OsRng)).collect::<Vec<Fr>>()
    }

    #[test]
    fn poseidon_padding_with_last_chunk_len_is_not_rate_multiples() {
        let mut poseidon = Poseidon::<Fr, T, RATE>::new(R_F, R_P);
        let number_of_permutation = 5;
        let number_of_inputs = RATE * number_of_permutation - 1;
        let inputs = gen_random_vec(number_of_inputs);

        poseidon.update(&inputs[..]);
        let result_0 = poseidon.squeeze();

        let spec = poseidon.spec.clone();
        let mut inputs = inputs.clone();
        inputs.push(Fr::one());
        assert!(inputs.len() % RATE == 0);
        let mut state = State::<Fr, T>::default();
        for chunk in inputs.chunks(RATE) {
            let mut inputs = vec![Fr::zero()];
            inputs.extend_from_slice(chunk);
            state.add_constants(&inputs.try_into().unwrap());
            spec.permute(&mut state)
        }
        let result_1 = state.result();

        assert_eq!(result_0, result_1);
    }

    #[test]
    fn poseidon_padding_with_last_chunk_len_is_rate_multiples() {
        let mut poseidon = Poseidon::<Fr, T, RATE>::new(R_F, R_P);
        let number_of_permutation = 5;
        let number_of_inputs = RATE * number_of_permutation;
        let inputs = (0..number_of_inputs)
            .map(|_| Fr::random(OsRng))
            .collect::<Vec<Fr>>();
        poseidon.update(&inputs[..]);
        let result_0 = poseidon.squeeze();

        let spec = poseidon.spec.clone();
        let mut inputs = inputs.clone();
        let mut extra_padding = vec![Fr::zero(); RATE];
        extra_padding[0] = Fr::one();
        inputs.extend(extra_padding);

        assert!(inputs.len() % RATE == 0);
        let mut state = State::<Fr, T>::default();
        for chunk in inputs.chunks(RATE) {
            let mut inputs = vec![Fr::zero()];
            inputs.extend_from_slice(chunk);
            state.add_constants(&inputs.try_into().unwrap());
            spec.permute(&mut state)
        }
        let result_1 = state.result();

        assert_eq!(result_0, result_1);
    }

    macro_rules! test_padding {
        ($T:expr, $RATE:expr) => {
            paste! {
                #[test]
                fn [<test_padding_ $T _ $RATE>]() {
                    for number_of_iters in 1..25 {
                        let mut poseidon = Poseidon::<Fr, $T, $RATE>::new(R_F, R_P);

                        let mut inputs = vec![];
                        for number_of_inputs in 0..=number_of_iters {
                            let chunk = (0..number_of_inputs)
                                .map(|_| Fr::random(OsRng))
                                .collect::<Vec<Fr>>();
                            poseidon.update(&chunk[..]);
                            inputs.extend(chunk);
                        }
                        let result_0 = poseidon.squeeze();

                        // Accept below as reference and check consistency
                        inputs.push(Fr::one());
                        let offset = inputs.len() % $RATE;
                        if offset != 0 {
                            inputs.extend(vec![Fr::zero(); $RATE - offset]);
                        }

                        let spec = poseidon.spec.clone();
                        let mut state = State::<Fr, $T>::default();
                        for chunk in inputs.chunks($RATE) {
                            // First element is zero
                            let mut round_inputs = vec![Fr::zero()];
                            // Round inputs must be T sized now
                            round_inputs.extend_from_slice(chunk);

                            state.add_constants(&round_inputs.try_into().unwrap());
                            spec.permute(&mut state)
                        }
                        let result_1 = state.result();
                        assert_eq!(result_0, result_1);
                    }
                }
            }
        };
    }

    test_padding!(3, 2);
    test_padding!(4, 3);
    test_padding!(5, 4);
    test_padding!(6, 5);
    test_padding!(7, 6);
    test_padding!(8, 7);
    test_padding!(9, 8);
    test_padding!(10, 9);
}
