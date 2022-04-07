use crate::{Spec, State};
use pairing::arithmetic::FieldExt;

/// Poseidon hasher that maintains state and inputs and yields single element
/// output when desired
#[derive(Debug, Clone)]
pub struct Poseidon<F: FieldExt, const T: usize, const RATE: usize> {
    state: State<F, T>,
    spec: Spec<F, T, RATE>,
    absorbing: Vec<F>,
}

impl<F: FieldExt, const T: usize, const RATE: usize> Poseidon<F, T, RATE> {
    /// Constructs a clear state poseidon instance
    pub fn new(r_f: usize, r_p: usize) -> Self {
        Self {
            spec: Spec::new(r_f, r_p),
            state: State::default(),
            absorbing: Vec::new(),
        }
    }

    /// Appends elements to the absorbation line
    pub fn update(&mut self, elements: &[F]) {
        // TODO: if `RATE` is full spawn a permutation here
        self.absorbing.extend_from_slice(elements);
    }

    /// Results a single element by absorbing already added inputs
    pub fn squeeze(&mut self) -> F {
        let mut padding_offset = 0;
        let input_elements = self.absorbing.clone();
        for chunk in input_elements.chunks(RATE) {
            // `padding_offset` is expected to be in [0, RATE-1]
            padding_offset = RATE - chunk.len();
            let mut chunk = chunk.to_vec();
            if padding_offset > 0 {
                chunk.push(F::one())
            }

            // Add new chunk of inputs for the next permutation cycle
            for (input_element, state) in chunk.iter().zip(self.state.0.iter_mut().skip(1)) {
                state.add_assign(input_element);
            }
            self.permute();
        }

        // TODO/FIX: can we avoid this in the on going transcript squeezing context?
        self.finalize_padding(padding_offset == 0);
        self.finalize()
    }

    fn permute(&mut self) {
        self.spec.permute(&mut self.state)
    }

    // If there is no room left for the padding sign add it to the state and
    // perform one more permutation
    fn finalize_padding(&mut self, must_perform: bool) {
        if must_perform {
            self.state.0[1].add_assign(F::one());
            self.permute();
        }
    }

    /// Flushes absorbation line and returns the result
    fn finalize(&mut self) -> F {
        self.absorbing.clear();
        self.state.result()
    }
}

#[test]
fn test_padding() {
    use group::ff::Field;
    use pairing::bn256::Fr;

    const R_F: usize = 8;
    const R_P: usize = 57;
    const T: usize = 5;
    const RATE: usize = 4;
    use rand::thread_rng;
    let mut rng = thread_rng();

    // w/o extra permutation
    {
        let mut poseidon = Poseidon::<Fr, T, RATE>::new(R_F, R_P);
        let number_of_permutation = 5;
        let number_of_inputs = RATE * number_of_permutation - 1;
        let inputs = (0..number_of_inputs)
            .map(|_| Fr::random(&mut rng))
            .collect::<Vec<Fr>>();
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

    // w/ extra permutation
    {
        let mut poseidon = Poseidon::<Fr, T, RATE>::new(R_F, R_P);
        let number_of_permutation = 5;
        let number_of_inputs = RATE * number_of_permutation;
        let inputs = (0..number_of_inputs)
            .map(|_| Fr::random(&mut rng))
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
}
