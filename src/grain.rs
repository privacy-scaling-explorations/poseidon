use crate::spec::MDSMatrix;
use halo2curves::group::ff::{FromUniformBytes, PrimeField};
use std::marker::PhantomData;

/// Grain initializes round constants and MDS matrix at given sponge parameters
pub(super) struct Grain<F: PrimeField, const T: usize, const RATE: usize> {
    bit_sequence: Vec<bool>,
    _field: PhantomData<F>,
}

impl<F: FromUniformBytes<64>, const T: usize, const RATE: usize> Grain<F, T, RATE> {
    pub(crate) fn generate(r_f: usize, r_p: usize) -> (Vec<[F; T]>, MDSMatrix<F, T, RATE>) {
        debug_assert!(T > 1 && T == RATE + 1);

        // Support only prime field construction
        const FIELD_TYPE: u8 = 1u8;
        // Support only \alpha s-box
        const SBOX_TYPE: u8 = 0;

        let field_size = F::NUM_BITS;
        let n_bytes = F::Repr::default().as_ref().len();
        assert_eq!((field_size as f32 / 8.0).ceil() as usize, n_bytes);
        assert_eq!(r_f % 2, 0);

        // Pseudo random number generation. See:
        // Initialization of the Grain LFSR Used for Parameter Generation
        // Supplementary Material Section F
        // https://eprint.iacr.org/2019/458.pdf
        let mut bit_sequence: Vec<bool> = Vec::new();
        append_bits(&mut bit_sequence, 2, FIELD_TYPE);
        append_bits(&mut bit_sequence, 4, SBOX_TYPE);
        append_bits(&mut bit_sequence, 12, field_size);
        append_bits(&mut bit_sequence, 12, T as u32);
        append_bits(&mut bit_sequence, 10, r_f as u16);
        append_bits(&mut bit_sequence, 10, r_p as u16);
        append_bits(&mut bit_sequence, 30, 0b111111111111111111111111111111u128);
        debug_assert_eq!(bit_sequence.len(), 80);

        let mut grain: Grain<F, T, RATE> = Grain {
            bit_sequence,
            _field: PhantomData,
        };

        for _ in 0..160 {
            grain.new_bit();
        }
        assert_eq!(grain.bit_sequence.len(), 80);

        let number_of_rounds = r_p + r_f;
        let constants = (0..number_of_rounds)
            .map(|_| {
                let mut round_constants = [F::ZERO; T];
                for c in round_constants.iter_mut() {
                    *c = grain.next_field_element();
                }
                round_constants
            })
            .collect::<Vec<[F; T]>>();

        let (mut xs, mut ys) = ([F::ZERO; T], [F::ZERO; T]);
        for x in xs.iter_mut() {
            *x = grain.next_field_element_without_rejection();
        }
        for y in ys.iter_mut() {
            *y = grain.next_field_element_without_rejection();
        }

        (constants, MDSMatrix::cauchy(&xs, &ys))
    }

    /// Credit: https://github.com/zcash/halo2/tree/main/halo2_gadgets/src/primitives/poseidon
    /// Returns the next field element from this Grain instantiation.
    pub(super) fn next_field_element(&mut self) -> F {
        // Loop until we get an element in the field.
        loop {
            let mut bytes = F::Repr::default();

            // Poseidon reference impl interprets the bits as a repr in MSB order, because
            // it's easy to do that in Python. Meanwhile, our field elements all use LSB
            // order. There's little motivation to diverge from the reference impl; these
            // are all constants, so we aren't introducing big-endianness into the rest of
            // the circuit (assuming unkeyed Poseidon, but we probably wouldn't want to
            // implement Grain inside a circuit, so we'd use a different round constant
            // derivation function there).
            let view = bytes.as_mut();
            for (i, bit) in self.take(F::NUM_BITS as usize).enumerate() {
                // If we diverged from the reference impl and interpreted the bits in LSB
                // order, we would remove this line.
                let i = F::NUM_BITS as usize - 1 - i;

                view[i / 8] |= if bit { 1 << (i % 8) } else { 0 };
            }

            if let Some(f) = F::from_repr_vartime(bytes) {
                break f;
            }
        }
    }

    /// Credit: https://github.com/zcash/halo2/tree/main/halo2_gadgets/src/primitives/poseidon
    /// Returns the next field element from this Grain instantiation, without
    /// using rejection sampling.
    pub(super) fn next_field_element_without_rejection(&mut self) -> F {
        let mut bytes = [0u8; 64];

        // Poseidon reference impl interprets the bits as a repr in MSB order, because
        // it's easy to do that in Python. Additionally, it does not use rejection
        // sampling in cases where the constants don't specifically need to be uniformly
        // random for security. We do not provide APIs that take a field-element-sized
        // array and reduce it modulo the field order, because those are unsafe APIs to
        // offer generally (accidentally using them can lead to divergence in consensus
        // systems due to not rejecting canonical forms).
        //
        // Given that we don't want to diverge from the reference implementation, we
        // hack around this restriction by serializing the bits into a 64-byte
        // array and then calling F::from_bytes_wide. PLEASE DO NOT COPY THIS
        // INTO YOUR OWN CODE!
        let view = bytes.as_mut();
        for (i, bit) in self.take(F::NUM_BITS as usize).enumerate() {
            // If we diverged from the reference impl and interpreted the bits in LSB
            // order, we would remove this line.
            let i = F::NUM_BITS as usize - 1 - i;

            view[i / 8] |= if bit { 1 << (i % 8) } else { 0 };
        }

        F::from_uniform_bytes(&bytes)
    }

    fn new_bit(&mut self) -> bool {
        // See supplementary material Section F. Step 2.
        // https://eprint.iacr.org/2019/458.pdf
        let new_bit = vec![62, 51, 38, 23, 13usize]
            .iter()
            .fold(self.bit_sequence[0], |acc, pos| {
                acc ^ self.bit_sequence[*pos]
            });
        assert_eq!(self.bit_sequence.len(), 80);
        self.bit_sequence.remove(0);
        self.bit_sequence.push(new_bit);
        new_bit
    }
}

impl<F: FromUniformBytes<64>, const T: usize, const RATE: usize> Iterator for Grain<F, T, RATE> {
    type Item = bool;

    fn next(&mut self) -> Option<Self::Item> {
        while !self.new_bit() {
            self.new_bit();
        }
        Some(self.new_bit())
    }
}

fn append_bits<T: Into<u128>>(vec: &mut Vec<bool>, n: usize, from: T) {
    let val = from.into();
    for i in (0..n).rev() {
        vec.push((val >> i) & 1 != 0);
    }
}
