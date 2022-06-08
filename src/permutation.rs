use crate::spec::{Spec, State};
use halo2curves::FieldExt;

impl<F: FieldExt, const T: usize, const RATE: usize> Spec<F, T, RATE> {
    /// Applies the Poseidon permutation to the given state
    pub fn permute(&self, state: &mut State<F, T>) {
        let r_f = self.r_f / 2;

        // First half of the full rounds
        {
            state.add_constants(&self.constants.start[0]);
            for round_constants in self.constants.start.iter().skip(1).take(r_f - 1) {
                state.sbox_full();
                state.add_constants(round_constants);
                self.mds_matrices.mds.apply(state);
            }
            state.sbox_full();
            state.add_constants(self.constants.start.last().unwrap());
            self.mds_matrices.pre_sparse_mds.apply(state)
        }

        // Partial rounds
        {
            for (round_constant, sparse_mds) in self
                .constants
                .partial
                .iter()
                .zip(self.mds_matrices.sparse_matrices.iter())
            {
                state.sbox_part();
                state.add_constant(round_constant);
                sparse_mds.apply(state);
            }
        }

        // Second half of the full rounds
        {
            for round_constants in self.constants.end.iter() {
                state.sbox_full();
                state.add_constants(round_constants);
                self.mds_matrices.mds.apply(state);
            }
            state.sbox_full();
            self.mds_matrices.mds.apply(state);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::State;
    use crate::spec::{tests::SpecRef, Spec};
    use group::ff::PrimeField;
    use halo2curves::bn256::Fr;
    use halo2curves::FieldExt;

    /// We want to keep unoptimized poseidion construction and permutation to
    /// cross test with optimized one
    impl<F: FieldExt, const T: usize, const RATE: usize> SpecRef<F, T, RATE> {
        fn permute(&self, state: &mut State<F, T>) {
            let (r_f, r_p) = (self.r_f / 2, self.r_p);

            for constants in self.constants.iter().take(r_f) {
                state.add_constants(constants);
                state.sbox_full();
                self.mds.apply(state);
            }

            for constants in self.constants.iter().skip(r_f).take(r_p) {
                state.add_constants(constants);
                state.sbox_part();
                self.mds.apply(state);
            }

            for constants in self.constants.iter().skip(r_f + r_p) {
                state.add_constants(constants);
                state.sbox_full();
                self.mds.apply(state);
            }
        }
    }

    #[test]
    fn cross_test() {
        use halo2curves::group::ff::Field;
        use rand_core::OsRng;
        use std::time::Instant;

        macro_rules! run_test {
            (
                $([$RF:expr, $RP:expr, $T:expr, $RATE:expr]),*
            ) => {
                $(
                    {
                        const R_F: usize = $RF;
                        const R_P: usize = $RP;
                        const T: usize = $T;
                        const RATE: usize = $RATE;
                        let mut state = State(
                            (0..T)
                                .map(|_| Fr::random(OsRng))
                                .collect::<Vec<Fr>>()
                                .try_into().unwrap(),
                        );
                        let spec = SpecRef::<Fr, T, RATE>::new(R_F, R_P);
                        let mut state_expected = state.clone();
                        spec.permute(&mut state_expected);

                        let spec = Spec::<Fr, T, RATE>::new(R_F, R_P);
                        let now = Instant::now();
                        {
                            spec.permute(&mut state);
                        }
                        let elapsed = now.elapsed();
                        println!("Elapsed: {:.2?}", elapsed);
                        assert_eq!(state_expected, state);
                    }
                )*
            };
        }
        run_test!([8, 57, 3, 2]);
        run_test!([8, 57, 4, 3]);
        run_test!([8, 57, 5, 4]);
        run_test!([8, 57, 6, 5]);
        run_test!([8, 57, 7, 6]);
        run_test!([8, 57, 8, 7]);
        run_test!([8, 57, 9, 8]);
        run_test!([8, 57, 10, 9]);
    }

    #[test]
    fn test_against_test_vectors() {
        // https://extgit.iaik.tugraz.at/krypto/hadeshash/-/blob/master/code/test_vectors.txt
        // poseidonperm_x5_254_3
        {
            const R_F: usize = 8;
            const R_P: usize = 57;
            const T: usize = 3;
            const RATE: usize = 2;

            let state = State(
                vec![0u64, 1, 2]
                    .into_iter()
                    .map(Fr::from)
                    .collect::<Vec<Fr>>()
                    .try_into()
                    .unwrap(),
            );

            let spec_ref = SpecRef::<Fr, T, RATE>::new(R_F, R_P);
            let mut state_0 = state.clone();

            spec_ref.permute(&mut state_0);
            let expected = vec![
                "7853200120776062878684798364095072458815029376092732009249414926327459813530",
                "7142104613055408817911962100316808866448378443474503659992478482890339429929",
                "6549537674122432311777789598043107870002137484850126429160507761192163713804",
            ];
            for (word, expected) in state_0.words().into_iter().zip(expected.iter()) {
                assert_eq!(word, Fr::from_str_vartime(expected).unwrap());
            }

            let spec = Spec::<Fr, T, RATE>::new(R_F, R_P);
            let mut state_1 = state;
            spec.permute(&mut state_1);
            assert_eq!(state_0, state_1);
        }

        // https://extgit.iaik.tugraz.at/krypto/hadeshash/-/blob/master/code/test_vectors.txt
        // poseidonperm_x5_254_5
        {
            const R_F: usize = 8;
            const R_P: usize = 60;
            const T: usize = 5;
            const RATE: usize = 4;

            let state = State(
                vec![0u64, 1, 2, 3, 4]
                    .into_iter()
                    .map(Fr::from)
                    .collect::<Vec<Fr>>()
                    .try_into()
                    .unwrap(),
            );

            let spec_ref = SpecRef::<Fr, T, RATE>::new(R_F, R_P);
            let mut state_0 = state.clone();

            spec_ref.permute(&mut state_0);
            let expected = vec![
                "18821383157269793795438455681495246036402687001665670618754263018637548127333",
                "7817711165059374331357136443537800893307845083525445872661165200086166013245",
                "16733335996448830230979566039396561240864200624113062088822991822580465420551",
                "6644334865470350789317807668685953492649391266180911382577082600917830417726",
                "3372108894677221197912083238087960099443657816445944159266857514496320565191",
            ];
            for (word, expected) in state_0.words().into_iter().zip(expected.iter()) {
                assert_eq!(word, Fr::from_str_vartime(expected).unwrap());
            }

            let spec = Spec::<Fr, T, RATE>::new(R_F, R_P);
            let mut state_1 = state;
            spec.permute(&mut state_1);
            assert_eq!(state_0, state_1);
        }
    }
}
