mod grain;
mod matrix;
mod permutation;
mod poseidon;
mod spec;
pub mod transcript;

use halo2::arithmetic::{BaseExt, FieldExt};
use num_bigint::BigUint as big_uint;
use num_traits::Num;
use std::ops::Shl;

pub(crate) use spec::{Spec, State};

pub(crate) fn decompose<Base: BaseExt, Scalar: FieldExt>(
    e: Base,
    number_of_limbs: usize,
    bit_len: usize,
) -> Vec<Scalar> {
    decompose_big(fe_to_big(e), number_of_limbs, bit_len)
}

fn native<Base: BaseExt, Scalar: FieldExt>(e: Base) -> Scalar {
    big_to_fe(fe_to_big(e) % modulus::<Scalar>())
}

fn modulus<F: FieldExt>() -> big_uint {
    big_uint::from_str_radix(&F::MODULUS[2..], 16).unwrap()
}

fn decompose_big<F: FieldExt>(e: big_uint, number_of_limbs: usize, bit_len: usize) -> Vec<F> {
    let mut e = e;
    let mask = big_uint::from(1usize).shl(bit_len) - 1usize;
    let limbs: Vec<F> = (0..number_of_limbs)
        .map(|_| {
            let limb = mask.clone() & e.clone();
            e = e.clone() >> bit_len;
            big_to_fe(limb)
        })
        .collect();

    limbs
}

fn big_to_fe<F: FieldExt>(e: big_uint) -> F {
    let modulus = modulus::<F>();
    let e = e % modulus;
    let mut bytes = e.to_bytes_le();
    bytes.resize(32, 0);
    let mut bytes = &bytes[..];
    F::read(&mut bytes).unwrap()
}

fn fe_to_big<F: BaseExt>(fe: F) -> big_uint {
    let mut bytes: Vec<u8> = Vec::new();
    fe.write(&mut bytes).unwrap();
    big_uint::from_bytes_le(&bytes[..])
}

fn sign<F: BaseExt>(fe: F) -> bool {
    let mut bytes: Vec<u8> = Vec::new();
    fe.write(&mut bytes).unwrap();
    (bytes[0] & 1) == 0
}
