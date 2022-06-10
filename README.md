# `poseidon`

`poseidon` is built to be used in SNARK and non native recursion friendly transcript for [appliedzkp/halo2](https://github.com/appliedzkp/halo2/).

[Poseidon hash function](https://eprint.iacr.org/2019/458.pdf) implmenetation is in line with the reference and the [test vectors](https://extgit.iaik.tugraz.at/krypto/hadeshash/-/tree/master/code). It also uses optimized constants and sparse MDS matrices to reduce number of multiplications. For now constants are calculated in construction time they are planned to be hardcoded once transcript design matures. Currently only supports variable length hashing with $\alpha = 5$ sbox. Some parts of Poseidon implementation are adapted or ported from:

* [filecoin-project/neptune](https://github.com/filecoin-project/neptune/tree/master/spec)
* [matter-labs/rescue-poseidon](https://github.com/matter-labs/rescue-poseidon)
* [zcash/halo2/halo2_gadgets](https://github.com/zcash/halo2/tree/main/halo2_gadgets)
* [dusk-network/Poseidon252](https://github.com/dusk-network/Poseidon252)

## Example usage

```rust
// Initialize a mutable hasher with constant capacity parameters 
// and number of rounds arguments. This will also generate matrices 
// and constants according to the specification
let mut hasher = Poseidon::<Fr, T, RATE>::new(number_of_full_rounds, number_of_half_rounds);

// In sake of the example we generate some dummy scalar inputs
let inputs = (0..number_of_inputs_0)
    .map(|_| Fr::random(&mut rng))
    .collect::<Vec<Fr>>();

// Feed inputs to the Absorption line
hasher.update(&inputs[..]);

// Yield your challange with squeeze function
let challenge_alpha = hasher.squeeze();

// Then again ...
let inputs = (0..number_of_inputs_1)
    .map(|_| Fr::random(&mut rng))
    .collect::<Vec<Fr>>();
hasher.update(&inputs[..]);
let challenge_beta = hasher.squeeze();

```
