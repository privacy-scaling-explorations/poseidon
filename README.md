# transcript

Work in progress SNARK and non native recursion friendly transcript for [appliedzkp/halo2](https://github.com/appliedzkp/halo2/).

[Poseidon hash function](https://eprint.iacr.org/2019/458.pdf) implmenetation is in line with the reference and the [test vectors](https://extgit.iaik.tugraz.at/krypto/hadeshash/-/tree/master/code). It also uses optimized constants and sparse MDS matrices to reduce number of multiplications. For now constants are calculated in construction time they are planned to be hardcoded once transcript design matures. Some parts of Poseidon implementation are adapted or ported from:

* [filecoin-project/neptune](https://github.com/filecoin-project/neptune/tree/master/spec)
* [matter-labs/rescue-poseidon](https://github.com/matter-labs/rescue-poseidon)
* [zcash/halo2/halo2_gadgets](https://github.com/zcash/halo2/tree/main/halo2_gadgets)
* [dusk-network/Poseidon252](https://github.com/dusk-network/Poseidon252)

Map:

- [x] Optimized permutation
- [x] Paddind and absorbing
- [ ] Transcript pipeline with non native inputs
- [ ] Transcript circuit will be under [appliedzkp/halo2wrong](https://github.com/appliedzkp/halo2wrong)
