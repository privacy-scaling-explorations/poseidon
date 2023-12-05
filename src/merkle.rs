use crate::{spec::MerkleMod, Spec, State};
use halo2curves::group::ff::{FromUniformBytes, PrimeField};

/// Stateless Merkle tree hasher where `RATE` is arity of the tree
#[derive(Debug, Clone)]
pub struct Merkle<F: PrimeField, const T: usize, const RATE: usize> {
    spec: Spec<F, T, RATE>,
}

impl<F: FromUniformBytes<64>, const T: usize, const RATE: usize> Merkle<F, T, RATE> {
    /// Constructs a clear state poseidon instance
    pub fn new(r_f: usize, r_p: usize) -> Self {
        Self {
            spec: Spec::new(r_f, r_p),
        }
    }

    /// Hashes the leafs to a node
    pub fn hash(&self, elements: &[F; RATE]) -> F {
        let mut state = State::new::<MerkleMod>();
        state
            .0
            .iter_mut()
            .skip(1)
            .zip(elements.iter())
            .for_each(|(s, e)| *s = *e);
        self.spec.permute(&mut state);
        state.result()
    }
}

#[cfg(test)]
mod tests {

    use crate::spec::MerkleMod;
    use crate::{MDSMatrix, Spec};
    use halo2curves::ff::{Field, PrimeField};
    use halo2curves::pasta::pallas::Scalar;

    fn hex<F: PrimeField<Repr = [u8; 32]>>(hex: &str) -> F {
        let repr = (0..32)
            .rev()
            .map(|i| u8::from_str_radix(&hex[i * 2..i * 2 + 2], 16).unwrap())
            .collect::<Vec<_>>();
        F::from_repr(repr.try_into().unwrap()).unwrap()
    }

    #[test]
    fn test_merkle() {
        const R_F: usize = 8;
        const R_P: usize = 55;
        const T: usize = 3;

        let mut state = crate::State::<Scalar, T>::new::<MerkleMod>();
        state.0[1] = Scalar::ONE;
        let spec = neptune_pallas_55_8_3_2(R_F, R_P);
        spec.permute(&mut state);

        // generated with:
        // https://github.com/lurk-lab/neptune/blob/28bed5841e4cec30d9e64e2fba529e0dba493457/src/poseidon.rs#L982
        let expect: Scalar =
            hex("10ee2cc4f218a2cc08c42185b0053d3ba214cec56acec60f83412a78df474d08");
        assert_eq!(state.result(), expect);
    }

    fn neptune_pallas_55_8_3_2(r_f: usize, r_p: usize) -> Spec<Scalar, 3, 2> {
        let mds = MDSMatrix::<Scalar, 3, 2>::from([
            [
                hex("2aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaac18465fd5bb87093b2d9f21600000001"),
                hex("3000000000000000000000000000000019b4f2bd072f7ea629353058c0000001"),
                hex("19999999999999999999999999999999a74f7064d0a1dd256b4f914066666667"),
            ],
            [
                hex("3000000000000000000000000000000019b4f2bd072f7ea629353058c0000001"),
                hex("19999999999999999999999999999999a74f7064d0a1dd256b4f914066666667"),
                hex("3555555555555555555555555555555571e57f7cb2a68cb89f906e9b80000001"),
            ],
            [
                hex("19999999999999999999999999999999a74f7064d0a1dd256b4f914066666667"),
                hex("3555555555555555555555555555555571e57f7cb2a68cb89f906e9b80000001"),
                hex("24924924924924924924924924924924a5df32d92a0bce3574bacf8092492493"),
            ],
        ]);
        let constants: &[Scalar] = &[
            hex("0a61a8defbacca36e4537ff2c84fa66ceee67c9645ac27346e72ab842b9d3f15"),
            hex("21e9cefa24b89d09f91b1e8a45df275b17292b4e1aaa49301e234771128165d5"),
            hex("02255d26879b6d0d49762a88d8d5d0879f01a3b06fa57c240a8119099cb02ea3"),
            hex("27aa46e0263ddc662845f98eb5db8645f2e8baeec680ddfc16f4dd3d288d8c5b"),
            hex("19f00d6d8a121f271c14bc471fdd3bc901907fb5e36e40f96a772c143cfb0ca0"),
            hex("3499f695552f2c56f23e96df7f8f3c4dd3be0cf50259f3ef7e65d565609eafe6"),
            hex("236d5cf282247aaf3941406ecdb907076f45f12f89b007168f7c7a2b0157e0f5"),
            hex("21a58e9955e3645b224d031a69ef9b23678bd7e144dd9214cf29a71c0a8426bf"),
            hex("1d5284a657f69b5c6ae0b7e714cb38bcd175ff96a11f8bd9c77d3a2f1613a743"),
            hex("3264d05b6c40c40c0d487222b6f385f69317a90cd6793532036125902138f079"),
            hex("1f00c3b0f6ebf427238bcaf2d615ca05e8cb7b9a6eb7c667d70587f1e7aee758"),
            hex("3ef895995e82a276f18232d1f5ad341f46e4119098b175038e540a80d926d400"),
            hex("106150f58dccb090de9cfee7219a9448af8127f93ba645903a93a4dc666f5808"),
            hex("3117bbb789d8a521034590fe35014bcaa3be1f4cd443985574249e62d15cf96f"),
            hex("206c2444abbaf221cc01b4d06b0915c3ab384102429c4c946846f58c08e45d2e"),
            hex("384c77b79a3f1cc004b6d28122541c92bf3f7ade7172c92b06f24e236980012f"),
            hex("1d641158354425746db2c72c52f05d7935191258c6fc7c419a6c29b6da32c1a3"),
            hex("02c302e7c79800112c508bffb262bea2e38a1eccabbcfa18ea7c313d50266e6f"),
            hex("2118be89bcf376e2cbd7741faa335b572f125299604678eca035e35114759ab4"),
            hex("33652630502c8a44b9d65f5e9815fe36c31023b4d8b7c3639f59447024b9f2cf"),
            hex("3f86a1801f46d7f30a14daa6fce4a6589986fad1b8ec46d42799bf8985852053"),
            hex("05a46f165678e990fbde85c5e5fb483ef720160b08b9eb343d2ce7db33d2a857"),
            hex("3b8f9db5bf15e2dc98d96bb31bc7c0d011d23fa38b1079ce1e0aa0a2177246c4"),
            hex("161e9f7d4ef018f47889bd6497d7f93c6b24d9d232dfc84b50315e6f520f17bc"),
            hex("264a98f297e5b58bf21a27d1c865dc75c914d49ede4b27168b365da04eb4152b"),
            hex("38802665b39df51c986bed3c9a386d77643c7a233c75eb2aa0d5dbef34c7a11e"),
            hex("2ca0307ccf7f7b3f7ad993a80b2973e375fd28c689ba05d90a63a7c908629be8"),
            hex("1f3072135ccf376c85d03d9b56361c721422706aba7c8612c9b647ff9c4035ea"),
            hex("3ce7859d585daada6ef747bc68544191c9eb9a59f79703d473424ee64276f9db"),
            hex("2d993835b6108c60aa78d1f8b241188a77256399101a131b848212b2db53c67b"),
            hex("1fdadd9c8da406323c5849b3989781fdd0438a961b2eeb7e079a37c72cd24cdb"),
            hex("359359b7346f649bb44b8b636523abad631b76a0de563282393a6c147cf4acab"),
            hex("15675df989854f06371913c37bdc1c006814e9cf4ac6124b11aeb774d42ac2b6"),
            hex("3670c3a64b4ba2694ca7e8f75ab7f7758e511a30e570d5a9f4e83afcf50f4ee0"),
            hex("37cf10985d6ba77eef963771323962a435cc9a0f3b3e3a8c11a5ac9e39db1395"),
            hex("27e4538f90ec0b009d0fba99491c34bf33d48fc90e3f1e331e31446a160192b7"),
            hex("1b711906c22ac9937281cd1ae0f2e82828f10f229eb5b027ab11f05968d9eb6f"),
            hex("0e0aeb0f0e759a2bbde9c7ca41ebfa2f3fdd3de93f4a1595195ed82b9b2f055e"),
            hex("29839536e73ae40362f720bda04e4e56f167d6d1d054109a6f4f11466086efdc"),
            hex("0dd33054f93fa501add5c7bd132146cce38b708774c505c93b604f3a5aa931d0"),
            hex("0734e8779b78f832b20642577c774662bde5b65af599e48382a6a2fbd555a72f"),
            hex("2b3ed2067e88294cacc8bca7b14dd4d15750b5ea5f076ca7c41854f5cd16a79d"),
            hex("1f6005f1c1e35f721df557808f34d2cb1ef353729be365a5ba3e679b2e928fba"),
            hex("29fcc8a49f9ace62888ad3b41c7cb12999d3a583dc70d473c2ce7d8871f179a0"),
            hex("31c4d94ce026bad4959787e2181dcc89a61462a70c1e175750528998cf81e018"),
            hex("3f26b82b525e9e40b511580861cc406bb79a074b8a98b32f718ba2c11362f807"),
            hex("0c22db7080672b164157ab4063e0fcc66c1e221b6e85212292cc6f22772459f8"),
            hex("1757f9fac9967eedbbb9a3256e3535d67c1e269a14a8fb091ace65e23661ac39"),
            hex("376de57ed93de38cd2f30a79fb29b92e4e144fcb881ef39426b8c98dd993f8ed"),
            hex("0dba6d221e3690786243d96ea8c49e84533b982a05d2fa88175bbf19c7c4c1d4"),
            hex("1da3918f5cd423cba96924cfad5663bea9752c74f2931f3b468f9851606155e9"),
            hex("2ac84677966e174c8c71dfe13fbfae2aee0e4d88b3c54135286526cfd9cb57bd"),
            hex("06fbb5301e4cffa534c8135147704fd7e72a364858597eb3483d4eafcaa9c226"),
            hex("04f4473b814bef668e7e4a22010dba7874e44c34da109c5acec3750c3ffb9365"),
            hex("207c2cba5430a95a8f3bd861cb5b0ee38b643eadf8bf659f7c765b067a1cf3cd"),
            hex("05969d41cc3bbaaa2e9d4d05b337149453c2028c1d1a4a71bebdd8b9ca1aac95"),
            hex("19548ff7a77670d86a925c8c3ed6d343f60337ee3bf60b37c3eee2bfad77aa9c"),
            hex("19cc9df4e33ed7741e3b72596a0489968be9253c3981e851727c1342eea142ff"),
            hex("1cee1165c26e9a5698fcb59ee514b2dd2921f6dfbb7d0b51cd841e09beaa7e97"),
            hex("28ad294f58c5c3e4b9c1aac0f1a014d7a7e0a03cfe4d42585e7b7036a6219619"),
            hex("3592ab44a12ea16648299f14738ea0817e1b563d3b8000a39cb16eb18b76414a"),
            hex("375513f5f15360c7f90824119b64d79f4e95a792e79bab0ae960235b33d34d13"),
            hex("16a315ece8d14bf6e327e1c2aa7f9764367fb501470a428cb43cbda39895ef7a"),
            hex("2af548fbe4c0e02176091ffbb6f3a1b9d1f313dd66643eebdeb9126c7e897ee8"),
            hex("0d58dbde546445276d5029bae57752e71ab406d91046f32fdb312592a8d7c86a"),
            hex("1de64bab1cae8dc35d54bec400e5f25128ad37a1069f8fce553b03ad0ed33fcc"),
            hex("2fa34fcd7ac076e02e9679ffd214ce9d68e180184449717a5926b68d82ef116a"),
            hex("0f2ec24a54a31240f34724ea1476973ff062af17d068b27512e241ce51351030"),
            hex("29d8d55cec33000e31b019a0dd8cb14bf9d915a75570bf5f9320b6c5d9d0da0d"),
            hex("37a368626670937de8b35c912a6b82f4aaf60ab63803649802434ca27a071eed"),
            hex("32f8ef20df6c486572088b70537db5ceaeb79ccd643e9e1c0290bb19350f572c"),
            hex("009df831fa9cedacaa065f7bdd5581482e6aad615eb472fd5850a4d0c7477e0c"),
            hex("227dee5b1c7f8c409f38d4980b01355d7b21d12e0610ac1317ac724296f8dc2b"),
            hex("2527405b6d64e04172b68288cda1cb70f848ce6e831d84558b8be56b19d9db33"),
            hex("3dfa278c4f8d8744448747e6340e908bad61b520d15377262173b78d6d3547b8"),
            hex("28470d01333eeb1857e7595d0778318d986c435b8c5c094e225b7d11ca1174ac"),
            hex("15fd96cd7931143ec393ba38b3b15e689b611dd3b99196677260cb3f7f8db9c0"),
            hex("3101540e1d7f0a923e53e1ef86def6b0a86e9f036a004be5ddde44df43fd7706"),
            hex("22e2471398685d19af0f408fd1b18cd9442ed279fffdfee27bcfe39d8b8f236e"),
            hex("35b59f926a693af114bddce789fca3b67c458e224e82be2aa914c79ed039ad85"),
            hex("0b2d458ee73f4e6054d60365fe40b6120861d7bef7c63c25204577a73448ba0c"),
            hex("02f95cf5081a78c48e317874ab8670a7075da766fc20154f99368269de3b94b5"),
            hex("307b41c27156ac0be05eca41a102dc8823b3ddbc1399dc8ce841c1f933caf041"),
            hex("02fc9786ab1b448457ac98c9314454dfc6c638ae8481c3200ecd17dd31d253a5"),
            hex("39176c107168eb6620ebfefff311bfe280abe3f48f8751485a1db30deaad785d"),
            hex("204f0ac20ece2cc11169e30799579971c5f02234708b45ea5f066214cccb06a6"),
            hex("15a6df21167207003961190a8236d7e7dcd6268085760c8563962e6f87d09acf"),
            hex("2c0e9f70303690b7e14cb877cd268d9781e0050668e9665f04ea09f94ab0817c"),
            hex("231dda54ee054a0afc0530a3064ff80b787e0b25d99596e37ba53aeb3cdffd79"),
            hex("274b78559e51e5426dc8ff296a065973bee858b2434639f66f58d466f9d99af4"),
            hex("2a0c027d819603f0735178ad5f7bf87ca5b9387ed0a3b7a8c1b7228c8e335e54"),
            hex("00fb6ba60efe32561f1d9d57f69c56fca34c01cf6a3b1ab9ccbcf2a2e647d1dc"),
            hex("2de99c53670cdb1ed5436b55ce0c1fe82628e2897b39f3c8f47fbc7ce454fc0d"),
            hex("35a1bad35f8a78a9e422036bc48dbf999bccbc5baf65a9867df972615989d00e"),
            hex("2784e7ef9ac462b56e24d1117b25343e45e37a625e3491614513030aa7f18968"),
            hex("11bcab21e844da9a2a9da9721a1140213ce8d1f5adbadf22d3bdb0393ec5a7e9"),
            hex("1aa8afc39349a7c2dd1dc466421a8891e04c456cc781e0f8bc7d64462213f45c"),
            hex("1cdba5a6ff825aff032ee88da1cc50cc4f4e71d012107432c172ac9894c225d8"),
            hex("15090a6658f804a428093e2f98736de433abe87a4b1dc8b2466a3733a62c6787"),
            hex("3a1864e0ca051a15ce8c00eae2043a7faccbbc2bc0c94629667ca55db404f474"),
            hex("2e2b4cba4fdd3698a9047f1db8a13338220e2777f5a1fd5bb83df6a7c4a649c3"),
            hex("3d316650b675ecb72159173b507f73db67e700372a8943ae3965da155ed9d012"),
            hex("08d5b91751a690a7545406c32eb2cc458a4611741c12a4c3e4b2ef5865202d92"),
            hex("1662a1e3946393e6cbd798b894f344ace1d2071489ad43e6d414e393e5baecc1"),
            hex("02f7b1099752ece5dd477e5e8d5d670d49fdba7cce0e1ef8510218a63de10aab"),
            hex("12768a65354495aac5e9331749539815747a065739e7626ba8df881b2d9324d4"),
            hex("3be432f5ad1432ac32aadcd4b72fc82e3be92e6a7e9a616cb0e0d0647cf18ca2"),
            hex("30c8f2e7c2779bbbd3b488077a65fc2ab45f2e7a4c2e5df46df4c57f5f2cc426"),
            hex("11ae3074dff2586f44849500f1194adf58afade0fc11d2058a4e3e93055adb15"),
            hex("1e043ac32449620006757805064d9582b3bb9741bf9822a5aa8748fa93e99bf6"),
            hex("1620a9895b24bf187c3cffdebf0d0a0b396716b978f4d8717e7c6a0f27d99297"),
            hex("0d314541a511fb03b53dac3ebc3d4ab633be5f2d055068f9ca2b02a6d67f60cd"),
            hex("12c0e69fefd4681508a8fa231d7a9726e1b5d774d4e96dfa91dc4a6dd721652b"),
            hex("02baa375ffe20b8935d97dc15148e4983fbdea18d7686dafb2532351894efe21"),
            hex("0af9eaaf42e0471c746350080980ac46b87130ae07312713233f6ebc52d497ec"),
            hex("3be4140e6cbe6f140d0a21fd77655d95bedbdee13267c380028d5aaac780f426"),
            hex("2adfc9f3088a4fb13abf3f0f7d662f6586df8d4f5acf0375850ba497295370b0"),
            hex("04ca9a07f4f4213ee91a6a4841a5e9e667e6a1a548d93ef0fb83736aca4d742e"),
            hex("01bc1f61aa7885a97b8a4d866897af47f5badb9e23e2aaa0375ae4f3550e5b4d"),
            hex("06235273e6aa4b38e9494ec78790e54f2098b2c4fe303e7ddd380539c85101fb"),
            hex("0f87d3b2ca4127ab6ce8dea7c1684871224b2ebf5bee8a335d233862caf97943"),
            hex("189291386264653eaf7613561cba1fa2e48ece765e6b02c1487949ee076a44ee"),
            hex("19eec3a14ff37a8314c64a9da7500329477cd106c381f888dfe0aa827f5b03cc"),
            hex("33dac6273285f29c40cd24ad24cbe1bfe81989a4ed17d2eda7f62f97a78cfd5c"),
            hex("2285f2061b71d339f94e7c9c1ca0c9159b53e583f7ceac91e0efb0644aa0c167"),
            hex("27676e54c3c06eb3ae2a5643e75218ff6569a9cfac9608f8a24960530efeef55"),
            hex("3836db482fb32f40c497ab3f794f7ac8f823e12da2bbeaa1d2f98b3b76389e12"),
            hex("2d11c4a9627ccd56a928232a08e738cdf3fb5e66bd8448a4c8d68591ef0a3499"),
            hex("0fdc0f588f62457ad4c296e5ab550e88076ceaf9cc55efe237cb0a0dbc4ed4f8"),
            hex("32429631cb6c4db1861f200f41559c6b4bf5be6438b5bff86cff68bc82aad3db"),
            hex("3b77098988c869548355e8bee5f218449fa02fd827bedb3cf4c18e45bb42704c"),
            hex("1735e3aa0cc0559603bee805b25456f9cef8d8911527fbcbea9554d814f2e686"),
            hex("30c05a5468568ef7ed2e889b8fad4b73239b0df6761da49ceea8ef58194ecd6f"),
            hex("2dd6e1e17bcc08bbcd712c26c3c19a7be009c6bce752e19f001f22c4da0731c6"),
            hex("269db0b70b999d11aa397b3e6ec274817fa8026cbc74c27fa3dfeb9b98c8c9a4"),
            hex("191d49d4da2a3d62ae722e318fe8b950b87ab1a82761fe7a452938c38c33a863"),
            hex("0ca2b6a600230cc9272a7fbbf378ac7332619e820c34a3e8bc56af7cecf684a8"),
            hex("053b5c9cf4def900ee1d329dc02648c043779f9f6aa89a050cd773dba622af9f"),
            hex("21ee54e7e62fecbbaed5bca466222f10385bdfee53a99ffc3787087a26946374"),
            hex("32e253d1bdfe40724140769883523beb8b18ca4b72250639bee7e3b11ba663b7"),
            hex("211de7393e8671c36f686b24dc86f19434f437beb5e1ccd118040e7ea5c103c3"),
            hex("0379d31b265da71af9a641ae545a9bc606454c3c01827b3f19282a5d5042d453"),
            hex("2be704bd24c29362d1a3b545faeafda0f26829b4193768af1bb418dbfb780faa"),
            hex("2480270ae8972dbb237187658b6025fd524fe0ad96509c8bd484af3f0f2e08de"),
            hex("35278936efe03adae52ab4b5ad992303a9905c3c13052cf7d2e7d056f360943d"),
            hex("288ab45752b11b0f7f712524e6ed123213d385fc6ad29dfb8f68c8b845c15522"),
            hex("33c466a10c364d1731408a6fe2a9e32ddc21ac93c0f1c54bbce4b6812f769e16"),
            hex("0d3cbc57306d8630b3374678b07ef62978e3e8c9b4dfe97077e4b85291949046"),
            hex("2795dfd778aebb427281424a31fc44bcb069a3e42e9c8118275632a8e7d8fdc1"),
            hex("27a4244d6875b63d4f7279111c71a4f6779f832d686c9789e76e86b58311c2da"),
            hex("329ebbfac553eec1be83e3fa4cb77a55fc09972997b9f73f1ba17015798de39e"),
            hex("257d992bf6c1f489bd563ea4842c772d19d3e3eaefdbc7eed1f7084cac3df5cb"),
            hex("1db4e47e8ba32fd203456816509b868a5d35fece6758e6dea8ce885feb2eadf6"),
            hex("1fb6c7a3bdfe5d7c896f7472654c711b9317abf2aa539e9d02c216c602728e1e"),
            hex("2343cf3afdc6f6094661c63a50a86630963f44b4b9ef81f82a0cd7b6289df8e7"),
            hex("09ce53dd8c572396566508e8ccddd9902ebb07b8454fa07c63e5cfd1f2986730"),
            hex("16d3f0ddbc2fd48c5c386626f0eb51f670a032582c1dc4d40e5cd7ed220a2e06"),
            hex("305269a3a6f6926ec1f5b70fa216016aa42a8a73da5421861e1b251f4c568abb"),
            hex("0e0102d913fff921047d137789321bab22181583c4cef84415ce6df2fcd844f4"),
            hex("2ad1a118294d7c6605ed0b49465df6c130ecdc199f1801b5be401d005cff9730"),
            hex("035898d75c0739ff36400fd9e1fbbf9f83bd215d3acff8edc3b0d9e1338a9a09"),
            hex("23bc2773cd5a858b0ae4dd30366567e8a1e5a64ce5c77eda21abc86e60ba547a"),
            hex("11dccb69a2eb42171c7249bcc6a3d4d4285d225558cce147c85f8221fbe8df38"),
            hex("12221ab31fe2148f5217e460e634ecf61b69fff150dd2e2044f9777274cb223c"),
            hex("114b366ab7fb2731a5e83f538ca315e49f52e71eca63bc4a953cdf14a3b3678b"),
            hex("1cfa8cf2006372ca7255dc327ccdc26c2c1d65b7161469ec419244e1597bf6ae"),
            hex("2d92ce61dc4bbce6b0c573fc6bf7370782614000477b241aa95b9af5d3cfa32b"),
            hex("3bad995c6fd0ae47e6a91405953945438b958fcc08021caf61fa8131179cd76d"),
            hex("3ba3af8f95e6f509799bc74228b4b6161660a0affa25fcf91bef65f9bf381b20"),
            hex("3063df8f291478469619b6dada5bd422fe1b8308019adc2a2ef1f3ad0348c67b"),
            hex("1666980303268177d51e0619e70ce7b49c21bc9cf28e1ecd12fe54154f124feb"),
            hex("2f8447d8fd45ac7186c08eb1e1d829f8a4cc67202fdf1d0b003f751eff8db3a3"),
            hex("1c66830e779755e38f815677c8e8e520fce9308c57c679a884a9ed4d9b6ebb90"),
            hex("3cc920c6dcd7cf5064965d07f6ed6969a87ea3869c57ee0fb96901a679eda0e0"),
            hex("21e2236fe90d17af96edb1aeab5301a18bb2e75804ac1b47afcfd5fbc8f90dce"),
            hex("1a6526263fa811c9f6530078f1fe72325f985a2926a3f44111aaa38a0dd21593"),
            hex("0fefc6aa310e9a63734049e8a76af7f943d76a1ab509505dc384dd7db04c2da7"),
            hex("186e39c5a90623656345d3150f17c92ede61eb9fb4d2aa52a1e75c9d1e27745c"),
            hex("304ef3e92ed0c929d8a16dea05d6b9dfb7fb1747fd3a5547cccaeba3a8701d63"),
            hex("3cfcf53bd0a2ba563e4595326b921cbedafc53be7ea0a22f9e36d07d1e81f770"),
            hex("046829d581cbe1e566055418ba4fc39ceb83815dc530f34af318e95c5fc0fabb"),
            hex("3540a5f1b2c6f6c45c022e2b8aac05f872c39ae4ad06e8451f763e77e5b3a213"),
            hex("35268d70cc609fa4e497678c1dc56f0735837617b8de58a7bab8062286efef73"),
            hex("222520ae7d3e5a507b7e33423f2da4b96b5ca761d7d330ccf80c6ccdc7680415"),
            hex("076fec06cd217955b375dc1c28ae5a0848ec7979c98e92381b443c6775c8a045"),
            hex("04dee46810b7dd8a74a4dd680176925723873c4cbc34ed5196af7f2cac519199"),
            hex("04debf311719afbf83f2a7a2b8af1359541c2174b86209f9fe43ba102907127c"),
            hex("137ac1f217a88d0a93dc5b5d9464e8e337a379a468d508d30815f1b6a5ca2a6f"),
            hex("2d79bb0ff9cf5c8931c9ea1259edc69aef16a707f933511cc6b821cf4d0f490e"),
        ];

        Spec::from(r_f, r_p, constants, &mds)
    }
}
