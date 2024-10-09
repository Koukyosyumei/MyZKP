```bash
circom sha256.circom --r1cs --wasm --sym
node sha256_js/generate_witness.js sha256_js/sha256.wasm input.json witness.wtns

snarkjs powersoftau new bn128 16 pot16_0000.ptau -v
snarkjs powersoftau contribute pot16_0000.ptau pot16_0001.ptau --name="First contribution" -v

snarkjs powersoftau prepare phase2 pot16_0001.ptau pot16_final.ptau -v
snarkjs groth16 setup sha256.r1cs pot16_final.ptau sha256_0000.zkey
snarkjs zkey contribute sha256_0000.zkey sha256_0001.zkey --name="1st Contributor Name" -v
snarkjs zkey export verificationkey sha256_0001.zkey verification_key.json

snarkjs groth16 prove circuit_0000.zkey witness.wtns proof.json public.json
snarkjs groth16 verify verification_key.json public.json proof.json"
```