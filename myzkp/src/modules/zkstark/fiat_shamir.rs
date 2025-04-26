use bincode;
use serde::Serialize;
use sha3::{
    digest::{ExtendableOutput, Update, XofReader},
    Shake256,
};


pub type SerializedObj = Vec<u8>;

pub struct FiatShamirTransformer {
    objects: Vec<Vec<SerializedObj>>,
    read_index: usize,
}

impl FiatShamirTransformer {
    pub fn new() -> Self {
        FiatShamirTransformer {
            objects: Vec::new(),
            read_index: 0,
        }
    }

    pub fn push(&mut self, obj: &Vec<SerializedObj>) {
        self.objects.push(obj.clone());
    }

    pub fn pull(&mut self) -> Vec<SerializedObj> {
        if self.read_index >= self.objects.len() {
            panic!("ProofStream: cannot pull object; queue empty.");
        }
        let obj = self.objects[self.read_index].clone();
        self.read_index += 1;
        obj
    }

    pub fn serialize(&self) -> SerializedObj {
        bincode::serialize(&self.objects).expect("Serialization failed")
    }

    pub fn deserialize(bb: &SerializedObj) -> Self {
        let objects: Vec<Vec<SerializedObj>> =
            bincode::deserialize(bb).expect("Deserialization failed");
        FiatShamirTransformer {
            objects,
            read_index: 0,
        }
    }

    pub fn prover_fiat_shamir(&self, num_bytes: usize) -> SerializedObj {
        let serialized = self.serialize();
        let mut hasher = Shake256::default();
        hasher.update(&serialized);
        let mut reader = hasher.finalize_xof();
        let mut result = vec![0u8; num_bytes];
        reader.read(&mut result);
        result
    }

    pub fn verifier_fiat_shamir(&self, num_bytes: usize) -> SerializedObj {
        let serialized =
            bincode::serialize(&self.objects[..self.read_index]).expect("Serialization failed");
        let mut hasher = Shake256::default();
        hasher.update(&serialized);
        let mut reader = hasher.finalize_xof();
        let mut result = vec![0u8; num_bytes];
        reader.read(&mut result);
        result
    }
}
