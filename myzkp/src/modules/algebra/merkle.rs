use sha3::{Digest, Sha3_256};

pub struct Merkle;
pub type MerkleRoot = Vec<u8>;
pub type MerklePath = Vec<Vec<u8>>;

impl Merkle {
    fn hash(data: &[u8]) -> Vec<u8> {
        let mut hasher = Sha3_256::new();
        hasher.update(data);
        hasher.finalize().to_vec()
    }

    // Computes the Merkle root of a given array.
    pub fn commit(leafs: &[Vec<u8>]) -> MerkleRoot {
        // assert!(leafs.len().is_power_of_two())
        if leafs.len() == 1 {
            return leafs[0].clone();
        } else {
            let mid = leafs.len() / 2;
            let left_commit = Self::commit(&leafs[..mid]);
            let right_commit = Self::commit(&leafs[mid..]);
            Self::hash(&[left_commit, right_commit].concat())
        }
    }

    // Computes the authentication path of an indicated leaf in the Merkle tree.
    pub fn open(index: usize, leafs: &[Vec<u8>]) -> MerklePath {
        // assert!(leafs.len().is_power_of_two())
        // assert!(index < leafs.len())

        if leafs.len() == 2 {
            return vec![leafs[1 - index].clone()];
        }

        let mid = leafs.len() / 2;
        if index < mid {
            let mut path = Self::open(index, &leafs[..mid]);
            path.push(Self::commit(&leafs[mid..]));
            path
        } else {
            let mut path = Self::open(index - mid, &leafs[mid..]);
            path.push(Self::commit(&leafs[..mid]));
            path
        }
    }

    // Verifies that a given leaf is an element of the committed vector at the given index.
    pub fn verify(root: &MerkleRoot, index: usize, path: &[Vec<u8>], leaf: &[u8]) -> bool {
        // assert!(index < (1 << path.len()))

        if path.len() == 1 {
            if index == 0 {
                root == &Self::hash(&[leaf, &path[0]].concat())
            } else {
                root == &Self::hash(&[&path[0], leaf].concat())
            }
        } else {
            let next_hash = if index % 2 == 0 {
                Self::hash(&[leaf, &path[0]].concat())
            } else {
                Self::hash(&[&path[0], leaf].concat())
            };
            Self::verify(root, index >> 1, &path[1..], &next_hash)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_merkle() {
        let leafs = vec![
            b"leaf1".to_vec(),
            b"leaf2".to_vec(),
            b"leaf3".to_vec(),
            b"leaf4".to_vec(),
        ];
        let root = Merkle::commit(&leafs);

        let index = 2;
        let proof = Merkle::open(index, &leafs);

        let is_valid = Merkle::verify(&root, index, &proof, &leafs[index]);
        assert!(is_valid);

        let is_not_valid = Merkle::verify(&root, index, &proof, &leafs[index + 1]);
        assert!(!is_not_valid);
    }
}
