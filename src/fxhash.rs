//! Dependency-free FxHash — the hash `rustc` itself uses — for `u64`-keyed maps and sets.
//!
//! Every hot container in this crate is keyed by [`crate::markers::Marker`] (a `u64`). The
//! std default hasher is SipHash-1-3, which is DoS-resistant but computes a full keyed
//! permutation per lookup; on a `u64` key that cost dominates the lookup itself. Marker keys
//! come from untrusted-input-independent digestion of the user's own reads, so collision
//! resistance buys us nothing here.
//!
//! FxHash is a single rotate-xor-multiply per 8 bytes. Swapping it in changes no marker
//! values and no on-disk format — `Marker` is still FNV-1a of the canonical tag (see
//! [`crate::markers::hash_bytes`]), so **databases built by earlier versions still load**.
//! Only the in-memory bucket assignment changes.

use std::hash::{BuildHasherDefault, Hasher};

/// Fractional part of the golden ratio scaled to 64 bits (rustc-hash's constant).
const SEED: u64 = 0x51_7c_c1_b7_27_22_0a_95;

#[derive(Default, Clone, Copy)]
pub struct FxHasher {
    hash: u64,
}

impl FxHasher {
    #[inline]
    fn add(&mut self, word: u64) {
        self.hash = (self.hash.rotate_left(5) ^ word).wrapping_mul(SEED);
    }
}

impl Hasher for FxHasher {
    #[inline]
    fn write(&mut self, bytes: &[u8]) {
        let mut chunks = bytes.chunks_exact(8);
        for c in &mut chunks {
            self.add(u64::from_ne_bytes(c.try_into().unwrap()));
        }
        let rem = chunks.remainder();
        if !rem.is_empty() {
            let mut buf = [0u8; 8];
            buf[..rem.len()].copy_from_slice(rem);
            self.add(u64::from_ne_bytes(buf));
        }
    }

    #[inline]
    fn write_u8(&mut self, i: u8) {
        self.add(i as u64)
    }
    #[inline]
    fn write_u16(&mut self, i: u16) {
        self.add(i as u64)
    }
    #[inline]
    fn write_u32(&mut self, i: u32) {
        self.add(i as u64)
    }
    #[inline]
    fn write_u64(&mut self, i: u64) {
        self.add(i)
    }
    #[inline]
    fn write_usize(&mut self, i: usize) {
        self.add(i as u64)
    }
    #[inline]
    fn finish(&self) -> u64 {
        self.hash
    }
}

pub type FxBuildHasher = BuildHasherDefault<FxHasher>;
pub type FxHashMap<K, V> = std::collections::HashMap<K, V, FxBuildHasher>;
pub type FxHashSet<T> = std::collections::HashSet<T, FxBuildHasher>;

/// `FxHashMap::default()` with room for `n` entries (avoids rehash churn on big digests).
pub fn map_with_capacity<K, V>(n: usize) -> FxHashMap<K, V> {
    FxHashMap::with_capacity_and_hasher(n, FxBuildHasher::default())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn behaves_as_a_map() {
        let mut m: FxHashMap<u64, u32> = FxHashMap::default();
        for i in 0..10_000u64 {
            m.insert(i.wrapping_mul(0x9e37_79b9_7f4a_7c15), i as u32);
        }
        assert_eq!(m.len(), 10_000);
        for i in 0..10_000u64 {
            assert_eq!(m.get(&i.wrapping_mul(0x9e37_79b9_7f4a_7c15)), Some(&(i as u32)));
        }
    }

    #[test]
    fn distinct_keys_mostly_distinct_hashes() {
        // Sanity: the finalizer must not collapse sequential keys into one bucket value.
        use std::collections::HashSet;
        let mut seen = HashSet::new();
        for i in 0..1000u64 {
            let mut h = FxHasher::default();
            h.write_u64(i);
            seen.insert(h.finish());
        }
        assert_eq!(seen.len(), 1000);
    }
}
