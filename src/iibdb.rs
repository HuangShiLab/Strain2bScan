//! Reader for Fast2bRAD-M binary tag files (`.iibdb` / `.iibsp`).
//!
//! Record layout (from Fast2bRAD-M `io_utils::write_binary_record`):
//!   `hash: u64 (little-endian)` · `id_len: u16 (LE)` · `id: id_len bytes (UTF-8)`
//!
//! Because the marker hash is **stored** in the file, we consume Fast2bRAD's own values
//! directly — there is no need to re-hash tags or reproduce its FxHash. This is the clean
//! interop path: build the species genome tag sets from genome `.iibdb`, and sample tag
//! counts from read `.iibsp`, all sharing Fast2bRAD's hashing.
//!
//! Genome-DB record ids look like `GCF_xxx|tag_index|scaffold|pos|0|1` (see
//! `build_quan_db.rs`); raw digests use `seqid:pos`. We expose the `gcf`/first id field so
//! callers can group tags by genome.
//!
//! NOTE: gzipped (`.gz`) inputs are not handled here (no flate2 dep in this prototype);
//! decompress first, or add flate2 in production.

use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufReader, Read};
use std::path::Path;

use crate::markers::Marker;

/// Streaming reader over `(hash, id)` records.
pub struct IibReader<R: Read> {
    reader: R,
}

impl<R: Read> IibReader<R> {
    pub fn new(reader: R) -> Self {
        Self { reader }
    }

    /// Read the next record, or `None` at clean EOF.
    pub fn next_record(&mut self) -> io::Result<Option<(Marker, String)>> {
        let mut hash_buf = [0u8; 8];
        if let Err(e) = self.reader.read_exact(&mut hash_buf) {
            if e.kind() == io::ErrorKind::UnexpectedEof {
                return Ok(None);
            }
            return Err(e);
        }
        let hash = u64::from_le_bytes(hash_buf);

        let mut len_buf = [0u8; 2];
        self.reader.read_exact(&mut len_buf)?;
        let len = u16::from_le_bytes(len_buf) as usize;

        let mut id_buf = vec![0u8; len];
        self.reader.read_exact(&mut id_buf)?;
        let id = String::from_utf8_lossy(&id_buf).into_owned();
        Ok(Some((hash, id)))
    }
}

pub fn open(path: &Path) -> io::Result<IibReader<BufReader<File>>> {
    if path.extension().and_then(|e| e.to_str()) == Some("gz") {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "gzip .iibdb not supported in prototype; decompress first or add flate2",
        ));
    }
    Ok(IibReader::new(BufReader::new(File::open(path)?)))
}

/// The genome/sample id of a record = the field before the first `|`, or before the
/// first `:` for raw digests (`seqid:pos`).
pub fn record_owner(id: &str) -> &str {
    if let Some((head, _)) = id.split_once('|') {
        head
    } else if let Some((head, _)) = id.rsplit_once(':') {
        head
    } else {
        id
    }
}

/// Load a sample's marker counts from a read `.iibsp` file.
pub fn sample_counts(path: &Path) -> io::Result<HashMap<Marker, u32>> {
    let mut r = open(path)?;
    let mut counts: HashMap<Marker, u32> = HashMap::new();
    while let Some((hash, _id)) = r.next_record()? {
        *counts.entry(hash).or_insert(0) += 1;
    }
    Ok(counts)
}

/// Load genome -> marker set from one or more `.iibdb` files (grouped by record owner id).
pub fn genome_marker_sets(
    paths: &[&Path],
) -> io::Result<HashMap<String, std::collections::HashSet<Marker>>> {
    let mut out: HashMap<String, std::collections::HashSet<Marker>> = HashMap::new();
    for p in paths {
        let mut r = open(p)?;
        while let Some((hash, id)) = r.next_record()? {
            out.entry(record_owner(&id).to_string())
                .or_default()
                .insert(hash);
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    fn encode(records: &[(u64, &str)]) -> Vec<u8> {
        let mut buf = Vec::new();
        for &(h, id) in records {
            buf.extend_from_slice(&h.to_le_bytes());
            let b = id.as_bytes();
            buf.extend_from_slice(&(b.len() as u16).to_le_bytes());
            buf.extend_from_slice(b);
        }
        buf
    }

    #[test]
    fn reads_records_roundtrip() {
        let data = encode(&[
            (42, "GCF_001|0|scaf1|100|0|1"),
            (42, "seqA:200"),
            (7, "GCF_002|0|s|9|0|1"),
        ]);
        let mut r = IibReader::new(Cursor::new(data));
        let mut recs = Vec::new();
        while let Some(x) = r.next_record().unwrap() {
            recs.push(x);
        }
        assert_eq!(recs.len(), 3);
        assert_eq!(recs[0], (42, "GCF_001|0|scaf1|100|0|1".to_string()));
    }

    #[test]
    fn owner_parsing() {
        assert_eq!(record_owner("GCF_001|0|scaf1|100|0|1"), "GCF_001");
        assert_eq!(record_owner("seqA:200"), "seqA");
        assert_eq!(record_owner("plain"), "plain");
    }

    #[test]
    fn sample_counts_accumulate() {
        let data = encode(&[(42, "a"), (42, "b"), (7, "c")]);
        // write to a temp file
        let dir = std::env::temp_dir();
        let path = dir.join("strain2bscan_iib_test.iibsp");
        std::fs::write(&path, &data).unwrap();
        let counts = sample_counts(&path).unwrap();
        assert_eq!(counts[&42], 2);
        assert_eq!(counts[&7], 1);
        let _ = std::fs::remove_file(path);
    }
}
