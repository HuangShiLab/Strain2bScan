//! Reporting and manifest parsing for the CLI.

use std::path::{Path, PathBuf};

use strain2bscan::identify::StrainCall;
use strain2bscan::markers::is_fastx_path;

/// What the sequence reader can parse.
pub const FASTX_HINT: &str = "expected .fa/.fasta/.fna/.fq/.fastq, optionally .gz. The reader \
     selects FASTA or FASTQ by extension, so a misnamed file would be parsed as the wrong \
     format with no error at all";

/// Resolve `--reads` and reject anything the sequence reader cannot parse.
pub fn reads_path(reads: &str) -> Result<PathBuf, String> {
    let path = PathBuf::from(reads);
    if !is_fastx_path(&path) {
        return Err(format!(
            "--reads {}: unrecognized sequence format — {FASTX_HINT}.",
            path.display()
        ));
    }
    Ok(path)
}

/// Pretty-print strain calls to stdout.
pub fn report(calls: &[StrainCall]) {
    if calls.is_empty() {
        println!("  (no strains passed thresholds)");
        return;
    }
    for c in calls {
        println!(
            "  {:<12} abundance={:>6.2}%  depth={:>7.3}x  coverage={:>6.2}%  support={:.0}",
            c.name,
            c.rel_abundance * 100.0,
            c.depth,
            c.coverage * 100.0,
            c.support
        );
    }
}

/// Write predictions as `name<TAB>abundance<TAB>coverage<TAB>support<TAB>depth<TAB>n_markers`.
pub fn write_pred_tsv(path: &Path, calls: &[StrainCall]) -> std::io::Result<()> {
    use std::io::Write;
    let mut w = std::fs::File::create(path)?;
    writeln!(
        w,
        "#cluster\tabundance\tcoverage\tsupport\tdepth\tn_markers"
    )?;
    for c in calls {
        writeln!(
            w,
            "{}\t{:.6}\t{:.4}\t{:.0}\t{:.4}\t{}",
            c.name, c.rel_abundance, c.coverage, c.support, c.depth, c.n_markers
        )?;
    }
    Ok(())
}

/// One row of the `batch` manifest.
pub struct ManifestSample {
    pub name: String,
    pub reads1: PathBuf,
    pub reads2: Option<PathBuf>,
}

/// Read the batch manifest: a CSV with header `sample,reads1[,reads2]`.
pub fn read_manifest(path: &Path) -> Result<Vec<ManifestSample>, String> {
    let text = std::fs::read_to_string(path)
        .map_err(|e| format!("cannot read manifest {}: {e}", path.display()))?;
    let dir = path.parent().unwrap_or_else(|| Path::new("."));
    let mut lines = text.lines();

    let header = lines
        .next()
        .ok_or_else(|| format!("manifest {} is empty", path.display()))?;
    let cols: Vec<&str> = header.split(',').map(str::trim).collect();
    if cols.len() < 2 || cols[0] != "sample" || cols[1] != "reads1" {
        return Err(format!(
            "manifest {}: header must be `sample,reads1[,reads2]`, got `{header}`",
            path.display()
        ));
    }
    if cols.len() > 2 && cols[2] != "reads2" {
        return Err(format!(
            "manifest {}: third column must be `reads2`, got `{}`",
            path.display(),
            cols[2]
        ));
    }

    let resolve = |p: &str| {
        let q = Path::new(p);
        if q.is_absolute() {
            q.to_path_buf()
        } else {
            dir.join(q)
        }
    };
    let mut samples = Vec::new();
    for (i, line) in lines.enumerate() {
        let lineno = i + 2;
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split(',').map(str::trim).collect();
        if !(2..=3).contains(&f.len()) {
            return Err(format!(
                "manifest {} line {lineno}: want `sample,reads1[,reads2]` (2-3 fields), got {} field(s)",
                path.display(),
                f.len()
            ));
        }
        if f[0].is_empty() {
            return Err(format!(
                "manifest {} line {lineno}: empty sample name",
                path.display()
            ));
        }
        if f[1].is_empty() {
            return Err(format!(
                "manifest {} line {lineno} (sample '{}'): empty reads1",
                path.display(),
                f[0]
            ));
        }
        let reads1 = resolve(f[1]);
        let reads2 = match f.get(2) {
            Some(&p) if !p.is_empty() => Some(resolve(p)),
            _ => None,
        };
        for p in [&reads1].into_iter().chain(reads2.iter()) {
            if !p.is_file() {
                return Err(format!(
                    "manifest {} line {lineno} (sample '{}'): reads file not found: {}",
                    path.display(),
                    f[0],
                    p.display()
                ));
            }
            if !is_fastx_path(p) {
                return Err(format!(
                    "manifest {} line {lineno} (sample '{}'): unrecognized sequence format: {} \
                     — {FASTX_HINT}.",
                    path.display(),
                    f[0],
                    p.display()
                ));
            }
        }
        if samples.iter().any(|s: &ManifestSample| s.name == f[0]) {
            return Err(format!(
                "manifest {} line {lineno}: duplicate sample name '{}'",
                path.display(),
                f[0]
            ));
        }
        samples.push(ManifestSample {
            name: f[0].to_string(),
            reads1,
            reads2,
        });
    }
    if samples.is_empty() {
        return Err(format!("manifest {}: no samples", path.display()));
    }
    Ok(samples)
}

#[cfg(test)]
mod tests {
    use super::reads_path;

    #[test]
    fn reads_path_rejects_unrecognized_formats() {
        for good in [
            "s.fq",
            "s.fastq",
            "s.fq.gz",
            "s.fastq.gz",
            "s.FQ.GZ",
            "s.fa",
            "s.fasta",
            "s.fna",
            "s.fna.gz",
        ] {
            assert!(reads_path(good).is_ok(), "{good} must be accepted");
        }
        for bad in ["s.txt", "s.fastq.bz2", "s.bam", "s", "s.fq.zst"] {
            let err = reads_path(bad)
                .expect_err("{bad} must be rejected")
                .to_string();
            assert!(
                err.contains("unrecognized sequence format"),
                "{bad}: unexpected error {err}"
            );
        }
    }
}
