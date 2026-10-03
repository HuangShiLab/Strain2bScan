//! CLI regression tests for the clap-derived parser.

use std::process::Command;

fn run(args: &[&str]) -> std::process::Output {
    Command::new(env!("CARGO_BIN_EXE_strain2bscan"))
        .args(args)
        .output()
        .unwrap()
}

#[test]
fn top_level_help_succeeds() {
    let out = run(&["--help"]);
    let stdout = String::from_utf8_lossy(&out.stdout);
    assert!(out.status.success(), "--help failed: {stdout}");
    assert!(
        stdout.contains("Usage:") && stdout.contains("strain2bscan"),
        "help should show usage: {stdout}"
    );
}

#[test]
fn profile_subcommand_help_succeeds() {
    let out = run(&["profile", "--help"]);
    let stdout = String::from_utf8_lossy(&out.stdout);
    assert!(out.status.success(), "profile --help failed: {stdout}");
    assert!(
        stdout.contains("Usage:") && stdout.contains("profile"),
        "help should show profile usage: {stdout}"
    );
}

#[test]
fn value_starting_with_dash_parsed_via_equals() {
    // A value starting with `--` must be passable as `--reads=--sample.fastq`.
    // The DB is missing, so the command fails, but it must NOT fail because of an
    // unknown argument.
    let out = run(&[
        "profile",
        "--reads=--sample.fastq",
        "--db=does_not_exist.tsv",
    ]);
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        !stderr.contains("unexpected argument"),
        "value starting with -- was misparsed as a flag: {stderr}"
    );
    assert!(
        stderr.contains("does_not_exist.tsv")
            || stderr.contains("No such file")
            || stderr.contains("cannot find"),
        "expected a missing-DB error, got: {stderr}"
    );
}
