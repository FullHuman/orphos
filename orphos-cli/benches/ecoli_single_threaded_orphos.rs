use criterion::{Criterion, SamplingMode, Throughput, criterion_group, criterion_main};
use std::hint::black_box;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::time::Duration;
use tempfile::NamedTempFile;

mod criterion_config;
use criterion_config::configure_criterion;

const ECOLI_FASTA: &str = "tests/data/ecoli.fasta";

fn orphos_binary_path() -> Result<PathBuf, Box<dyn std::error::Error>> {
    let bench_exe = std::env::current_exe()?;
    let target_dir = bench_exe
        .parent()
        .and_then(Path::parent)
        .ok_or("unable to resolve target directory from benchmark executable")?;

    Ok(target_dir.join(format!("orphos{}", std::env::consts::EXE_SUFFIX)))
}

fn ensure_orphos_binary(orphos_bin: &Path) -> Result<(), Box<dyn std::error::Error>> {
    if orphos_bin.exists() {
        return Ok(());
    }

    let output = Command::new("cargo")
        .arg("build")
        .arg("--release")
        .arg("-p")
        .arg("orphos-cli")
        .arg("--bin")
        .arg("orphos")
        .output()?;

    if !output.status.success() {
        return Err(format!(
            "failed to build orphos binary: {}",
            String::from_utf8_lossy(&output.stderr)
        )
        .into());
    }

    Ok(())
}

fn run_orphos_cli(
    orphos_bin: &Path,
    input_file: &str,
    output_file: &str,
) -> Result<Duration, Box<dyn std::error::Error>> {
    let start = std::time::Instant::now();

    let output = Command::new(orphos_bin)
        .env("RAYON_NUM_THREADS", "1")
        .arg("-i")
        .arg(input_file)
        .arg("-o")
        .arg(output_file)
        .arg("-f")
        .arg("gff")
        .arg("-p")
        .arg("single")
        .output()?;
    let duration = start.elapsed();

    if !output.status.success() {
        return Err(format!(
            "orphos-cli failed: {}",
            String::from_utf8_lossy(&output.stderr)
        )
        .into());
    }

    Ok(duration)
}

fn benchmark_ecoli_single_threaded_orphos(c: &mut Criterion) {
    if !Path::new(ECOLI_FASTA).exists() {
        eprintln!("Warning: {} not found, skipping benchmark", ECOLI_FASTA);
        return;
    }

    let orphos_bin = match orphos_binary_path().and_then(|path| {
        ensure_orphos_binary(&path)?;
        Ok(path)
    }) {
        Ok(path) => path,
        Err(e) => {
            eprintln!("Warning: unable to prepare orphos binary: {}", e);
            return;
        }
    };

    let file_size = std::fs::metadata(ECOLI_FASTA).map(|m| m.len()).unwrap_or(0);
    let mut group = c.benchmark_group("ecoli_single_threaded");
    group.throughput(Throughput::Bytes(file_size));
    group.sampling_mode(SamplingMode::Flat);

    group.bench_function("orphos", |b| {
        b.iter_custom(|iters| {
            let mut total_duration = Duration::new(0, 0);
            for _ in 0..iters {
                let output_file = NamedTempFile::new().unwrap();
                let duration = run_orphos_cli(
                    black_box(&orphos_bin),
                    black_box(ECOLI_FASTA),
                    output_file.path().to_str().unwrap(),
                )
                .unwrap_or_else(|e| {
                    eprintln!("orphos-cli benchmark failed: {}", e);
                    Duration::from_secs(0)
                });
                total_duration += duration;
            }
            total_duration
        });
    });

    group.finish();
}

criterion_group!(
    name = benches;
    config = configure_criterion();
    targets = benchmark_ecoli_single_threaded_orphos
);
criterion_main!(benches);
