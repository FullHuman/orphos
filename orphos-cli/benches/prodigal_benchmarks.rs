use criterion::{BenchmarkId, Criterion, Throughput, black_box, criterion_group, criterion_main};
use std::path::Path;
use std::process::Command;
use std::time::Duration;
use tempfile::NamedTempFile;

mod criterion_config;
use criterion_config::configure_criterion;

// Helper function to run original prodigal
fn run_original_prodigal(
    input_file: &str,
    output_file: &str,
    threads: Option<usize>,
) -> Result<Duration, Box<dyn std::error::Error>> {
    let start = std::time::Instant::now();

    let mut cmd = Command::new("prodigal");
    cmd.arg("-i")
        .arg(input_file)
        .arg("-o")
        .arg(output_file)
        .arg("-f")
        .arg("gff")
        .arg("-p")
        .arg("single");

    // Original prodigal doesn't have threading, but we include this for consistency
    if threads.is_some() {
        // Original prodigal is single-threaded
    }

    let output = cmd.output()?;
    let duration = start.elapsed();

    if !output.status.success() {
        return Err(format!(
            "Original prodigal failed: {}",
            String::from_utf8_lossy(&output.stderr)
        )
        .into());
    }

    Ok(duration)
}

// Helper function to run pyrodigal
fn run_pyrodigal(
    input_file: &str,
    output_file: &str,
    threads: Option<usize>,
) -> Result<Duration, Box<dyn std::error::Error>> {
    let start = std::time::Instant::now();

    let mut cmd = Command::new("pyrodigal");
    cmd.arg("-i")
        .arg(input_file)
        .arg("-o")
        .arg(output_file)
        .arg("-f")
        .arg("gff")
        .arg("-p")
        .arg("single");

    if let Some(t) = threads {
        cmd.arg("-j").arg(t.to_string());
    }

    let output = cmd.output()?;
    let duration = start.elapsed();

    if !output.status.success() {
        return Err(format!(
            "Pyrodigal failed: {}",
            String::from_utf8_lossy(&output.stderr)
        )
        .into());
    }

    Ok(duration)
}

// Helper function to run orphos
fn run_orphos_cli(
    input_file: &str,
    output_file: &str,
    threads: Option<usize>,
) -> Result<Duration, Box<dyn std::error::Error>> {
    let start = std::time::Instant::now();

    // Set thread count via environment variable if specified
    if let Some(t) = threads {
        unsafe { std::env::set_var("RAYON_NUM_THREADS", t.to_string()) };
    }

    let mut cmd = Command::new("cargo");
    cmd.arg("run")
        .arg("--release")
        .arg("--bin")
        .arg("orphos-cli")
        .arg("--")
        .arg("-i")
        .arg(input_file)
        .arg("-o")
        .arg(output_file)
        .arg("-f")
        .arg("gff")
        .arg("-p")
        .arg("single");

    let output = cmd.output()?;
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

// Get file size for throughput calculation
fn get_file_size(path: &str) -> u64 {
    std::fs::metadata(path).map(|m| m.len()).unwrap_or(0)
}

fn benchmark_single_threaded(c: &mut Criterion) {
    let test_files = [
        ("ecoli", "tests/data/ecoli.fasta"),
        ("salmonella", "tests/data/salmonella.fasta"),
    ];

    for (name, input_file) in test_files {
        if !Path::new(input_file).exists() {
            eprintln!("Warning: {} not found, skipping benchmark", input_file);
            continue;
        }

        let file_size = get_file_size(input_file);
        let mut group = c.benchmark_group(format!("{}_single_threaded", name));
        group.throughput(Throughput::Bytes(file_size));

        // Check if original prodigal is available
        if Command::new("prodigal").arg("--help").output().is_ok() {
            group.bench_function("original_prodigal", |b| {
                b.iter_custom(|iters| {
                    let mut total_duration = Duration::new(0, 0);
                    for _ in 0..iters {
                        let output_file = NamedTempFile::new().unwrap();
                        let duration = run_original_prodigal(
                            black_box(input_file),
                            output_file.path().to_str().unwrap(),
                            None,
                        )
                        .unwrap_or_else(|e| {
                            eprintln!("Original prodigal benchmark failed: {}", e);
                            Duration::from_secs(0)
                        });
                        total_duration += duration;
                    }
                    total_duration
                });
            });
        } else {
            eprintln!("Warning: Original prodigal not found, skipping original_prodigal benchmark");
        }

        // Check if pyrodigal is available
        if Command::new("pyrodigal").arg("--help").output().is_ok() {
            group.bench_function("pyrodigal", |b| {
                b.iter_custom(|iters| {
                    let mut total_duration = Duration::new(0, 0);
                    for _ in 0..iters {
                        let output_file = NamedTempFile::new().unwrap();
                        let duration = run_pyrodigal(
                            black_box(input_file),
                            output_file.path().to_str().unwrap(),
                            Some(1),
                        )
                        .unwrap_or_else(|e| {
                            eprintln!("Pyrodigal benchmark failed: {}", e);
                            Duration::from_secs(0)
                        });
                        total_duration += duration;
                    }
                    total_duration
                });
            });
        } else {
            eprintln!("Warning: Pyrodigal not found, skipping pyrodigal benchmark");
        }

        // Benchmark orphos-cli (single-threaded)
        group.bench_function("prodigal_cli", |b| {
            b.iter_custom(|iters| {
                let mut total_duration = Duration::new(0, 0);
                for _ in 0..iters {
                    let output_file = NamedTempFile::new().unwrap();
                    let duration = run_orphos_cli(
                        black_box(input_file),
                        output_file.path().to_str().unwrap(),
                        Some(1),
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
}

#[allow(dead_code)]
fn benchmark_multi_threaded(c: &mut Criterion) {
    let test_files = [
        ("ecoli", "tests/data/ecoli.fasta"),
        ("salmonella", "tests/data/salmonella.fasta"),
    ];

    let thread_counts = [2, 4, 8];

    for (name, input_file) in test_files {
        if !Path::new(input_file).exists() {
            eprintln!("Warning: {} not found, skipping benchmark", input_file);
            continue;
        }

        let file_size = get_file_size(input_file);

        for &threads in &thread_counts {
            let mut group = c.benchmark_group(format!("{}_{}_threads", name, threads));
            group.throughput(Throughput::Bytes(file_size));

            // Benchmark pyrodigal (multi-threaded) if available
            if threads == 1 && Command::new("pyrodigal").arg("--help").output().is_ok() {
                group.bench_with_input(
                    BenchmarkId::new("pyrodigal", threads),
                    &threads,
                    |b, &threads| {
                        b.iter_custom(|iters| {
                            let mut total_duration = Duration::new(0, 0);
                            for _ in 0..iters {
                                let output_file = NamedTempFile::new().unwrap();
                                let duration = run_pyrodigal(
                                    black_box(input_file),
                                    output_file.path().to_str().unwrap(),
                                    Some(threads),
                                )
                                .unwrap_or_else(|e| {
                                    eprintln!("Pyrodigal benchmark failed: {}", e);
                                    Duration::from_secs(0)
                                });
                                total_duration += duration;
                            }
                            total_duration
                        });
                    },
                );
            }

            // Benchmark original prodigal
            if threads == 1 && Command::new("prodigal").arg("--help").output().is_ok() {
                group.bench_with_input(
                    BenchmarkId::new("original_prodigal", threads),
                    &threads,
                    |b, &_threads| {
                        b.iter_custom(|iters| {
                            let mut total_duration = Duration::new(0, 0);
                            for _ in 0..iters {
                                let output_file = NamedTempFile::new().unwrap();
                                let duration = run_original_prodigal(
                                    black_box(input_file),
                                    output_file.path().to_str().unwrap(),
                                    None,
                                )
                                .unwrap_or_else(|e| {
                                    eprintln!("Original prodigal benchmark failed: {}", e);
                                    Duration::from_secs(0)
                                });
                                total_duration += duration;
                            }
                            total_duration
                        });
                    },
                );
            }

            // Benchmark orphos-cli (multi-threaded)
            group.bench_with_input(
                BenchmarkId::new("prodigal_cli", threads),
                &threads,
                |b, &threads| {
                    b.iter_custom(|iters| {
                        let mut total_duration = Duration::new(0, 0);
                        for _ in 0..iters {
                            let output_file = NamedTempFile::new().unwrap();
                            let duration = run_orphos_cli(
                                black_box(input_file),
                                output_file.path().to_str().unwrap(),
                                Some(threads),
                            )
                            .unwrap_or_else(|e| {
                                eprintln!("orphos-cli benchmark failed: {}", e);
                                Duration::from_secs(0)
                            });
                            total_duration += duration;
                        }
                        total_duration
                    });
                },
            );

            group.finish();
        }
    }
}

fn benchmark_scaling_analysis(c: &mut Criterion) {
    let input_file = "tests/data/ecoli.fasta";

    if !Path::new(input_file).exists() {
        eprintln!(
            "Warning: {} not found, skipping scaling analysis",
            input_file
        );
        return;
    }

    let file_size = get_file_size(input_file);
    let mut group = c.benchmark_group("scaling_analysis");
    group.throughput(Throughput::Bytes(file_size));
    group.measurement_time(Duration::from_secs(45));
    group.sample_size(10);

    let thread_counts = [1, 4, 8];

    for &threads in &thread_counts {
        // Pyrodigal scaling (if available)
        if threads == 1 && Command::new("pyrodigal").arg("--help").output().is_ok() {
            group.bench_with_input(
                BenchmarkId::new("pyrodigal_scaling", threads),
                &threads,
                |b, &threads| {
                    b.iter_custom(|iters| {
                        let mut total_duration = Duration::new(0, 0);
                        for _ in 0..iters {
                            let output_file = NamedTempFile::new().unwrap();
                            let duration = run_pyrodigal(
                                black_box(input_file),
                                output_file.path().to_str().unwrap(),
                                Some(threads),
                            )
                            .unwrap_or_else(|e| {
                                eprintln!("Pyrodigal scaling benchmark failed: {}", e);
                                Duration::from_secs(0)
                            });
                            total_duration += duration;
                        }
                        total_duration
                    });
                },
            );
        }

        // orphos-cli scaling
        group.bench_with_input(
            BenchmarkId::new("prodigal_cli_scaling", threads),
            &threads,
            |b, &threads| {
                b.iter_custom(|iters| {
                    let mut total_duration = Duration::new(0, 0);
                    for _ in 0..iters {
                        let output_file = NamedTempFile::new().unwrap();
                        let duration = run_orphos_cli(
                            black_box(input_file),
                            output_file.path().to_str().unwrap(),
                            Some(threads),
                        )
                        .unwrap_or_else(|e| {
                            eprintln!("orphos-cli scaling benchmark failed: {}", e);
                            Duration::from_secs(0)
                        });
                        total_duration += duration;
                    }
                    total_duration
                });
            },
        );
    }

    group.finish();
}

criterion_group!(
    name = benches;
    config = configure_criterion();
    targets = benchmark_single_threaded,
    // benchmark_multi_threaded,
    benchmark_scaling_analysis
);
criterion_main!(benches);
