# Contributing to Orphos

Thank you for helping improve Orphos. Contributions of bug reports, tests,
documentation, performance improvements, and code are all welcome.

## Before you start

- Search the existing issues before opening a new one.
- For a substantial feature or behavior change, open an issue first so the
  approach can be discussed before implementation.
- Never include private genomic data, credentials, or generated build artifacts
  in a commit.

## Development setup

Orphos requires Rust 1.89 or newer. Clone your fork and create a branch from
`main`:

```console
git clone https://github.com/<your-user>/orphos.git
cd orphos
git switch -c <short-description>
```

Build and test the Rust workspace:

```console
cargo build --all-features
cargo test --all-features
```

Before submitting a pull request, run the checks used by CI:

```console
cargo fmt --all -- --check
cargo clippy --all-targets --all-features -- -D warnings
cargo test --all-features
```

When changing a language binding, also test that binding:

```console
# Python (requires maturin)
cd orphos-python
maturin develop
python -c "import orphos; print(orphos.__version__)"

# WebAssembly/Node.js (requires wasm-pack and Node.js)
cd ../orphos-wasm
npm run build:node
```

Add or update tests for behavior changes. Changes to gene prediction or output
formats should be checked against representative biological data and, where
appropriate, the existing snapshots and reference annotations.

## Pull requests

Keep each pull request focused on one change. In its description:

- explain the problem and the chosen solution;
- link any related issues;
- describe the tests you ran and the test data used;
- include before-and-after benchmarks for performance-sensitive changes;
- call out breaking changes and document how users should migrate; and
- update `README.md`, package documentation, and `CHANGELOG.md` when relevant.

Use clear, imperative commit messages. Conventional prefixes such as `feat:`,
`fix:`, `docs:`, and `chore:` are encouraged when they make the history easier
to scan.

By contributing, you agree that your contribution is licensed under the
project's [GPL-3.0-or-later license](LICENSE).

## Release process (maintainers)

Releases use [Semantic Versioning](https://semver.org/). The steps below assume
the new version is `0.3.0`: use the version without `v` in package manifests and
the pre-release workflow, and use `v0.3.0` for the Git tag and GitHub release.

Publishing to crates.io, PyPI, and npm cannot be undone or overwritten with the
same version. Complete the validation and dry runs before publishing.

### 1. Prepare the release

Start from an up-to-date `main` branch with a clean working tree:

```console
git switch main
git pull --ff-only
git status --short
```

Choose the next version according to Semantic Versioning, then update every
version source:

- `orphos-core/Cargo.toml` — package version;
- `orphos-cli/Cargo.toml` — package version and the `orphos-core` dependency
  version;
- `orphos-wasm/Cargo.toml` — package version;
- `orphos-wasm/package.json` — npm package version;
- `orphos-python/Cargo.toml` — package version; and
- `orphos-python/pyproject.toml` — Python package version.

Regenerate the npm and Cargo lockfiles after editing the manifests:

```console
npm install --package-lock-only --prefix orphos-wasm
cargo check --all-features
cargo check --manifest-path orphos-python/Cargo.toml
```

Confirm that `Cargo.lock`, `orphos-python/Cargo.lock`, and
`orphos-wasm/package-lock.json` contain the new local package versions. Do not
edit generated `orphos-wasm/pkg` or `orphos-wasm/pkg-node` files by hand.

Move the release's entries in `CHANGELOG.md` into a heading with the version and
release date:

```markdown
## [0.3.0] - YYYY-MM-DD
```

Summarize user-visible additions, changes, fixes, and breaking changes. Commit
these updates in a release pull request and merge it into `main` only after all
required CI checks pass.

### 2. Validate the release commit

From the release commit, run:

```console
cargo fmt --all -- --check
cargo clippy --all-targets --all-features -- -D warnings
cargo test --all-features
cargo build --release -p orphos-cli
```

Run any binding-specific checks relevant to the release. Then, in GitHub
Actions, run **Pre-Release Validation** on `main` with `version` set to `0.3.0`
(without `v`). It verifies version consistency, Rust builds on the supported
operating systems, the WASM build, the Python build, and the changelog entry.

Before continuing, verify:

- the validation workflow and all CI checks are green;
- the release commit on `main` contains every intended change;
- the package metadata and changelog use the same version; and
- the working tree is clean.

### 3. Create the tag and draft GitHub release

Tag the exact validated commit:

```console
git switch main
git pull --ff-only
git tag -a v0.3.0 -m "Release v0.3.0"
git show --stat v0.3.0
git push origin v0.3.0
```

Pushing the tag starts the **Release** workflow. It creates a draft GitHub
release and attaches CLI archives and SHA-256 checksums for Linux, macOS, and
Windows. Wait for every matrix job to finish, inspect the generated notes, and
confirm that all expected assets are attached.

Keep the GitHub release as a draft until package publishing succeeds.

### 4. Publish packages

Run each manual publishing workflow from the tagged release commit on `main`.
Use the available dry-run or test option first:

1. **Publish to crates.io** — select `orphos-core` and enable `dry_run`, then
   publish `orphos-core` with `dry_run` disabled. After the new core version is
   visible on crates.io, repeat those two runs for `orphos-cli`. A CLI dry run
   cannot resolve a new `orphos-core` version until the core crate has actually
   been published. For subsequent retries, publish only the crate that has not
   succeeded yet.
2. **Publish to PyPI** — optionally publish to TestPyPI first by enabling
   `test_pypi`, then run again with `test_pypi` disabled for the production
   index.
3. **Publish to npm** — enable `dry_run`, then run again with `dry_run`
   disabled to publish `@fullhuman/orphos`.

These workflows require the repository publishing credentials and trusted
publisher configuration to already be set up.

### 5. Publish and verify

Once all package workflows have succeeded, review and publish the draft GitHub
release. Verify the public artifacts:

```console
cargo search orphos-core
cargo search orphos-cli
python -m pip index versions orphos
npm view @fullhuman/orphos version
```

Perform a clean installation and a small smoke test for the interfaces changed
in the release. Finally, update the Homebrew formula and Bioconda recipe, if
applicable, using the release archive checksums, and monitor the release
workflows and issue tracker for regressions.

If a package publish fails, fix and rerun only the failed workflow. Never move
or replace a public release tag, and never reuse a version already accepted by
a package registry; prepare a new patch release instead.
