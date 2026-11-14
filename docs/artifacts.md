# Artifact qualification and release boundaries

The Cargo package version is the single version authority. Ordinary CPython
builds enable the `abi3` feature targeting the consumer minimum. Free-threaded
builds use `abi3t` where supported and an interpreter-specific legacy wheel
below its ABI floor. Every Python wheel explicitly enables both `sparse-direct`
and `l1`. Python bindings are enabled by the `python` feature. Rust consumers
can use the solver without enabling Python.

Candidate wheels are installed into temporary external environments for consumer
validation. Persistent consumer registry declarations and standalone locks wait
until an authorized version exists in the registry. No publishable sibling path
is substituted. The monorepo itself does not become a Cargo or uv workspace.

CI Rust and Python tests and release wheel builds use the same native platform
matrix. CI covers every configured ordinary and free-threaded Python
interpreter. Release builds the ordinary stable-ABI and two free-threaded wheel
variants on each platform:

| Platform       | Runner             |
| -------------- | ------------------ |
| Linux x86-64   | `ubuntu-latest`    |
| Linux ARM64    | `ubuntu-24.04-arm` |
| Windows x86-64 | `windows-latest`   |
| Windows ARM64  | `windows-11-arm`   |
| macOS ARM64    | `macos-latest`     |

Builds run directly on the runners without containers or cross-compilation.
Python test coverage is authoritative in the CI matrix; release wheel build
interpreters and ABI variants are authoritative in the release matrix. Python
tests run in CI, not in release wheel jobs. The Rust source crate is packaged,
verified, and uploaded once on `ubuntu-latest`; crate packaging has no platform
matrix.

Linux wheel compatibility follows the runner's glibc and resulting wheel tags;
these builds do not promise older manylinux compatibility. Musl, ARMv7, PowerPC,
and Windows x86 wheels are not built. Unexecuted platform jobs do not establish
support. HiGHS compilation requires CMake, a C++ compiler, and a discoverable
libclang shared library on the native runner.

Repository creation, registry identities, trusted-publisher setup, tags,
releases, and uploads require separate authorization. Cargo publication
repackages source; reconstruct tested immutable source and compare repackaging
before obtaining credentials. Query registry state after upload timeouts before
retrying. Complete the full Python artifact matrix before first upload; verify
current registry file addition and identical-retry policies at release time.
Conflicting artifacts stop coordinated publication. Publish jobs must not run
untrusted source after acquiring credentials. The existing CI and release
workflows retain separate jobs for lint, Rust tests, Python tests, artifact
builds, and publication. Their matrices are authoritative. Release builds one
ordinary stable-ABI wheel per platform at the consumer-minimum interpreter; the
ordinary build matrix has only the platform dimension. The existing CI gate runs
the Python tests before automatically dispatching a tagged release.
Free-threaded jobs upload the stable-ABI or legacy interpreter-specific wheel as
appropriate. The sdist carries the repository toolchain file, fixtures, schemas,
and documentation; artifact inspection checks crate/sdist core source
correspondence and wheel contents before publication.

The release workflow runs from a tag matching the Cargo package version. It
builds and verifies the core crate and checks source correspondence against the
sdist before either registry job can run. The crate publisher reconstructs and
compares the tested archive before requesting credentials. Publishing uses
`--no-verify` because package verification has already completed before OIDC
authentication.

Configure trusted publishers for the actual upstream repository and
`release.yml`: the PyPI publisher uses the `pypi` GitHub environment, and the
crates.io publisher uses `crates-io`. Only those jobs receive OIDC permission.
The root Cargo package supplies both the Rust library and optional Python
bindings. No persistent registry token is required. Registry-side publisher
configuration, environment protection, and any initial crate bootstrap remain
setup prerequisites; adding workflow files does not configure those services.

On a partial release, inspect registry filenames and digests before retrying.
Cargo can time out after uploading a crate; check its registry checksum against
the tested archive before deciding whether publication succeeded. Never accept
conflicting artifacts as a successful retry.

The implementation is independently written under AGPL. No application solver
source or research document is redistributed. Dependency licenses remain
individually reviewable through Cargo manifests and the lockfile.
