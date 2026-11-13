# Artifact qualification and release boundaries

The Cargo workspace package version is the single version authority. Ordinary
CPython builds enable the binding `abi3` feature targeting the consumer minimum.
Free-threaded builds omit that feature and target individual interpreters. The
binding crate is unpublished. The core crate has no Python dependency.

Candidate wheels are installed into temporary external environments for consumer
validation. Persistent consumer registry declarations and standalone locks wait
until an authorized version exists in the registry. No publishable sibling path
is substituted. The monorepo itself does not become a Cargo or uv workspace.

Build one sdist, then build each release wheel from its extracted source. Verify
installed fixtures outside the checkout, artifact contents, wheel tags, and core
source correspondence. CI covers supported glibc Linux, macOS, and Windows
runners, ordinary CPython and free-threaded interpreters. Unexecuted platform
jobs do not establish support. Musl and other interpreters are unqualified.

Repository creation, registry identities, trusted-publisher setup, tags, releases,
and uploads require separate authorization. Cargo publication repackages source;
reconstruct tested immutable source and compare repackaging before obtaining
credentials. Query registry state after upload timeouts before retrying. Complete
the full Python artifact matrix before first upload; verify current registry file
addition and identical-retry policies at release time. Conflicting artifacts stop
coordinated publication. Publish jobs must not run untrusted source after acquiring
credentials. The normal CI and release workflows mirror OpenSpeleo Core's separation of lint,
Rust tests, Python tests, platform wheel builds, and publication. Their matrices
are authoritative. Supplemental jobs retain free-threaded coverage. Every wheel
is built from the single sdist, which carries the repository toolchain file.
Ordinary builds enable the stable ABI on every interpreter in the matrix; only
the consumer-minimum row uploads each platform's wheel, avoiding duplicate
filenames. Free-threaded jobs upload interpreter-specific wheels.

The release workflow runs from a tag matching the Cargo workspace version. It
builds and verifies the core crate and checks source correspondence against the
sdist before either registry job can run. The crate publisher reconstructs and
compares the tested archive before requesting credentials. Publishing uses
`--no-verify` because package verification has already completed before OIDC
authentication.

Configure trusted publishers for the actual upstream repository and
`release.yml`: the PyPI publisher uses the `pypi` GitHub environment, and the
crates.io publisher uses `crates-io`. Only those jobs receive OIDC permission.
The Python binding crate remains unpublished. No persistent registry token is
required. Registry-side publisher configuration, environment protection, and any
initial crate bootstrap remain setup prerequisites; adding workflow files does
not configure those services.

On a partial release, inspect registry filenames and digests before retrying.
Cargo can time out after uploading a crate; check its registry checksum against
the tested archive before deciding whether publication succeeded. Never accept
conflicting artifacts as a successful retry.

The implementation is independently written under AGPL. No application solver
source or research document is redistributed. Dependency licenses remain
individually reviewable through Cargo manifests and the lockfile.
