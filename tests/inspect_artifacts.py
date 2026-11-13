"""Check packaged core correspondence and emit a machine-readable artifact record."""

import argparse
import hashlib
import json
import tarfile
import zipfile
from pathlib import Path, PurePosixPath


def members(path):
    with tarfile.open(path) as archive:
        return {
            PurePosixPath(member.name)
            .relative_to(PurePosixPath(member.name).parts[0])
            .as_posix(): archive.extractfile(member).read()
            for member in archive.getmembers()
            if member.isfile()
        }


def inspect(crate, sdist, wheel):
    rust, source = members(crate), members(sdist)
    for path in (
        "rust-toolchain.toml",
        "tests/fixtures/controlled_tree.json",
        "schemas/input.json",
        "schemas/result.json",
        "docs/numerics.md",
        "docs/interchange.md",
        "docs/artifacts.md",
    ):
        if path not in source:
            raise ValueError(
                f"sdist is missing required corpus/schema/documentation: {path}"
            )
    corpus_path = "tests/fixtures/corpus.json"
    if corpus_path not in source:
        raise ValueError("sdist is missing the persisted fixture manifest")
    for case in json.loads(source[corpus_path])["cases"]:
        for filename in (case["input"], case["reference"], *case["results"].values()):
            if "tests/fixtures/" + filename not in source:
                raise ValueError(f"sdist is missing a corpus artifact: {filename}")
    for path in (
        "LICENSE",
        "src_rust/lib.rs",
        "src_rust/model.rs",
        "src_rust/prepare.rs",
        "src_rust/solve.rs",
        "src_rust/json.rs",
        "src_rust/normalize.rs",
        "tests/core.rs",
        "tests/fixtures/controlled_tree.json",
        "examples/fixture_results.rs",
    ):
        if rust[path] != source[path]:
            raise ValueError(f"crate/sdist source mismatch: {path}")
    with zipfile.ZipFile(wheel) as archive:
        files = set(archive.namelist())
        required = {
            "cave_closure/__init__.py",
            "cave_closure/__init__.pyi",
            "cave_closure/_native/__init__.pyi",
            "cave_closure/py.typed",
        }
        if not required.issubset(files):
            raise ValueError("wheel is missing Python interface files")
        if not any(
            name.startswith("cave_closure/_native") and name.endswith((".so", ".pyd"))
            for name in files
        ):
            raise ValueError("wheel is missing its compiled extension")
        if not any(name.endswith("/licenses/LICENSE") for name in files):
            raise ValueError("wheel is missing its license")
        wheel_metadata = next(
            name for name in files if name.endswith(".dist-info/WHEEL")
        )
        tags = [
            line.removeprefix("Tag: ")
            for line in archive.read(wheel_metadata).decode().splitlines()
            if line.startswith("Tag: ")
        ]
    return {
        "artifacts": [
            {"filename": p.name, "sha256": hashlib.sha256(p.read_bytes()).hexdigest()}
            for p in (crate, sdist, wheel)
        ],
        "core_source_correspondence": "verified",
        "wheel_tags": tags,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--crate", type=Path, required=True)
    parser.add_argument("--sdist", type=Path, required=True)
    parser.add_argument("--wheel", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.write_text(
        json.dumps(inspect(args.crate, args.sdist, args.wheel), indent=2)
    )
    print("Artifact contents, core source correspondence, and wheel tags verified")


if __name__ == "__main__":
    main()
