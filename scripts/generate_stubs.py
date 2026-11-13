"""Generate formatted Python stubs, or check them without modifying sources."""

import argparse
import difflib
import os
import subprocess
import sys
import sysconfig
import tempfile
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="fail if stubs differ")
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    env = os.environ.copy()
    env["PYO3_PYTHON"] = sys.executable
    if sys.platform == "linux":
        library_dir = sysconfig.get_config_var("LIBDIR")
        env["LD_LIBRARY_PATH"] = os.pathsep.join(
            filter(None, [library_dir, env.get("LD_LIBRARY_PATH")])
        )

    with tempfile.TemporaryDirectory(prefix="cave-closure-stubs-") as directory:
        subprocess.run(
            [
                "cargo",
                "run",
                "--locked",
                "--features",
                "generate-stubs",
                "--bin",
                "stub_gen",
                "--",
                "--output-dir",
                directory,
            ],
            cwd=root,
            env=env,
            check=True,
        )
        generated = Path(directory) / "cave_closure/_native/__init__.pyi"
        for command in (["check", "--fix"], ["format"]):
            subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "ruff",
                    *command,
                    "--config",
                    str(root / "pyproject.toml"),
                    str(generated),
                ],
                cwd=root,
                check=True,
            )
        target = root / "src_python/cave_closure/_native/__init__.pyi"
        content = generated.read_text()
        if not args.check:
            target.write_text(content)
            return 0
        current = target.read_text() if target.exists() else ""
        if current == content:
            print("Generated Python stubs are up to date.")
            return 0
        sys.stdout.writelines(
            difflib.unified_diff(
                current.splitlines(keepends=True),
                content.splitlines(keepends=True),
                fromfile=str(target.relative_to(root)),
                tofile="regenerated stubs",
            )
        )
        print("Run: uv run --frozen python scripts/generate_stubs.py")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
