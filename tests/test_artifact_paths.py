"""Archive member names must not depend on the host filesystem's path rules."""

import io
import tarfile
from pathlib import PurePosixPath, PureWindowsPath

import inspect_artifacts
import pytest


@pytest.mark.parametrize("host_path", [PurePosixPath, PureWindowsPath])
def test_tar_members_keep_posix_names_and_exact_bytes(tmp_path, monkeypatch, host_path):
    archive_path = tmp_path / "source.tar.gz"
    expected = {
        "tests/fixtures/controlled_tree.json": b'{"schema_version": 1}\n',
        "crates/cave_closure/src/solve.rs": b"// source\r\n",
        "LICENSE": b"license\n",
    }
    with tarfile.open(archive_path, "w:gz") as archive:
        directory = tarfile.TarInfo("package/tests/fixtures")
        directory.type = tarfile.DIRTYPE
        archive.addfile(directory)
        for name, content in expected.items():
            member = tarfile.TarInfo(f"package/{name}")
            member.size = len(content)
            archive.addfile(member, io.BytesIO(content))

    # Reproduce Windows Path semantics without needing a Windows filesystem.
    # Reading the archive itself still uses tmp_path's real filesystem path.
    monkeypatch.setattr(inspect_artifacts, "Path", host_path)
    assert inspect_artifacts.members(archive_path) == expected
