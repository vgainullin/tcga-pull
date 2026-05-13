"""Installer for the official NCI gdc-client binary.

NCI publishes prebuilt zips at gdc.cancer.gov. They're zip-in-a-zip and have no
ARM Linux build. MD5s are pinned per-platform — if NCI re-publishes v2.3, the
install will refuse rather than silently install a different binary.
"""

from __future__ import annotations

import hashlib
import platform
import shutil
import subprocess
import sys
import tempfile
import zipfile
from dataclasses import dataclass
from pathlib import Path

import requests

GDC_CLIENT_VERSION = "2.3"
_BASE = "https://gdc.cancer.gov/system/files/public/file"


@dataclass(frozen=True)
class _Build:
    url: str
    md5: str
    arch_note: str


# Keyed by (sys.platform, machine). The macOS-14 build is published with `x64`
# in the filename but is actually arm64-native — verified with `file(1)`.
_BUILDS: dict[tuple[str, str], _Build] = {
    ("darwin", "arm64"): _Build(
        f"{_BASE}/gdc-client_2.3_OSX_x64-py3.8-macos-14.zip",
        "56cca3594fa5fb47bc8297f5b6fd0e20",
        "macOS arm64 (Apple Silicon)",
    ),
    ("darwin", "x86_64"): _Build(
        f"{_BASE}/gdc-client_2.3_OSX_x64-py3.8-macos-12.zip",
        "fee6a557d16a6c1a9388bd859224e638",
        "macOS x86_64 (Intel)",
    ),
    ("linux", "x86_64"): _Build(
        f"{_BASE}/gdc-client_2.3_Ubuntu_x64-py3.8-ubuntu-20.04.zip",
        "18591d74de07cdcd396dab71c52663da",
        "Linux x86_64",
    ),
    ("win32", "AMD64"): _Build(
        f"{_BASE}/gdc-client_2.3_Windows_x64-py3.8-windows-2019.zip",
        "525ce44bb5f3f0624066b906c7dbdaf4",
        "Windows x64",
    ),
}


class InstallError(RuntimeError):
    """User-actionable install failure."""


def _detect_build() -> _Build:
    key = (sys.platform, platform.machine())
    build = _BUILDS.get(key)
    if build is None:
        raise InstallError(
            f"No prebuilt gdc-client available for {sys.platform}/{platform.machine()}. "
            f"NCI publishes binaries only for: "
            f"{', '.join(b.arch_note for b in _BUILDS.values())}. "
            f"Build from source: https://github.com/NCI-GDC/gdc-client"
        )
    return build


def _md5(path: Path) -> str:
    h = hashlib.md5()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def install_gdc_client(dest_dir: Path, *, force: bool = False) -> Path:
    """Download, verify, and install gdc-client. Returns the installed binary path."""
    build = _detect_build()
    binary_name = "gdc-client.exe" if sys.platform == "win32" else "gdc-client"
    dest = dest_dir / binary_name

    if dest.exists() and not force:
        raise InstallError(f"{dest} already exists. Pass force=True (or --force) to overwrite.")

    dest_dir.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory() as tmp_str:
        tmp = Path(tmp_str)
        outer_zip = tmp / "outer.zip"

        with requests.get(build.url, stream=True, timeout=60) as r:
            r.raise_for_status()
            with outer_zip.open("wb") as fh:
                for chunk in r.iter_content(chunk_size=1 << 20):
                    fh.write(chunk)

        actual = _md5(outer_zip)
        if actual != build.md5:
            raise InstallError(
                f"MD5 mismatch: expected {build.md5}, got {actual}. "
                f"NCI may have re-published v{GDC_CLIENT_VERSION}. "
                f"Refusing to install."
            )

        # NCI ships a zip wrapping a zip wrapping the binary.
        with zipfile.ZipFile(outer_zip) as zf:
            inner_name = zf.namelist()[0]
            zf.extract(inner_name, tmp)
        inner_zip = tmp / inner_name

        with zipfile.ZipFile(inner_zip) as zf:
            bin_name_in_zip = zf.namelist()[0]
            zf.extract(bin_name_in_zip, tmp)
        extracted = tmp / bin_name_in_zip

        shutil.move(str(extracted), dest)

    if sys.platform != "win32":
        dest.chmod(0o755)

    # Verify it runs.
    try:
        out = subprocess.run(
            [str(dest), "--version"], capture_output=True, text=True, timeout=10, check=True
        )
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired, OSError) as e:
        raise InstallError(f"Installed binary at {dest} but it failed to run: {e}") from e
    if GDC_CLIENT_VERSION not in (out.stdout + out.stderr):
        raise InstallError(
            f"Installed binary reports unexpected version: {out.stdout!r} {out.stderr!r}"
        )

    return dest
