#!/usr/bin/env python3
"""Generate the deterministic SHA-256 manifest for source and documentation."""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path
from typing import Iterable, Sequence


ROOT_FILES = (
    ".github/workflows/ci.yml",
    ".gitignore",
    ".zenodo.json",
    "CHANGELOG.md",
    "CITATION.cff",
    "FiniteNativeCarryOperator.lean",
    "LICENSE",
    "README.md",
    "lake-manifest.json",
    "lakefile.toml",
    "lean-toolchain",
    "requirements.txt",
)
SOURCE_DIRECTORIES = (
    "FiniteNativeCarryOperator",
    "certification",
    "docs",
    "laboratory",
    "scripts",
    "tests",
)
EXCLUDED_SOURCE_DIRECTORIES = (
    Path("docs/notes"),
)


def source_paths(root: Path) -> Iterable[Path]:
    for relative in ROOT_FILES:
        path = root / relative
        if path.is_file():
            yield path
    for directory in SOURCE_DIRECTORIES:
        base = root / directory
        if not base.exists():
            continue
        for path in base.rglob("*"):
            relative = path.relative_to(root)
            if (
                path.is_file()
                and "__pycache__" not in path.parts
                and not any(relative.is_relative_to(excluded) for excluded in EXCLUDED_SOURCE_DIRECTORIES)
            ):
                yield path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def build_manifest(root: Path) -> str:
    unique = sorted(set(source_paths(root)), key=lambda path: path.relative_to(root).as_posix())
    return "".join(
        f"{sha256(path)}  {path.relative_to(root).as_posix()}\n" for path in unique
    )


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("audit/SOURCE_SHA256.txt"))
    args = parser.parse_args(argv)
    root = Path(__file__).resolve().parents[1]
    output = args.output if args.output.is_absolute() else root / args.output
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(build_manifest(root), encoding="utf-8")
    print(f"source manifest written to {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
