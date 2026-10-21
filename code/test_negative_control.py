#!/usr/bin/env python3
"""Hostile semantic check for the exact finite SIS verifier.

The test creates a temporary copy requiring every jet determinant to vanish.
Both ordinary and optimized Python must reject it.  The reviewed source file is
never modified.
"""
from pathlib import Path
import subprocess
import sys
import tempfile


SOURCE = Path(__file__).with_name("verify_sis.py")
ORIGINAL = "require(det==expected and det!=0, f'jet determinant identity failed at d={d}')"
MUTATED = "require(det==0, f'hostile control: zero jet determinant required at d={d}')"


def rejected(command, script):
    completed = subprocess.run(
        [sys.executable, *command, str(script)],
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    return completed.returncode != 0 and "VerificationError" in completed.stderr


def main():
    source = SOURCE.read_text(encoding="utf-8")
    if source.count(ORIGINAL) != 1:
        raise RuntimeError("expected exactly one determinant obligation to mutate")
    hostile = source.replace(ORIGINAL, MUTATED)
    with tempfile.TemporaryDirectory(prefix="sis-negative-control-") as tmp:
        script = Path(tmp) / SOURCE.name
        script.write_text(hostile, encoding="utf-8")
        if not rejected([], script):
            raise RuntimeError("ordinary Python accepted the corrupted determinant claim")
        if not rejected(["-O"], script):
            raise RuntimeError("optimized Python accepted the corrupted determinant claim")
    print("PASS: ordinary and optimized Python reject the corrupted determinant claim")


if __name__ == "__main__":
    main()
