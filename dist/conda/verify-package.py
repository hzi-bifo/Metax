#!/usr/bin/env python3
"""Check that a built conda package carries the exact binary that was staged.

Any rewrite by the packaging step would break the ad-hoc code signature that
arm64 macOS requires at exec time, so this runs on every build.

usage: verify-package.py <package.conda> <staged-binary>
"""
import hashlib
import io
import sys
import tarfile
import zipfile

import zstandard


def packaged_binary(pkg: str) -> bytes:
    with zipfile.ZipFile(pkg) as z:
        member = next(n for n in z.namelist() if n.startswith("pkg-"))
        blob = z.read(member)
    stream = zstandard.ZstdDecompressor().stream_reader(io.BytesIO(blob))
    with tarfile.open(fileobj=stream, mode="r|") as tar:
        for entry in tar:
            if entry.name.lstrip("./") == "bin/metax":
                handle = tar.extractfile(entry)
                if handle is None:
                    raise SystemExit("bin/metax is not a regular file")
                return handle.read()
    raise SystemExit("bin/metax not found in package")


def main() -> int:
    pkg, staged = sys.argv[1], sys.argv[2]
    got = hashlib.sha256(packaged_binary(pkg)).hexdigest()
    with open(staged, "rb") as fh:
        want = hashlib.sha256(fh.read()).hexdigest()
    if got != want:
        print(f"FAIL packaging altered the binary\n  staged   {want}\n  packaged {got}")
        return 1
    print(f"  verified identical: {got[:16]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
