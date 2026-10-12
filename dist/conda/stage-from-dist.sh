#!/usr/bin/env bash
# Copy release binaries into staging/<conda-platform>/ and refresh the
# third-party notices that ship in the package metadata.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
DIST=${DIST:-$(cd "$HERE/.." && pwd)}
VERSION=${VERSION:-2.0.3}
REPO=${REPO:-$(cd "$HERE/../.." && pwd)}

map_linux_64=x86_64-unknown-linux-gnu
map_linux_aarch64=aarch64-unknown-linux-gnu
map_osx_arm64=aarch64-apple-darwin

cp "$REPO/THIRD_PARTY_NOTICES.md" "$HERE/THIRD_PARTY_NOTICES.md"
cp "$REPO/LICENSE" "$HERE/LICENSE"

for p in linux-64 linux-aarch64 osx-arm64; do
  var="map_${p//-/_}"; triple="${!var}"
  src="$DIST/metax-v$VERSION-$triple/metax"
  [ -f "$src" ] || { echo "skip $p (no $src)"; continue; }
  mkdir -p "$HERE/staging/$p"
  cp "$src" "$HERE/staging/$p/metax"
  chmod 755 "$HERE/staging/$p/metax"
  echo "staged $p from $triple"
done
