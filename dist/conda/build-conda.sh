#!/usr/bin/env bash
# Build conda packages for the three supported platforms from prebuilt binaries,
# and optionally upload them to a personal anaconda.org channel.
#
#   ./build-conda.sh                          # all staged platforms
#   ./build-conda.sh linux-64                 # just one
#   UPLOAD=1 CHANNEL=my-lab ./build-conda.sh  # build and publish
#
# Expects each binary at staging/<platform>/metax alongside the documents that
# ship with it. stage-from-dist.sh populates that from dist/.
set -euo pipefail

HERE=$(cd "$(dirname "$0")" && pwd)
PLATFORMS=${*:-"linux-64 linux-aarch64 osx-arm64"}
OUT=${OUT:-$HERE/output}
RATTLER=${RATTLER:-rattler-build}
HOST="$(uname -s)-$(uname -m)"

# Both licence texts ship inside every package; stage-from-dist.sh copies them
# here from the repository root.
for f in LICENSE THIRD_PARTY_NOTICES.md; do
  [ -f "$HERE/$f" ] || { echo "ERROR: $HERE/$f is missing; run stage-from-dist.sh." >&2; exit 1; }
done

for p in $PLATFORMS; do
  bin="$HERE/staging/$p/metax"
  [ -f "$bin" ] || { echo "ERROR: missing $bin" >&2; exit 1; }

  # Only the package matching this host can execute its own test.
  test_flag="--no-test"
  case "$HOST:$p" in
    Linux-x86_64:linux-64|Linux-aarch64:linux-aarch64|Darwin-arm64:osx-arm64)
      test_flag="" ;;
  esac

  echo ">>> $p ${test_flag:+}${test_flag:+(cross-built: test skipped)}"
  "$RATTLER" build --recipe "$HERE/recipe.yaml" \
      --target-platform "$p" $test_flag --output-dir "$OUT"

  pkg=$(ls -t "$OUT/$p"/metax-*.conda | head -1)
  python3 "$HERE/verify-package.py" "$pkg" "$bin"
done

echo
echo "=== built packages ==="
find "$OUT" -name '*.conda' | sort

if [ "${UPLOAD:-0}" = "1" ]; then
  : "${CHANNEL:?set CHANNEL to your anaconda.org user or organisation}"
  echo
  echo ">>> uploading to anaconda.org/$CHANNEL"
  # Requires ANACONDA_API_TOKEN in the environment, or a prior `anaconda login`.
  find "$OUT" -name '*.conda' -exec anaconda upload --user "$CHANNEL" --label main {} +
fi
