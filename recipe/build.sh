#!/usr/bin/env bash
set -euxo pipefail

mkdir -p "${PREFIX}/bin"
cp "${SRC_DIR}/metax" "${PREFIX}/bin/metax"
chmod +x "${PREFIX}/bin/metax"
