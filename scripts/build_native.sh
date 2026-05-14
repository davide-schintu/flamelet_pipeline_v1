#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PYTHON_BIN="${PYTHON_BIN:-$(command -v python)}"

echo "[native] Python: ${PYTHON_BIN}"

echo "[native] building RFThermo"
cd "${ROOT_DIR}/flamelet_pipeline/rfthermo"
"${PYTHON_BIN}" setup.py build_ext --inplace
mkdir -p build
cp RFThermo*.so build/

echo "[native] building PDF convolution module"
cd "${ROOT_DIR}/flamelet_pipeline/cpp_pdf"
make python PYTHON_BIN="${PYTHON_BIN}"

echo "[native] done"
