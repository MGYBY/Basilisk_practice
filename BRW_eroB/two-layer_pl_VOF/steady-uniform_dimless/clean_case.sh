#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
cd "$ROOT"
rm -f steady_uniform_relaxation perfs
rm -rf profiles history analysis __pycache__
rm -f verification/compile.stdout verification/compile.stderr \
      verification/run.stdout verification/run.stderr \
      verification/runtime.txt verification/postprocess.log \
      verification/automated_checks.txt
