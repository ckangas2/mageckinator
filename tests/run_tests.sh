#!/usr/bin/env bash
# Unit tests: one planted failure per gate, plus the passing case.
# Seconds; no MAGeCK, Docker, R or network.
set -euo pipefail
cd "$(dirname "$0")/.."
python3 -m pytest tests -q "$@"
