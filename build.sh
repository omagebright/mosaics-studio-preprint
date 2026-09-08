#!/bin/sh
# Build both documents with tectonic. Fails on the first error.
set -e
cd "$(dirname "$0")"
python3 make_tables.py > /dev/null
tectonic -X compile main.tex
tectonic -X compile si.tex
