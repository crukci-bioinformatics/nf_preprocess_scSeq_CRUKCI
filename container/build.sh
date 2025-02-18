#!/bin/bash

set -ex

# Version is the version of CellRanger in the container.
version=9.0.1

sudo docker build --tag cruk_ci_preprocess_scseq:${version} .
sudo singularity build cruk_ci_preprocess_scseq-${version}.sif docker-daemon://cruk_ci_preprocess_scseq:${version}
