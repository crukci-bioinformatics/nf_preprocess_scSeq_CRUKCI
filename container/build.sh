#!/bin/bash

set -ex

# Version is the version of CellRanger in the container.
version=8.0.1

sudo docker build --tag cruk_ci_preprocess_scSeq:${version} .
singularity build cruk_ci_preprocess_scSeq-${version}.sif docker-daemon://cruk_ci_preprocess_scSeq:${version}
