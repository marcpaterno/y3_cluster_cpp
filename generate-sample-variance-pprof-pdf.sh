#!/usr/bin/env bash

pushd /cosmosis/cosmosis-standard-library/y3_cluster_cpp
google-pprof --pdf test/sample_variance test/sample_variance_output/dump.txt > ${1:-sample-variance-pprof.pdf}
popd
