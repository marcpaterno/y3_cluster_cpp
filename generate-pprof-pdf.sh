#!/usr/bin/env bash

pushd /cosmosis/cosmosis-standard-library/y3_cluster_cpp
google-pprof --pdf test/trivial_gamma_t dump.txt > ${1:-google-dump.pdf}
popd
