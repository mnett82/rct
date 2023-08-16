#!/bin/bash
set eux -o pipefail
mkdir -p ./data/
wget http://ann-benchmarks.com/mnist-784-euclidean.hdf5 -P ./data/