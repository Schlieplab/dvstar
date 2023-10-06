#!/bin/bash

# apptainer can't handle the circular submodule in CLI11 tests
rm -rf submodules/CLI11/tests/mesonTest/subprojects/CLI11

mkdir -p artifacts

sudo apptainer build --force vlmc-from-kmers.sif vlmc-from-kmers.def

apptainer exec --bind $PWD/artifacts/:/artifacts vlmc-from-kmers.sif sh -c 'cp /vlmc-from-kmers/build/dvstar /artifacts'
