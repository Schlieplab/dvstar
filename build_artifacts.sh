#!/bin/bash

sudo apptainer build --force vlmc-from-kmers.sif vlmc-from-kmers.def

apptainer exec --bind $PWD/artifacts/:/artifacts vlmc-from-kmers.sif sh -c 'cp /vlmc-from-kmers/build/dvstar /artifacts'
