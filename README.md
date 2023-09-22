# Building and comparing a variable-length Markov chain from k-mers

This repository contains the code for constructing a variable-length Markov chain (VLMC) from a lexicographically
sorted list of _k_-mers. This enables VLMCs to be constructed on even very large genomes, as well as collections
of sequences.

This repository deals mainly with the task of constructing the VLMC, for scoring of sequences, or most dissimilarities
between
two VLMCs, please see the [pst-classifier-seqan repository](https://github.com/Schlieplab/PstClassifierSeqan), which
accepts the output form this program.

## Publication

Publication pending...

## Compilation

To compile the program, please first download all submodules:

```shell script
git submodule update --init --recursive
```

### Container solution / apptainer (previously singularity)

We provide an [apptainer container](https://apptainer.org/). The definition file needs to be modified to
include the path to kmc, see the fifth and sixth line. The container is built by running:

```shell script
sudo apptainer build vlmc-from-kmers.sif vlmc-from-kmers.def
```

### Manually

#### Dependencies:

* intel TBB (`tbb` in brew, `libtbb-dev` on debian).
* boost iostreams (`boost` in brew, `libboost-dev` on debian).
* hdf5 (`hdf5` in brew, `libhdf5-dev` on debian).
* eigen3 (`eigen` in brew, `libeigen3-dev` on debian.)

Create a build directory, configure with cmake and build:

```shell script
cmake -DCMAKE_BUILD_TYPE=Release -S . -B build
cmake --build build
```

#### On MacOS
You may need to install a different c++ compiler than the one apple provides. This can, for example, be done through [brew](https://brew.sh/) or [MacPorts](https://www.macports.org/). This is tested with gcc (version 13) so that would be my recommendation. You may also need to provide the path to the installed compiler in the cmake command, e.g.:

```shell script
CC=/opt/homebrew/bin/gcc-13 CXX=/opt/homebrew/bin/g++-13 cmake -DCMAKE_BUILD_TYPE=Release -S . -B build
cmake --build build
```

__Note that kmc3 needs to be installed separately, and be in the current working directory, or in your path.__
Both `kmc` and `kmc_tools` are required.

This provides an executable `dvstar`, which can be used as follows:

```shell
% ./dvstar --help
Construction and comparisons of variable-length Markov chains with the aid of a k-mer counter.
Usage: ./dvstar [OPTIONS]

Options:
  -h,--help                   Print this help message and exit
  -m,--mode ENUM:value in {bic->3,build->0,build-from-kmc-db->4,dissimilarity->5,dump->2,reprune->6,score->1} OR {3,0,4,5,2,6,1}
                              Program mode, 'build', 'build-from-kmc-db', 'dump', 'score', 'reprune', or 'dissimilarity'.  For build-from-kmc-db, the kmc db needs to include all k-mers (not in canonical form), with no minimum count cutoff.  The length of the k-mers needs to be set to 1 more than the maximum depth of the VLMC.  The kmc db also has to be sorted.
  --dissimilarity ENUM:value in {dvstar->0,penalized-dvstar->1} OR {0,1}
                              Dissimilarity type, either 'dvstar',  or 'penalized-dvstar'.
  --estimator ENUM:value in {kullback-leibler->0,peres-shields->1} OR {0,1}
                              Estimator for the pruning of the VLMC, either 'kullback-leibler',  or 'peres-shields'.
  -p,--fasta-path TEXT        Path to fasta file.  Required for 'build' and 'score' modes.
  --in-path TEXT              Path to saved tree file or kmc db file.  Required for 'build-from-kmc-db', 'dump', 'score', and 'dissimilarity' modes.  For 'build-from-kmc-db', the kmc db file needs to be supplied without the file extension.
  --to-path TEXT              Path to saved tree file.  Required for 'dissimilarity' mode.
  -o,--out-path TEXT          Path to output file.  The VLMCs are stored as binary, and can be read by the 'dump' or 'score' modes.  Required for 'build' and 'dump' modes.
  -t,--temp-path TEXT         Path to temporary folder for the external memory algorithms.  For good performance, this needs to be on a local machine.  For sorting, at least 2GB will be allocated to this path.  Defaults to ./tmp
  -c,--min-count INT          Minimum count required for every k-mer in the tree.
  -k,--threshold FLOAT        Kullback-Leibler threshold.
  -d,--max-depth INT          Maximum depth/length for included k-mers.
  -a,--pseudo-count-amount FLOAT
                              Size of pseudo count for probability estimation. See e.g. https://en.wikipedia.org/wiki/Additive_smoothing .
  -i,--in-or-out-of-core ENUM:value in {external->0,internal->1,hash->2} OR {0,1,2}
                              Specify 'internal' for in-core or 'external' for out-of-core memory model.  Out of core is slower, but is not memory bound.
  --adjust-for-sequencing-errors
                              Give this flag to adjust the estimator parameters and min counts for the sequencing depth and error rates of a sequencing dataset. See --sequencing-depth and --sequencing-error-rate for parameters.
  --sequencing-depth FLOAT    If --adjust-for-sequencing-errors is given, this parameter is used to alter the estimator parameters to reflect that many k-mers will be --sequencing-depth times more frequent.
  --sequencing-error-rate FLOAT
                              If --adjust-for-sequencing-errors is given, this parameter is used to alter to estimate the number of k-mers that will be missing due to sequencing errors.
```

For example, to construct a VLMC, run:

```shell
./dvstar --fasta-path NC_022098.1.fasta --threshold 3.9075 --max-depth 4 --min-count 100 --out-path NC_022098.1.bintree --temp-path tmp
```

To view the contents of the VLMC, run:

```shell
./dvstar --mode dump --in-path NC_022098.1.bintree
```

To compute the dvstar similarity between two VLMCs:

```shell
./dvstar --mode dissimilarity --dissimilarity dvstar --in-path NC_022098.1.bintree --to-path NC_022098.1.bintree
```

To directly compute vlmcs and their dissimilarities from a directory of vlmcs:

```shell
./dvstar --mode dissimilarity-fasta --fasta-path directory-with-multiple-fastas --out-path output-path
```

The output of this is distances in the phylip-format, which can then be provided to phylogenetic tools to construct trees (e.g., [`RapidNJ`](https://github.com/somme89/rapidNJ) or [`decentree`](https://github.com/iqtree/decenttree)).

To build a VLMC directly from a KMC db, ensure that the kmc parameters `-ci1` and `-cs4294967295`(or other large number)
are used. Also ensure that the `-k` parameter is set to 1 larger than the max depth parameter given to `dvstar`.

```shell
./dvstar --mode build-from-kmc-db --in-path kmc_db --out-path kmc_db.bintree --max-depth 4 --min-count 100
```

This approach also allows you to add new sequences/reads to an existing kmc db and then rerun the vlmc construction.
This can save computation time, especially when dealing with large sequence collections. The kmc command to run to add
the result of two kmc dbs would be:

```shell
./kmc_tools simple kmc_db1 kmc_db2 union new_db_path -ocsum
```

## Comments on parameter selection

There are three parameters to take into account, "threshold", "max-depth", and "min-count" when building VLMCs with this method.

Roughly, increasing the "threshold" gives a more general VLMC (and vice versa). When performing parameter selection using the BIC or similar methods, we've seen that the "threshold" parameter increases with sequence length, such that more pruning is better with longer genomes. However, this also leads to quite general models, which are not suitable for all applications. In practice, we've had success with either a threshold of 3.9075, 1.2, or 0.1.

We've had some success with building neighbor-joining trees based on the VLMCs and the dvstar dissimilarity with a "max-depth" of 9 and "min-count" of 10, but I would probably tune these parameters depending on how sensitive the VLMCs need to be. Increasing "max-depth" and decreasing "min-count" gives a more specific VLMC.

Note also the "pseudo-count-amount" parameter, which by default adds 1 to each next-symbol count (present or not). In deep branches, this can introduce more noise than might be preferable, and thus, setting the "pseudo-count-amount" to 0.01 might give better results.

## Visualisation

To produce a graphical representation of the VLMCs, see the [`visualiser.py`](visualiser.py) script.
It needs the `plotly` and `typer` python packages to run, which can be installed e.g.
through `pip install --user plotly typer`.

```shell
python visualiser.py NC_022098.1.bintree
```

Produces an interactive visualisation (should open in your browser) as well as a pdf and png file in the `images`
directory.

## Headers

If, for some reason, you wanted to include the code in some other project, this directory can be included with CMAKE as
`add_subdirectory(/path/to/vlmc-from-kmers)`, or alternatively, include the headers
in `/path/to/vlmc-from-kmers/include`.
Note that with the second option, the submodules need to be included manually.
