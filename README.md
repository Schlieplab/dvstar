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
* eigen3 (`eigen` in brew, `libeigen3-dev` on debian).
* cmake (`cmake` in brew, `cmake` on debian).

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
Both `kmc` and `kmc_tools` are required. [Both are available on their github releases page](https://github.com/refresh-bio/KMC/releases).

This provides an executable `dvstar`, which can be used as follows:

```shell
% dvstar --help
Construction and comparisons of variable-length Markov chains with the aid of a k-mer counter.
Usage: dvstar [OPTIONS] SUBCOMMAND

Options:
  -h,--help                   Print this help message and exit

Subcommands:
  build                       Constructs a single VLMC from the provided fasta file.
  build-from-kmc-db           Constructs a single VLMC from the provided kmc db.
  dump                        Dumps the contents of a .bintree file to raw text, either to a .txt file or stdout.
  score                       Computes the negative log-likelihood of a collection of fasta files for VLMCs.
  bic                         Runs BIC and AIC to give insight into parameter choice of the VLMC for the given fasta file.
  dissimilarity               Computes the dissimilarity between a collection of VLMCs.
  dissimilarity-fasta         Computes VLMCs from the fasta files and the dissimilarity between the resulting VLMCs.
  reprune                     Runs the pruning steps of the VLMC construction to make a given VLMC more general.
  size                        Prints the size of the VLMC. The size is determined by the sum of the length of all (terminal) k-mers in the VLMCs. A k-mer is considered terminal if any branch ends with the k-mer, meaning there is no more specific k-mer in the VLMC.
```

For example, to construct a VLMC, run:

```shell
dvstar build --threshold 3.9075 --max-depth 4 --min-count 100 --temp-path tmp NC_022098.1.fasta NC_022098.1.bintree
```

To view the contents of the VLMC, run:

```shell
dvstar dump NC_022098.1.bintree
```

To compute the dvstar similarity between two VLMCs:

```shell
dvstar dissimilarity NC_022098.1.bintree  NC_022098.1.bintree
```

To directly compute vlmcs and their dissimilarities from a directory of vlmcs:

```shell
dvstar dissimilarity-fasta  directory-with-multiple-fastas output-path
```

The output of this is distances in the phylip-format, which can then be provided to phylogenetic tools to construct trees (e.g., [`RapidNJ`](https://github.com/somme89/rapidNJ) or [`decentree`](https://github.com/iqtree/decenttree)).

To build a VLMC directly from a KMC db, ensure that the kmc parameters `-ci1` and `-cs4294967295`(or other large number)
are used. Also ensure that the `-k` parameter is set to 1 larger than the max depth parameter given to `dvstar`.

```shell
dvstar build-from-kmc-db --max-depth 4 --min-count 100 kmc_db kmc_db.bintree
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
