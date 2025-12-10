# Fast construction of variable-length Markov chains (VLMC) from k-mers

This repository contains the code for constructing a variable-length Markov chain (VLMC) from a lexicographically
sorted list of _k_-mers. This enables VLMCs to be constructed very efficiently even for very large genomes, as well as collections
of sequencing data.

This repository deals mainly with the task of constructing the VLMC, for scoring of sequences, or most dissimilarities
between
two VLMCs, please see the [pst-classifier-seqan repository](https://github.com/Schlieplab/PstClassifierSeqan), which
accepts the output form this progam.

## Publication

Publication pending...

## Installation

### Installing KMC3
__Note that KMC3 needs to be installed separately, and both executables `kmc` and `kmc_tools` are required to be in the current working directory, or in your path.__ [`kmc` and `kmc_tools are available on their github releases page of the KMC github](https://github.com/refresh-bio/KMC/releases).


### Container solution / apptainer (previously singularity)

We provide an [apptainer container](https://apptainer.org/). After downloading the `vlmc-from-kmers.sif` `vlmc-from-kmers.def`, the definition file `vlmc-from-kmers.def` needs to be modified to
use the correct path to kmc, see the fifth and sixth line, specifically
```shell script
    ../kmc/KMC3.linux/kmc /kmc/kmc
    ../kmc/KMC3.linux/kmc_tools /kmc/kmc_tools
```

The container is built by running:

```shell script
sudo apptainer build vlmc-from-kmers.sif vlmc-from-kmers.def
```

### Using Spack 
We provide an initial `package.py` to be used with the Spack package manager. Activate the separately installed KMC3 package (and possibly other packages
from the dependencies below)

### Compiling from Source (recommended)

To compile the program, please first make sure that in addition to KMC3 the following dependencies are installed.

#### Dependencies of dvstar 
* Intel TBB (`tbb` in brew, `libtbb-dev` on debian).
* Boost iostreams (`boost` in brew, `libboost-dev` on debian).
* hdf5 (`hdf5` in brew, `libhdf5-dev` on debian).
* eigen3 (`eigen` in brew, `libeigen3-dev` on debian).
* cmake (`cmake` in brew, `cmake` on debian).

This can be accomplished on a Linux system with package manager `apt` with the command  
```shell script
sudo apt install libtbb-dev libboost-dev libhdf5-dev libeigen3-dev cmake
```
Make sure that `apt` and your build tools are up-to-date, e.g. by
```shell script
sudo apt update
sudo apt install build-essential
```

#### Clone the dvstar repository including submodules

```shell script
git clone https://github.com/Schlieplab/dvstar.git
cd dvstar
git submodule update --init --recursive
```

#### Build
Create a build directory, configure with cmake and build:

```shell script
cmake -DCMAKE_BUILD_TYPE=Release -S . -B build
cmake --build build
```

Install dvstar to you prefered location by using 
```shell script
sudo cp build/dvstar /usr/local/bin/
sudo chmod +x /usr/local/bin/dvstar
```


#### Build On MacOS
You may need to install a different c++ compiler than the one Apple provides. This can, for example, be done through [brew](https://brew.sh/) or [MacPorts](https://www.macports.org/). This is tested with gcc (version 13) so that would be my recommendation. You may also need to provide the path to the installed compiler in the cmake command, e.g.:

```shell script
CC=/opt/homebrew/bin/gcc-13 CXX=/opt/homebrew/bin/g++-13 cmake -DCMAKE_BUILD_TYPE=Release -S . -B build
cmake --build build
```


# Running dvstar
After successful compilation you can run the dvstar executable.

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
