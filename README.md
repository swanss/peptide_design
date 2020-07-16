# peptide_design

This repository incorporates several tools that enable the de novo design of peptides to bind a given protein. The process begins with the generation of **seeds**, which are small segments of protein backbone arranged in the space around the selected binding site. Given these seeds, the programs enable the construction of a **seed graph** describing the overlaps between these seeds at various RMSD thresholds, the scoring of these seeds by various metrics, and the construction of **paths** that thread through the graph and optimize the score function.

## Build Instructions

Before building, adjust the `makefile` variable `MSTDIR` to be the path at which MST is stored on your system. The default location is in the parent directory of `peptide_design`.

* `make all` - builds all programs
* `make test` - builds only programs in the `tests` directory
* `make bin/[executable name]` - builds the specific executable with its dependencies
* `make clean` - removes build intermediates and products
* `make python` - builds the Python library (see below)

## Python Library

A Python library called `peptide_design` is provided using [Boost Python](https://www.boost.org/doc/libs/1_70_0/libs/python/doc/html/index.html), a library that enables bindings between Python and C++ code. The Makefile instructions for the library require Python 3.8 (instructions using a `conda` environment are given below), and the library requires also building the `MST` Python library, which can be done by following the same install instructions below.  

Building the Python library for the first time is a somewhat convoluted process. The steps below are tested on MacOS, and should also work on Linux/Anthill.

1. Download and extract Boost (at the time of writing, the latest version is [1.73.0](https://www.boost.org/users/history/version_1_73_0.html)).

2. Unzip the download and move the `boost_1_73_0` folder to `/usr/local`. Enter the directory with  `cd /usr/local/boost_1_73_0`. Note: if you need to install at a custom location, see below.

3. Activate your Python environment with `conda`. The Makefile requires that you use an environment with Python 3.8.

4. Set an environment variable to point to your Python path: 
   
   ```
   $ PYTHON_TO_USE=$(which python)
   ```

5. Run the bootstrap script to setup the Boost build engine: 
   
   ```
   $ ./bootstrap.sh --prefix=/usr/local --show-libraries --with-python=$PYTHON_TO_USE --with-libraries=python
   ```

6. Finally, run the install script:
   
   ```
   $ ./b2 install --with-python
   ```

7. Now you should be able to run `make python` on either the MST repo or this repo, to build a shared library that incorporates the Python symbols.

### Help, things went wrong

The most common issues with the Python library seem to be:

1. **Python interpreter mismatch.** Make sure that you use the same Python 3.8 interpreter when you install Boost as when you build the library. If you don't, you may see errors like `python3.8-config: Command not found` or `'pyconfig.h' file not found`. This may also cause import errors when trying to load the module into Python. 
2. **Not finding all required symbols.** Building the Python library involves three steps: (1) building the `libpeptide_design.a` library containing all the symbols in the project, (2) building the `python.o` file containing the Boost Python bindings, and (3) building the `peptide_design.so` shared object, which can be imported into Python. Building the `python.o` object requires including the Boost Python headers (which should be in `/usr/local/include` following the steps above), as well as the Python development headers (provided by `python3.8-config --includes`). Building the shared object requires including the MST and `peptide_design` C++ libraries, the Boost Python library (should be in `/usr/local/lib`), and the Python development library. Absences of any of these files at the correct places can introduce errors, so you may need to look at the instructions printed by the Makefile and compare them to where the files exist on your system. Then you can adjust the Makefile for your specific configuration.
3. **Boost installed at custom location.** If you don't have permission to install at `/usr/local`, you'll need to specify an alternate location. This can be controlled when running `bootstrap.sh`by including `--prefix=/path/to/install`. Later on, when building the shared library/importing it to a python session, the system will need to know where you installed boost. This can be achieved by modifying the following environment variables
   `CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH:/path/to/install/include`
   `LIBRARY_PATH=$LIBRARY_PATH:/path/to/install/lib`
   `LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/path/to/install/lib`
   If you add these to `.bash_profile`, then your environment variables will always include these paths.

When things go wrong, it may be helpful to look at the symbols included in each intermediate file of the build process. The `nm` tool on MacOS/Linux can be helpful for this purpose - simply call `nm [file]` on a `.o`, `.a`, or `.so` file and pipe the output into a `grep`. The single-letter code next to each symbol can be a helpful clue - see [the man page](https://sourceware.org/binutils/docs/binutils/nm.html) for more info on these symbol types.  

Another strange error you may come across while building Boost is the following Clang message: `unknown argument -fcoalesce-templates`. I found that [this solution](https://alice-talk.web.cern.ch/t/o2-build-failed-at-boost-due-to-unknown-argument/545/2) fixed the issue.

# Programs

Below is a description of the main programs in the repository and the parameters they take, in the order that they might be run in the pipeline.

### `generateSeeds`

This program generates seeds in the space around selected residues.

Seeds are generated by 1) defining **fragments** (discontinous segments of protein backbone) from the protein residues near the binding site, 2) searching these for matches in a database of protein structures, and 3) searching within each match protein for interacting segments. The segments are extracted from the match protein and placed in the space around the target protein. Their exact pose is determined by the transformation that aligns the original fragment to the match protein. In summary, this process is referred to as "TERM Extension".

The parameters to this program are as follows:

- `--pdb` *(required)*: a peptide-protein complex structure in the `.pdb` format.

- `--peptide` : the chain ID of the peptide in the `.pdb` file, e.g. "B". Note: must be a single chain.

- `--sel`: an alternative to `--peptide`, a selection string specifying the protein residues that should be used to define fragments.

- `--flanking_res` *(required)*: when defining fragments from a set of residues, will include this many residues +/- the selected residues. (recommended: 2)

- `--config` *(required)*: the configuration file that specifies the FASSTDB/rotamer library (same the configuration file used by dTERMen).

- `--match_req`: The maximum number of matches per fragment. If selected, will use a greedy algoritm to search for the largest possible protein fragments (defined by number of residues) that still have this many matches.

- 

### `findOverlaps`

This program performs an all-to-all pairwise comparison of the provided seeds, and assembles a list of mappings of residues between seeds such that the seeds can be considered to connect from one to the other. At most one overlap is provided per pair of seeds, and by default overlaps can occur at the ends or in the middle of a seed. `findOverlaps` makes use of the `FuseCandidateFinder` class (defined in `findfuseable.h/cpp`), and modifying the initialization of this object in `findOverlaps.cpp` can allow for different behaviors.

The parameters to this program are as follows:

* `--files` *(required)*: a text file containing a list of seed structures, relative to the `data` path parameter. This text file can be easily generated by running the following command in the `data` directory, if the only PDB files it contains are seeds: `find . -name *.pdb`.
* `--data` *(required)*: path to the directory in which seed structures are stored.
* `--out` *(required)*: path to a directory into which to write overlaps.
* `--overlapSize`: number of residues that must overlap between two residues (default 2).
* `--overlapRMSD`: maximum RMSD between sets of residues to be considered overlapping (default 1.0). This is not an adaptive threshold, so it should be passed in tandem with `overlapSize`.
* `--minCosAngle`: minimum cosine angle between Ca vectors allowed in a potential overlap (default -1.0, or no restriction). This helps to ensure that overlaps are pointing in the same direction, especially if the RMSD threshold is high.
* `--batch`: the batch number (possible values range from 1 to `numBatches`)
* `--numBatches`: the total number of batches. The `batch` and `numBatches` parameters are used to divide up the pairwise comparisons if the program is being run as an array job on a cluster.

### `buildSeedGraphs`

This program builds an adjacency list for the overlaps from the output of `findOverlaps`. It also finds connected neighborhoods in the graph and writes them out in "chunks," which helps to optimize the `scoreSeeds` step. The adjacency list can be built in two different ways: "bond", meaning an edge signifies a possible bond between the two connected residues; and "same", in which an edge signifies a direct spatial overlap between the connected residues. The "bond" mode is therefore useful for constructing paths of residues, while the "same" mode is useful for reducing downstream computation on similar residues. *Note:* because the "same" mode is used mainly for caching purposes, you may want to run `findOverlaps` with a more stringent RMSD cutoff to generate the "same" graph, then use a more permissive cutoff for the "bond" graph.

The parameters to `buildSeedGraph` are as follows:

* `--overlaps` *(required)*: path to a directory containing overlaps from the previous step.
* `--seeds` *(required)*: text file listing all seeds, i.e. the `--files` argument from above.  
* `--data` *(required)*: directory in which seed structures are stored, as above.
* `--out` *(required)*: directory into which to write seed graphs.
* `--adj`, "If 'same' (default) then write graphs where adjacencies correspond to equivalent residues; if 'bond', write graphs where adjacencies are potential bonds.
* `--chunkSize`: the number of residues required to write out a chunk file (default is 100 if adj = 'same', 15 otherwise). The larger this number, the fewer files will be written out. For adj = 'same', it should be adjusted so that the number of files is roughly a multiple of the number of parallel jobs used in `scoreSeeds`.

The output of this program consists of a `whole_graph.txt` file defining the adjacencies in the entire graph; some `subgraph_[number].txt` files containing connected subgraphs larger than a certain threshold; and `chunk_[number].txt` files containing smaller connected subgraphs.

### `scoreSeeds`

This program scores seed residues in an easily parallelizable way. In its current implementation, this is by far the performance bottleneck of the pipeline, since pair-fragment FASST searches take 60-90 seconds each. Therefore, this program implements a few strategies for reducing the amount of necessary computation to score all seed residues.

The parameters to this program are as follows:

* `--target` *(required)*: the path to the target PDB structure.
* `--graphs`: directory containing seed graphs generated with adj = 'same'. These are used to cluster similar residues and reduce the number of redundant computations.
* `--data` *(required)*: directory in which seed structures are stored.
* `--contacts`: directory into which to write contact counts, if in contact counting mode; or from which to read the counts if in scoring mode. Required for `mode` = 1.
* `--out`: directory into which to write seed scores. *Note*: it is recommended that you create the directory before starting the programs, to avoid race conditions between parallel executions.
* `--mode`: 0 to score seeds, or 1 to count contacts for each target residue. Contact counting was used in an earlier version of the scoring function, but it is now considered unnecessary.
* `--fasstDB`: path to the FASST database to be used in scoring searches. A default database path on the `anthill` cluster is used if this parameter is not provided.
* `--batch`: the batch number (possible values range from 1 to `numBatches`).
* `--numBatches`: the number of batches, corresponding to the number of parallel jobs if using an array on the cluster.
* `--numSeedFlank`: the number of residues on either side of the seed residue to use for scoring (default 2).
* `--numTargetFlank`: the number of residues on either side of the target residue to use for scoring (default 2).
* `--append`: if true, read any existing files numbered `batchIndex, batchIndex + numBatches, batchIndex + 2 * numBatches, ...` and determine which seeds have already been scored, then create a new file that covers residues not already scored for this batch index.

The output of this program consists of a `seed_scores_[batchIndex].csv` file, containing the seed residue (in the format `seed_name:residue_index`) and the seed score. The value is set to a very high number if the seed is not designable or if it clashes with the target. The program also outputs a file called `frag_scores_[batchIndex].csv` in a subdirectory `[outPath]/all_fragment_scores`, which contains all the components of the score for each fragment that was scored.

### `findPaths`

This program uses a dynamic programming algorithm to find the best-scoring paths in a seed graph. It then fuses these paths and writes out a PDB structure for each one. The parameters it takes are as follows:

* `--mode` *(required)*: the algorithm to use for finding paths. Currently only "DP" (the default) is supported.
* `--target` *(required)*: path to the original protein target PDB file. This is used to provide an anchor for the fuser.
* `--graphs` *(required)*: path to a graph file to traverse (generated with adj = bond).
* `--data`: directory in which seed structures are stored.
* `--scores`: directory containing seed scores, output by `scoreSeeds`.
* `--out`: path to a directory into which the best seed path structures will be written.
* `--numPaths`: number of paths to generate (linearly impacts running time - default 500).

The output of `findPaths` in the `out` directory consists of a `paths.txt` file containing the best-scoring sequences of residues as well as their estimated scores; and `fusion_[number].pdb` files containing the fused structures in the same order as specified in `paths.txt`.

## Future Directions

1. I am currently generating a clustered database using a variation of Craig's `structgen` repo, which should greatly speed up seed scoring tasks. It would be best to bring the `ClusterDatabase` class into this repository and adapt its search function to ensure the correct residue-to-residue alignment.
2. The `findPaths` program would benefit from a randomized approach, i.e. sampling random paths from the graph and scoring them post-hoc. This is because the fusing process with high RMSD cutoffs tends to alter sequence scores dramatically, and the so-called "optimal" paths output by the DP algorithm may no longer be the best.
3. If the randomized approach is implemented, better score functions are needed to rank the generated paths, measuring quantities such as compactness and designability.
