# peptide_design

This repository incorporates several tools that enable the de novo design of peptides to bind a given protein. The process begins with the generation of  **seeds**, which are small segments of protein backbone arranged in the space around the selected binding site. Given these seeds, the programs enable the construction of a **seed graph** describing geometric overlaps between these seeds. Paths can be sampled from the graph and fused to form peptide backbones. Finally, paths can be scored and selected for sequence design.

## Build Instructions

Before building, adjust the `makefile` variable `MSTDIR` to be the path at which MST is stored on your system. The default location is in the parent directory of `peptide_design`.

* `make all` - builds all programs
* `make test` - builds only programs in the `tests` directory
* `make bin/[executable name]` - builds the specific executable with its dependencies
* `make clean` - removes build intermediates and products
* `make python` - builds the Python library (see below)

## Main pipeline

See `peptide_design/example/` for an example of how to use this pipeline to design peptide backbones.

Before starting, you will need to make a configuration file, which will be reused throughout the process. This will provide the path to 1) a FASST file, i.e. a database of structures that have been processed and can be searched and 2) the backbone-dependent rotamer library (which can be found in the MST repo). ex: `peptide_design/example/input_files/singlechain.configfile`

### `generateSeeds`

`peptide_design/example/01_generateSeeds/run_generateSeeds.sh`

The details of how seeds are generated are controlled through the params file. Anything that is not included in this file
will be set to the default value, with the exception of **config_file** which must always be provided.

```
fragment_type ADAPTIVE_LENGTH #fragments will grow in length until they have less matches than match_req
max_rmsd 1.2 #the cutoff for defining matches to fragments
flanking_res 3 #the max number of flanking residues that can be added when defining a fragment
match_req 5000 #the number of matches that must be found for each target residue
adaptive_rmsd 0 #if 1, will scale the max RMSD cutoff based on fragment size/complexity
seq_const NONE #if NONE, no sequence constraint is applied when searching for matches
config_file /scratch/users/swans/config/singlechain.configfile
seed_flanking_res 2 #the number of flanking residues that are included around every central seed residue
allow_sidechain_clash 0 #if 1, will allow seeds that clash with target sidechains
relSASA_cutoff -1.0 #this threshold defines 'surface residues' that are used to generate seeds (ignored if -1.0)
verbose 1
```
When complete, all seeds will be stored at `dir/output/extendedfragments.bin`

### `findOverlaps`

`peptide_design/example/02_findOverlaps/run_findoverlaps.sh`

Depending on the number of seeds, finding overlaps can be slow. Options are provided to distribute the work via job arrays.

### `buildSeedGraph`

`peptide_design/example/02_findOverlaps/run_buildSeedGraph.sh`

### `samplePaths`

### `countContacts`

### `buildPathDistanceMatrix`

### `scoreStructures`

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

- `--params_file` *(required)*: a file specifying all parameter values 

- `--config` *(required)*: the configuration file that specifies the FASSTDB/rotamer library (same the configuration file used by dTERMen).

`params_file` formatting:

Each line consists of two fields: the name and value of each parameter separated by a space. Parameters may be provided in any order. Example:

```
param_name1 param_value1
param_name2 param_value2
```

Parameters that must be manually specified:

- `config_file` An additional parameters file specifying the location of the protein structure database and rotamer library. Example:

```
fasstdb=/path/to/dbrotlib=/path/to/db
rotlib=/path/to/db
```

Parameters with default values that may be tweaked

- `fragment_type` The method that is used to generate binding site fragments. May be either CEN_RES, ALL_COMBINATIONS, ADAPTIVE_SIZE, ADAPTIVE_LENGTH (default value), or COMPLEXITY_SCAN.

- `cd_threshold` The threshold used to define sidechain-sidechain interactions. May be between [0,1], default value = 0.01.

- `int_threshold` The threshold used to define sidechain-backbone interactions. May be between [0,1] , default value = 0.01.

- `bbInteraction_cutoff` The cutoff used to define backbone-backbone interactions. Default value = 3.5 angstroms

- `max_rmsd` When `adaptive_RMSD` is set to off, then this is just the RMSD cutoff used to define structural matches to a binding site fragment.  Default value = 1.2 angstroms

- `flanking_res` This is the maximum number of residues +/- a binding site residue that could be included in the fragment. Default value = 2

- `match_req` The number of matches required for each binding site fragment. Only applies when the `ADAPTIVE_X` `fragment_type` modes are used. Fragments will grow until they can not incorporate more residues without having less than the required number of matches. Not applied when less than 0. Default value = -1

- `adaptive_rmsd` When true, then the max RMSD cutoff is scaled by a factor that is determined by the size and topology of the fragment. See [https://doi.org/10.1073/pnas.1607178113](https://doi-org.libproxy.mit.edu/10.1073/pnas.1607178113) for more information. Default = 0.

- `seq_const` Constrains the structural matches to those that match the sequence of the query fragment. Modes include NONE (default), CEN_RES_ONLY, or ALL_RES

- `seed_flanking_res` The number of residues +/- a central seed residue that are included as structural context in the seed. Default = 2.

- `homology_cutoff` The sequence identity cutoff used to determine whether a structural match should be considered a homolog. The sequence of the aligned region and up to 30 residues +/- this region are considered. Default = 0.5

- `allow_sidechain_clash` Controls whether seeds (which consist of backbone atoms only) are allowed to clash with the sidechain atoms of the protein. If false, then clashing seeds will be removed. Default = 1.

- `verbose` Controls whether extra information is output. Useful for debugging. Default = 0.

```
fasstdb=/path/to/db
rotlib=/path/to/db
```

- `fragment_type` *(required)*: when defining fragments from a set of residues, will include this many residues +/- the selected residues. (recommended: 2)

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

