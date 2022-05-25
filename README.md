# peptide_design

[![Paper DOI : 10.1002/pro.4322](https://badgen.net/badge/Protein%20Science%20DOI/10.1002%2Fpro.4322/black)](https://doi.org/10.1002/pro.4322) &nbsp; [![Zenodo](https://badgen.net/badge/zenodo/FASST%20DB%20download/red)](https://zenodo.org/record/6569429) &nbsp; [![Twitter](https://badgen.net/badge/icon/Sebastian%20Swanson?icon=twitter&label)](https://twitter.com/SassSeabass)

This repository incorporates several tools that enable the de novo design of peptides to bind a given protein. The process begins with the generation of  **seeds**, which are small segments of protein backbone arranged in the space around the selected binding site. Given these seeds, the programs enable the construction of a **seed graph** describing geometric overlaps between these seeds. Paths can be sampled from the graph and fused to form peptide backbones. Finally, paths can be scored and selected for sequence design.

## Build Instructions

Peptide design has the following dependencies

### Mosaist

1. Clone the repo https://github.com/Grigoryanlab/Mosaist
2. Follow the instructions to install

### FreeSASA

1. Go to http://freesasa.github.io/ , download the latest tarball into the same directory containing the `peptide_design` repo and extract.
2. Run `./configure --disable-threads --disable-xml --disable-json` 
3. Run `make all`. You should see that the library `libfreesasa.a` is created.

### Building `peptide_design`

Update the `peptide_design` makefile so that:

1. `MSTDIR` is the path to Mosaist
2. `SASADIR` is the path to FreeSASA 

The default location for both is in the parent directory of `peptide_design`.

* `make all` - builds all programs
* `make test` - builds only programs in the `tests` directory
* `make bin/[executable name]` - builds the specific executable with its dependencies
* `make clean` - removes build intermediates and products
* `make python` - builds the Python library (see below)

## Main pipeline

See `peptide_design/example/` for an example of how to use this pipeline to design peptide backbones.

Before starting, you will need to make a configuration file, which will be reused throughout the process. This will provide the path to 1) a FASST file, i.e. a database of structures that have been processed and can be searched and 2) the backbone-dependent rotamer library (which can be found in the MST repo). ex: `peptide_design/example/input_files/singlechain.configfile`. NOTE: FASST databases are available for [download](https://zenodo.org/record/6569429).

### generateSeeds

Generates segments of protein backbone, or 'interface seeds', around a target protein.

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

### findOverlaps

Finds overlaps between all pairs of seeds.

`2_findOverlaps/run_findoverlaps.sh`

Depending on the number of seeds, finding overlaps can be slow. Options are provided to distribute the work via job arrays.

### buildSeedGraph

Builds a graph describing seed residues and their potential connections.

`3_buildSeedGraph/run_buildSeedGraph.sh`

### samplePaths

Samples random paths from a seed graph and fuses the residues together into a peptide backbone structure.

`4_samplePaths/run_samplePaths.sh`

By default this will not find contacts between the designed peptide backbones and the target (as this is fairly slow), but this can be modified with the `--countContacts` option.

### buildPeptideRMSDMatrix

Computes RMSD between peptide backbones in a parallelizable manner and builds a complete distance matrix for downstream analysis.

This program must be run three times to generate a complete distance matrix.

`5_buildPeptideRMSDMatrix/1_makePeptideBin/run_buildPeptideRMSDMatrix_1.sh`
Combines all peptide structures into a single binary file.

`5_buildPeptideRMSDMatrix/2_computeRMSD/run_buildPeptideRMSDMatrix_2.sh`
Computes the RMSD between all pairs of peptide structures, with support for job arrays.

`5_buildPeptideRMSDMatrix/3_buildMatrix/run_buildPeptideRMSDMatrix_3.sh`
Combines the output of each job into a single distance matrix.

### scoreStructures

Scores the interface formed between a set of peptides and the target protein. 

`6_scoreStructures/run_scoreStructures.sh`

The score for each residue of each peptide backbone structure is written to the `structure_scores_*.tsv*` file. The score is defined such that negative scores are more favorable. In practice, we generally take the average over all residues in the peptide to get the overall score.

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
