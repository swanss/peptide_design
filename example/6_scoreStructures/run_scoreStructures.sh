#!/bin/bash
#
#SBATCH -a 1-50
#SBATCH -J scoreStructures
#SBATCH -o scoreStructures.%A_%a.out
#SBATCH -e scoreStructures.%A_%a.err
#SBATCH -t 1-00:00:00
#SBATCH --mem=16G
#SBATCH -p defq
#SBATCH -n 1

peptide_design=/scratch/users/swans/MST_workspace/peptide_design

target=../input_files/1LB6_A__.pdb
structures=../5_scoreStructures #this should be a path to a directory containing PDB files of peptide structures to be scored which contains a file "structures.list" where each line is the filename of one of the structures 
out=scores
config=../input_files/multichainDB.configfile
worker=$SLURM_ARRAY_TASK_ID
numWorkers=50

SECONDS=0

srun $peptide_design/bin/scoreStructures --target $target --structures $structures --out $out  --config $config --worker $worker --numWorkers $numWorkers

ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo $ELAPSED

exit 0

