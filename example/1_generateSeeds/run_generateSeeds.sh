#!/bin/bash
#
#SBATCH -J generateseeds
#SBATCH -o generateseeds.%J.out
#SBATCH -e generateseeds.%J.err
#SBATCH -t 1-00:00:00
#SBATCH --mem=20G
#SBATCH -p defq

peptide_design=/scratch/users/swans/MST_workspace/peptide_design
pdb=../input_files/1LB6_A__.pdb
params_file=genSeeds.params

SECONDS=0

srun $peptide_design/bin/generateSeeds --pdb $pdb --params_file $params_file

ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo $ELAPSED

exit 0
