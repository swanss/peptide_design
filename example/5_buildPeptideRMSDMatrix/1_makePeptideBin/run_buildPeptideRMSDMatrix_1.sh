#!/bin/bash
#
#SBATCH -J buildPeptideRMSDMatrix
#SBATCH -o buildPeptideRMSDMatrix.%A.out
#SBATCH -e buildPeptideRMSDMatrix.%A.err
#SBATCH -p defq
#SBATCH -t 1-00:00:00
#SBATCH -n 1
#SBATCH --mem=5G

peptide_design=/scratch/users/swans/MST_workspace/peptide_design

# get a list of all peptide structures
ls ../../4_samplePaths/path_structures/*_fused-path_*.pdb > peptideStructures.list
list=peptideStructures.list

SECONDS=0

srun $peptide_design/bin/buildPeptideRMSDMatrix --list $list

ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo $ELAPSED

exit 0

