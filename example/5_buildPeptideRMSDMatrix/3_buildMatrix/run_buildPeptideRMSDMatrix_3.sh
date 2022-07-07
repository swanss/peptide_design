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

bin=../1_makePeptideBin/*bin

#get distances list
ls ../2_computeRMSD/distances*.csv > distances.list
distanceList=distances.list

SECONDS=0
srun $peptide_design/bin/buildPeptideRMSDMatrix --bin $bin --distanceList $distanceList
ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo $ELAPSED

exit 0

