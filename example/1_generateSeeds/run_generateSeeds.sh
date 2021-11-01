#!/bin/bash
#
#SBARCH -J generateSeeds
#SBATCH -o generateSeeds.%J.out
#SBATCH -e generateSeeds.%J.err
#SBATCH -t 1-00:00:00
#SBATCH --mem=20G
#SBATCH -p defq

peptide_design=/scratch/users/swans/MST_workspace/peptide_design

targetPDB=../input_files/1LB6_A__.pdb
paramsFile=genSeeds.params

SECONDS=0

srun $peptide_design/bin/generateSeeds --targetPDB $targetPDB --paramsFile $paramsFile

ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo $ELAPSED

exit 0
