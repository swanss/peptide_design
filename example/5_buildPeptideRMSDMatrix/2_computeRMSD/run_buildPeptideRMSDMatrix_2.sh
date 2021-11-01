#!/bin/bash
#
#SBATCH -a 1-20
#SBATCH -J buildPeptideRMSDMatrix
#SBATCH -o buildPeptideRMSDMatrix.%A_%a.out
#SBATCH -e buildPeptideRMSDMatrix.%A_%a.err
#SBATCH -p defq
#SBATCH -t 1-00:00:00
#SBATCH -n 1
#SBATCH --mem=5G

peptide_design=/scratch/users/swans/MST_workspace/peptide_design

bin=../1_makePeptideBin/*bin
numWorkers=20
worker=$SLURM_ARRAY_TASK_ID

SECONDS=0

srun $peptide_design/bin/buildPeptideRMSDMatrix --bin $bin --numWorkers $numWorkers --worker $worker

ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo $ELAPSED

exit 0

