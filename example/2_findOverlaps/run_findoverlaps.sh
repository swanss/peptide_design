#!/bin/bash
#
#SBATCH -a 1-20
#SBATCH -J findoverlaps
#SBATCH -o findoverlaps.%A_%a.out
#SBATCH -e findoverlaps.%A_%a.err
#SBATCH -p defq
#SBATCH -t 0-03:00:00
#SBATCH -n 1
#SBATCH --mem=4G

peptide_design=/scratch/users/swans/MST_workspace/peptide_design

numWorkers=20
echo "job ID: "$SLURM_ARRAY_TASK_ID

seedBin=../1_generateSeeds/output/extendedfragments.bin
dir=./output
overlapSize=4
deviation=1.0
minCosAngle=0.5
batchSize=10000

out=$dir/overlaps$SLURM_ARRAY_TASK_ID

#make a directory for the output, if it doesn't already exist
if [ ! -d $dir ]; then
        mkdir $dir
fi

SECONDS=0

$peptide_design/bin/findOverlaps --seedBin $seedBin --out $out --overlapSize $overlapSize --deviation $deviation --minCosAngle $minCosAngle --worker $SLURM_ARRAY_TASK_ID  --numWorkers $numWorkers --batchSize $batchSize

ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo $ELAPSED

exit 0

