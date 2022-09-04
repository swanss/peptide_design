#!/bin/bash
#
#SBATCH -J findOverlaps
#SBATCH -o findOverlaps.%A.out
#SBATCH -e findOverlaps.%A.err
#SBATCH -p defq
#SBATCH -t 1-0:00:00
#SBATCH -n 1
#SBATCH --mem=4G

peptide_design=/scratch/users/swans/MST_workspace/peptide_design

seedBin=../1_generateSeeds/output/extendedfragments.bin
out=overlaps
overlapSize=4
deviation=1.0
minCosAngle=0.5
batchSize=10000

SECONDS=0

$peptide_design/bin/findOverlaps --seedBin $seedBin --out $out --overlapSize $overlapSize --deviation $deviation --minCosAngle $minCosAngle --batchSize $batchSize

ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo $ELAPSED

exit 0

