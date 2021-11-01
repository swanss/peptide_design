#!/bin/bash
#
#SBATCH -J buildSeedGraph 
#SBATCH -o buildSeedGraph.%A.out
#SBATCH -e buildSeedGraph.%A.err
#SBATCH -p defq
#SBATCH -t 1-00:00
#SBATCH -n 1 
#SBATCH --mem=40G

peptide_design=/scratch/users/swans/MST_workspace/peptide_design

seedBin=../1_generateSeeds/output/extendedfragments.bin
out="1LB6_seedGraph"
overlaps=../2_findOverlaps/output/

SECONDS=0

srun $peptide_design/bin/buildSeedGraph --seedBin $seedBin --out $out --overlaps $overlaps

ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo $ELAPSED

exit 0

