#!/bin/bash
#
#SBATCH -J dTERMen
#SBATCH -o dTERMen.%A.out
#SBATCH -e dTERMen.%A.err
#SBATCH -t 1-00:00
#SBATCH -n 1
#SBATCH --mem=10G
#SBATCH -p defq

#this is just provided an example to help you get started

MST_install=""

p="" # the PDB structure of the designed peptide and target protein
c=../input_files/dTERMen.configfile #see Mosaist/fasstDB to construct a fasstDB for dTERMen
s="chain ?" # select the peptide by chain
o="test" # the name of the complex

SECONDS=0

echo $MST_install/MST/bin/design --p $p --c $c --s "$s" --o $o
srun /scratch/users/swans/MST_workspace/MST/bin/design --p $p --c $c --s "$s" --o $o

ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo $ELAPSED

exit 0

