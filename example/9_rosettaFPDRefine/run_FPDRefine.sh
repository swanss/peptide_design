#!/bin/bash
#
#SBATCH -J fpdRefine
#SBATCH -o fpdRefine.%A.out
#SBATCH -e fpdRefine.%A.err
#SBATCH -p defq
#SBATCH -t 2-00:00
#SBATCH -n 1
#SBATCH --mem=8G

complex=$1
nstruct=200
peptide_chain="B"
rec="A"
rep_ramp_cycles=$2
flags=$PWD/refine_flags

SECONDS=0

mkdir output

echo srun /scratch/users/swans/rosetta/rosetta_src_2018.33.60351_bundle/main/source/bin/FlexPepDocking.default.linuxgccrelease -in:file:s $complex -native $complex -nstruct $nstruct -flexPepDocking:peptide_chain $peptide_chain -flexPepDocking:receptor_chain $rec -rep_ramp_cycles $rep_ramp_cycles @ $flags
srun /scratch/users/swans/rosetta/rosetta_src_2018.33.60351_bundle/main/source/bin/FlexPepDocking.default.linuxgccrelease -in:file:s $complex -native $complex -nstruct $nstruct -flexPepDocking:peptide_chain $peptide_chain -flexPepDocking:receptor_chain $rec -rep_ramp_cycles $rep_ramp_cycles @ $flags

ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo $ELAPSED

exit 0
