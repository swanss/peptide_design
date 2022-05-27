#!/bin/bash
#
#SBATCH -J relax 
#SBATCH -o relax.%A.out
#SBATCH -e relax.%A.err
#SBATCH -p defq
#SBATCH -t 1-00:00
#SBATCH -n 1
#SBATCH --mem=8G

#Runs fast relax as per https://www.rosettacommons.org/docs/latest/application_documentation/structure_prediction/relax
rosetta_bin_dir=/scratch/users/swans/rosetta/rosetta_src_2018.33.60351_bundle/main/source/bin

protein= #path to structure
nstruct=10 #select the structure with the lowest energy

SECONDS=0

echo srun $rosetta_bin_dir/relax.default.linuxgccrelease -relax:default_repeats 5 -in:file:s $protein -nstruct $nstruct  
srun $rosetta_bin_dir/relax.default.linuxgccrelease -relax:default_repeats 5 -in:file:s $protein -nstruct $nstruct  

ELAPSED="Elapsed: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo $ELAPSED

exit 0
