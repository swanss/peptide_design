#!/bin/bash
#
#SBATCH -J submitFPD
#SBATCH -o submitFPD_repeat.%A.out
#SBATCH -e submitFPD_repeat.%A.err
#SBATCH -p defq
#SBATCH -t 1-00:00
#SBATCH -n 1
#SBATCH --mem=2G

proj_dir=$PWD
structures_list=selectedRelaxedPeptideDesigns.txt

mkdir jobs
cd jobs
while IFS= read -r path; do
	pdb=${path##*/}
	pdb_name=${pdb%.pdb}
	
	echo $pdb_name
	mkdir $pdb_name
	cd $pdb_name
	for c in 10 15 20; do
	#for c in 10 15 20 25; do
		dirname=$c"_outerCycles"
		mkdir $dirname
		cd $dirname
		echo sbatch $proj_dir/run_FPDRefine.sh $path $c
		sbatch $proj_dir/run_FPDRefine.sh $path $c
		cd ..
	done
	cd ..
done < $proj_dir/$structures_list
exit 0

