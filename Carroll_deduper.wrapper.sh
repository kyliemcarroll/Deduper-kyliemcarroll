#!/bin/bash

#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --cpus-per-task=16
#SBATCH --time=33:33:33
#SBATCH --job-name=DEDUPERRRR
#SBATCH --output=deduper_%j.out

input="/projects/bgmp/kymc/bioinfo/Bi624/C1_SE_uniqAlign.sorted.sam"
output="/projects/bgmp/kymc/bioinfo/Bi624/C1_SE_uniqReads.sorted.sam"
umi="/projects/bgmp/kymc/bioinfo/Bi624/STL96.txt"

/usr/bin/time -v python /projects/bgmp/kymc/bioinfo/Bi624/Carroll_deduper.py -f $input -o $output -u $umi