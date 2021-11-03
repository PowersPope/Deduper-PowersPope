#!/bin/bash
#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --job-name=deduper
#SBATCH --output=deduper-%j-sbatch.out
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=2-00:00:00




/usr/bin/time -v python python_scripts/Deduper_nolist.py -f datasets/C1_SE_uniqAlign_sort.sam -u STL96.txt -o C1_SE_uniqAlign_sort
