#!/bin/bash
#
#SBATCH --time=0-12:00:00
#SBATCH --job-name=stat
#SBATCH --output=%x-%j.out  #where to write terminal output
#SBATCH --mem=16G           #number of RAM per node
#SBATCH -N 1                #number of nodes
#SBATCH -n 4               #number of cores

python3 /home/zhenlun/scratch/medras_actual_new/newattempt.py >> cluster_single.txt
