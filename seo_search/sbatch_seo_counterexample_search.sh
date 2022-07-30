#!/bin/bash

#SBATCH -J seo_counterexample_search_strictly_monotonic_2_only
#SBATCH -a 000-035
#SBATCH --cpus-per-task 1
#SBATCH --output=slurmreports_strictly_monotonic_2/slurm-%A_%2a.out
#SBATCH -p IGIcuda4

SATID=`printf %03d $SLURM_ARRAY_TASK_ID`

#OUTDIR=/dev/null

date
../target/release/seo_search_counterexample -s $SATID -n 5 -p 0.8 -l 50000 

date
