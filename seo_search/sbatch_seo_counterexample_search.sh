#!/bin/bash

#SBATCH -J seo_counterexample_exact_search_small
#SBATCH -a 000-049
#SBATCH --cpus-per-task 1
#SBATCH --output=slurmreports_exact_small/slurm-%A_%2a.out
#SBATCH -p IGIcuda4

SATID=`printf %03d $SLURM_ARRAY_TASK_ID`

#OUTDIR=/dev/null

date
../target/release/seo_search_counterexample -s $SATID -n 10 -p 0.8 -l 1000000 -x

date
