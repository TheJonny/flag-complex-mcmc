#!/bin/bash

#SBATCH -J seo_counterexample_search_clique_strictly_monotonic
#SBATCH -a 000-035
#SBATCH --cpus-per-task 1
#SBATCH --output=slurmreports_clique_strictly_monotonic/slurm-%A_%2a.out
#SBATCH -p IGIcuda4

SATID=`printf %03d $SLURM_ARRAY_TASK_ID`

#OUTDIR=/dev/null

date
../target/release/seo_search_counterexample -s $SATID -n 7 -p 2.0 -l 50000 -x

date
