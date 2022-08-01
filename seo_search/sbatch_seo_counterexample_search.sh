#!/bin/bash

#SBATCH -J seo_counterexample_search_flip_only_once
#SBATCH -a 000-035
#SBATCH --cpus-per-task 1
#SBATCH --output=slurmreports_flip_only_once/slurm-%A_%2a.out
#SBATCH -p IGIcuda4

SATID=`printf %03d $SLURM_ARRAY_TASK_ID`

#OUTDIR=/dev/null

date
../target/release/seo_bt_flip_only_once -s $SATID -n 7 -p 0.8

date
