#!/bin/bash
#SBATCH -J ga_ffield
#SBATCH -p bulldozer_64c
#SBATCH -N 1
#SBATCH -t 00:00:00
#SBATCH --comment qe
# The 12 chromosomes of a generation are evaluated in parallel on one node.
export LMP=/home/mino/02.programs/lammps-23Jun2022_RexPoN_101/src/lmp_mpi
../../scripts/run_ga.sh 40000
