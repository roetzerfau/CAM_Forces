#!/bin/bash –l
#
#SBATCH --nodes=10
#SBATCH --time=24:00:00
#SBATCH --mail-user=maximilian.roetzer@fau.de
#SBATCH --job-name=test_run_mcmcm
#SBATCH --export=none

unset SLURM_EXPORT_ENVmodule 
module load python


python python_examples/mcmc_2parameter_2D3D_tests.py 
