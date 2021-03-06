#!/bin/bash
# Job name:
#SBATCH --job-name=coping
#
# Account:
#SBATCH --account=fc_coping
#
# Partition:
#SBATCH --partition=savio2
#
# Tasks per node:
# SBATCH --ntasks-per-node=4
#
# Processors per node:
#SBATCH --cpus-per-task=1
#
# Mail type:
#SBATCH --mail-type=all
#
# Mail user:
#SBATCH --mail-user=psmits@berkeley.edu
#
# Wall clock limit:
#SBATCH --time=500:00:00
#
# Working directory
#SBATCH --workdir=/home/users/psmits/coping/src/
#
## Command(s) to run:
bash trait_model_wrap.sh
