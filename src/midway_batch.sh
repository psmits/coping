#!/bin/bash
#SBATCH --job-name=example_sbatch
#SBATCH --output=example_sbatch
#SBATCH --error=example_sbatch
#SBATCH --time=00:05:00
#SBATCH --partition=broadwl
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=16
#SBATCH --mem-per-cpu=2000


# Execute the program
bash /psmits/coping/src/trait_model_wrap.sh
