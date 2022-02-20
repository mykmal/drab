#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64g
#SBATCH --time=30:00
#SBATCH --tmp=100g
#SBATCH -p amdlarge,amd512,amd2tb,ram256g,ram1t
#SBATCH --mail-type=ALL
#SBATCH --mail-user=scientist@university.edu
#SBATCH -o logs/%j.out
#SBATCH -e logs/%j.err

DREX="/path/to/drex"

TISSUE="Whole_Blood"

module load R/4.1.0

cd ${DREX}

Rscript code/SplitData.R ${TISSUE}

gzip expression_matrices/*.bed

