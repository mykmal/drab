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

module load R/4.1.0

cd ${DREX}

for f in expression_covariates/*.txt
do
cat ${f} | head -n1 | tr "\t" "\n" | tail -n+2 > ${f}.HEADER
done

Rscript code/ProcessCovariates.R
