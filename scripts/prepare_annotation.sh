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

cd ${DREX}

awk '(NR > 6) && ($1 ~ /^chr[1-22]/) && ($3 == "gene") {print $16 "\t" $10 "\t" $1 "\t" $4 "\t" $5 "\t" $14}' reference/gencode.v26.GRCh38.genes.gtf | \
tr -d ";\"" > reference/gene_annotation.txt

awk '$6 == "protein_coding" {print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5}' reference/gene_annotation.txt > reference/pc_gene_annotation.txt
