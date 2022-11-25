#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --mem=64g
#SBATCH --time=1:00:00
#SBATCH --tmp=10g
#SBATCH --partition=msismall,msilarge,msibigmem
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=malak039@umn.edu
#SBATCH -o logs/%j.out
#SBATCH -e logs/%j.err

module load R/4.1.0

cd ${DRAB}

awk '(NR > 6) && ($1 ~ /^chr[1-9]/) && ($3 == "gene") {print $16 "\t" $10 "\t" $1 "\t" $4 "\t" $5 "\t" $14}' raw/gencode.v26.GRCh38.genes.gtf | \
          tr -d ";\"" > annotations/all_genes.txt

awk '$6 == "protein_coding" {print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5}' annotations/all_genes.txt > annotations/protein_coding_genes.txt

./plink --vcf raw/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.vcf.gz \
          --const-fid \
          --make-bed \
          --out TEMP1

Rscript --vanilla util/extract_phenotypes.R

./plink --bfile TEMP1 \
          --keep phenotypes.txt \
          --make-bed \
          --out TEMP2

./plink --bfile TEMP2 \
          --autosome \
          --make-bed \
          --out TEMP3

./plink --bfile TEMP3 \
          --snps-only \
          --make-bed \
          --out TEMP4

./plink --bfile TEMP4 \
          --geno 0 \
          --make-bed \
          --out TEMP5

./plink --bfile TEMP5 \
          --hwe 1e-6 midp \
          --make-bed \
          --out TEMP6

./plink --bfile TEMP6 \
          --maf 0.01 \
          --make-bed \
          --out genotypes/dosages

rm TEMP*

Rscript --vanilla util/GTEx_to_plink.R

