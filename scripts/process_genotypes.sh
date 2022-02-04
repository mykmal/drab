#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64g
#SBATCH --time=2:00:00
#SBATCH --tmp=100g
#SBATCH -p amdlarge,amd512,amd2tb,ram256g,ram1t
#SBATCH --mail-type=ALL
#SBATCH --mail-user=scientist@university.edu
#SBATCH -o %j.out
#SBATCH -e %j.err

export DREX="/path/to/drex"

module load R/4.1.0

export PATH=${PATH}:${DREX}/plink

cd ${DREX}

Rscript code/ExtractEUR.R

plink --vcf genotypes/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.vcf.gz \
          --const-fid \
          --vcf-require-gt \
          --biallelic-only strict \
          --make-bed \
          --out TEMP1

plink --bfile TEMP1 \
          --keep reference/EUR_samples.txt \
          --make-bed \
          --out TEMP2

plink --bfile TEMP2 \
          --autosome \
          --snps-only just-acgt \
          --make-bed \
          --out TEMP3

awk -F '\t' '/G_C|C_G|A_T|T_A/ {print $2}' TEMP3.bim > reference/ambiguous_snps.txt

plink --bfile TEMP3 \
          --exclude reference/ambiguous_snps.txt \
          --make-bed \
          --out TEMP4

awk -F '\t' '(NR>1) && ($7 ~ /^rs/) {print $1 "\t" $7}' reference/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt > reference/varid_rsid_map.txt

plink --bfile TEMP4 \
          --update-name reference/varid_rsid_map.txt \
          --make-bed \
          --out TEMP5

awk '$2 !~ /^rs/ {print $2}' TEMP5.bim > reference/missing_rsids.txt

plink --bfile TEMP5 \
          --exclude reference/missing_rsids.txt \
          --make-bed \
          --out TEMP6

plink --bfile TEMP6 \
          --maf 0.01 \
          --make-bed \
          --out genotypes/dosages_processed

rm TEMP*
