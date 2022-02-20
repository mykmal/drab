#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64g
#SBATCH --time=1:00:00
#SBATCH --tmp=100g
#SBATCH -p amdlarge,amd512,amd2tb,ram256g,ram1t
#SBATCH --mail-type=ALL
#SBATCH --mail-user=scientist@university.edu
#SBATCH -o logs/%j.out
#SBATCH -e logs/%j.err

DREX="/path/to/drex"

module load R/4.1.0

export PATH=${PATH}:${DREX}/plink

cd ${DREX}

for f in expression_covariates/*.txt
do
cat ${f} | head -n1 | tr "\t" "\n" | tail -n+2 > ${f}.HEADER
done

Rscript code/ProcessCovariates.R

awk '(NR > 6) && ($1 ~ /^chr[1-22]/) && ($3 == "gene") {print $16 "\t" $10 "\t" $1 "\t" $4 "\t" $5 "\t" $14}' reference/gencode.v26.GRCh38.genes.gtf | \
tr -d ";\"" > reference/gene_annotation.txt

awk '$6 == "protein_coding" {print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5}' reference/gene_annotation.txt > reference/pc_gene_annotation.txt

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
          --geno 0 \
          --make-bed \
          --out TEMP3

plink --bfile TEMP3 \
          --autosome \
          --snps-only just-acgt \
          --make-bed \
          --out TEMP4

awk -F '\t' '/G_C|C_G|A_T|T_A/ {print $2}' TEMP4.bim > reference/ambiguous_snps.txt

plink --bfile TEMP4 \
          --exclude reference/ambiguous_snps.txt \
          --make-bed \
          --out TEMP5

awk -F '\t' '(NR>1) && ($7 ~ /^rs/) {print $1 "\t" $7}' reference/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt > reference/varid_rsid_map.txt

plink --bfile TEMP5 \
          --update-name reference/varid_rsid_map.txt \
          --make-bed \
          --out TEMP6

awk '$2 !~ /^rs/ {print $2}' TEMP6.bim > reference/missing_rsids.txt

plink --bfile TEMP6 \
          --exclude reference/missing_rsids.txt \
          --make-bed \
          --out TEMP7

plink --bfile TEMP7 \
          --hwe 1e-6 midp include-nonctrl \
          --make-bed \
          --out TEMP8

plink --bfile TEMP8 \
          --maf 0.01 \
          --make-bed \
          --out genotypes/dosages_processed

rm TEMP*

