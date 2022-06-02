#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64g
#SBATCH --time=2:00:00
#SBATCH --tmp=10g
#SBATCH --partition=agsmall,ag2tb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=scientist@university.edu
#SBATCH -o logs/%j.out
#SBATCH -e logs/%j.err

module load R/4.1.0

export PATH=${PATH}:${DREX}/plink

cd ${DREX}

awk '(NR > 6) && ($1 ~ /^chr[1-22]/) && ($3 == "gene") {print $16 "\t" $10 "\t" $1 "\t" $4 "\t" $5 "\t" $14}' raw/gencode.v26.GRCh38.genes.gtf | \
tr -d ";\"" > annotations/all_genes.txt

awk '$6 == "protein_coding" {print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5}' annotations/all_genes.txt > annotations/protein_coding_genes.txt

plink --vcf raw/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.vcf.gz \
          --const-fid \
          --vcf-require-gt \
          --biallelic-only strict \
          --make-bed \
          --out TEMP1

Rscript --vanilla scripts/ExtractEUR.R

plink --bfile TEMP1 \
          --keep EUR_samples.txt \
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

awk -F '\t' '/G_C|C_G|A_T|T_A/ {print $2}' TEMP4.bim > ambiguous_snps.txt

plink --bfile TEMP4 \
          --exclude ambiguous_snps.txt \
          --make-bed \
          --out TEMP5

awk -F '\t' '(NR > 1) && ($7 ~ /^rs/) {print $1 "\t" $7}' raw/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt > varid_rsid_map.txt

plink --bfile TEMP5 \
          --update-name varid_rsid_map.txt \
          --make-bed \
          --out TEMP6

awk '$2 !~ /^rs/ {print $2}' TEMP6.bim > missing_rsids.txt

plink --bfile TEMP6 \
          --exclude missing_rsids.txt \
          --make-bed \
          --out TEMP7

plink --bfile TEMP7 \
          --hwe 1e-6 midp include-nonctrl \
          --make-bed \
          --out TEMP8

plink --bfile TEMP8 \
          --maf 0.01 \
          --make-bed \
          --out genotypes/dosages

rm TEMP*
rm EUR_samples.txt
rm ambiguous_snps.txt
rm varid_rsid_map.txt
rm missing_rsids.txt

Rscript --vanilla scripts/GTEx2plink.R

