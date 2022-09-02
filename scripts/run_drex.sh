#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --mem=64g
#SBATCH --time=50:00:00
#SBATCH --tmp=100g
#SBATCH --partition=agsmall,aglarge,ag2tb
#SBATCH --mail-type=ALL
#SBATCH --mail-user=malak039@umn.edu
#SBATCH -o logs/%j.out
#SBATCH -e logs/%j.err

module load R/4.1.0

cd ${DREX}

mkdir ${SLURM_JOB_ID}

Rscript --vanilla scripts/split_data.R ${TISSUE_A} ${TISSUE_B} ${SLURM_JOB_ID}

read -r HEADER_A < expression/${TISSUE_A}.expression_matrix.txt
read -r HEADER_B < expression/${TISSUE_B}.expression_matrix.txt

while read -r NAME ID CHR START END; do

if ( [[ ${HEADER_A} != *"${ID}"* ]] || [[ ${HEADER_B} != *"${ID}"* ]] ); then
printf "WARNING: expression data not found for ${ID}. Skipping gene.\n"
continue
fi

mkdir ${SLURM_JOB_ID}/${NAME}

CHR=$(echo ${CHR} | tr -d "chr")

((START=${START}-500000))
if (( ${START} < 0 )); then
START=0
fi

((END=${END}+500000))

awk 'NR > 1 {print $1 "\t" $2}' ${SLURM_JOB_ID}/part1.expression_matrix.txt > ${SLURM_JOB_ID}/part1.individuals.txt
awk 'NR > 1 {print $1 "\t" $2}' ${SLURM_JOB_ID}/part2.expression_matrix.txt > ${SLURM_JOB_ID}/part2.individuals.txt
awk 'NR > 1 {print $1 "\t" $2}' ${SLURM_JOB_ID}/part3.expression_matrix.txt > ${SLURM_JOB_ID}/part3.individuals.txt

./plink --bfile genotypes/dosages \
          --silent \
          --allow-no-sex \
          --chr ${CHR} \
          --from-bp ${START} \
          --to-bp ${END} \
          --pheno ${SLURM_JOB_ID}/part1.expression_matrix.txt \
          --pheno-name ${ID} \
          --keep ${SLURM_JOB_ID}/part1.individuals.txt \
          --make-bed \
          --out ${SLURM_JOB_ID}/${NAME}/part1

./plink --bfile genotypes/dosages \
          --silent \
          --allow-no-sex \
          --chr ${CHR} \
          --from-bp ${START} \
          --to-bp ${END} \
          --pheno ${SLURM_JOB_ID}/part2.expression_matrix.txt \
          --pheno-name ${ID} \
          --keep ${SLURM_JOB_ID}/part2.individuals.txt \
          --make-bed \
          --out ${SLURM_JOB_ID}/${NAME}/part2

./plink --bfile genotypes/dosages \
          --silent \
          --allow-no-sex \
          --chr ${CHR} \
          --from-bp ${START} \
          --to-bp ${END} \
          --pheno ${SLURM_JOB_ID}/part3.expression_matrix.txt \
          --pheno-name ${ID} \
          --keep ${SLURM_JOB_ID}/part3.individuals.txt \
          --make-bed \
          --out ${SLURM_JOB_ID}/${NAME}/part3

if ( [ ! -f ${SLURM_JOB_ID}/${NAME}/part1.bed ] || [ ! -f ${SLURM_JOB_ID}/${NAME}/part2.bed ] || [ ! -f ${SLURM_JOB_ID}/${NAME}/part3.bed ] ); then
printf "Unable to extract genotype data for ${ID}. Skipping gene.\n"
rm -rf ${SLURM_JOB_ID}/${NAME}
continue
fi

Rscript --vanilla scripts/DREX.R ${SLURM_JOB_ID} ${TISSUE_A} ${TISSUE_B} ${GENES} ${NAME} ${ID}

rm -rf ${SLURM_JOB_ID}/${NAME}

done < annotations/${GENES}.txt

rm -rf ${SLURM_JOB_ID}

