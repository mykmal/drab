#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --mem=64g
#SBATCH --time=90:00:00
#SBATCH --tmp=10g
#SBATCH --partition=msismall,msilarge,msibigmem
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=malak039@umn.edu
#SBATCH -o logs/%j.out
#SBATCH -e logs/%j.err

module load R/4.1.0

cd ${DRAB}

printf "DRAB runtime parameters:\n"
printf "CONTEXT_A = ${CONTEXT_A}\n"
printf "CONTEXT_B = ${CONTEXT_B}\n"
printf "GENES = ${GENES}\n"
printf "BOOT = ${BOOT}\n\n"

mkdir ${SLURM_JOB_ID}

Rscript --vanilla src/split_data.R ${CONTEXT_A} ${CONTEXT_B} ${SLURM_JOB_ID}

DA_SIZE=$(wc -l < ${SLURM_JOB_ID}/part1.expression.txt)
if (( ${DA_SIZE} < 10 )); then
printf "Error: insufficient data. The training sample size for ${CONTEXT_A} is ${DA_SIZE}.\n"
rm -rf ${SLURM_JOB_ID}
exit 1
fi

DB_SIZE=$(wc -l < ${SLURM_JOB_ID}/part2.expression.txt)
if (( ${DB_SIZE} < 10 )); then
printf "Error: insufficient data. The training sample size for ${CONTEXT_B} is ${DB_SIZE}.\n"
rm -rf ${SLURM_JOB_ID}
exit 1
fi

DT_SIZE=$(wc -l < ${SLURM_JOB_ID}/part3.expression.txt)
if (( ${DT_SIZE} < 10 )); then
printf "Error: insufficient data. The testing sample size is ${DT_SIZE}.\n"
rm -rf ${SLURM_JOB_ID}
exit 1
fi

awk 'NR > 1 {print $1 "\t" $2}' ${SLURM_JOB_ID}/part1.expression.txt > ${SLURM_JOB_ID}/part1.individuals.txt
awk 'NR > 1 {print $1 "\t" $2}' ${SLURM_JOB_ID}/part2.expression.txt > ${SLURM_JOB_ID}/part2.individuals.txt
awk 'NR > 1 {print $1 "\t" $2}' ${SLURM_JOB_ID}/part3.expression.txt > ${SLURM_JOB_ID}/part3.individuals.txt

read -r HEADER_A < expression/${CONTEXT_A}.expression.txt
read -r HEADER_B < expression/${CONTEXT_B}.expression.txt

while read -r NAME ID CHR START END ETC; do

if ( [[ ${HEADER_A} != *"${ID}"* ]] || [[ ${HEADER_B} != *"${ID}"* ]] ); then
printf "WARNING: expression data not found for ${ID}. Skipping gene.\n"
continue
fi

mkdir ${SLURM_JOB_ID}/${NAME}

CHR=$(echo ${CHR} | tr -d "chr")

START=$(( ${START} - 500000 ))
if (( ${START} < 0 )); then
START=0
fi

END=$(( ${END} + 500000 ))

./plink --bfile genotypes/dosages \
          --silent \
          --allow-no-sex \
          --chr ${CHR} \
          --from-bp ${START} \
          --to-bp ${END} \
          --pheno ${SLURM_JOB_ID}/part1.expression.txt \
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
          --pheno ${SLURM_JOB_ID}/part2.expression.txt \
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
          --pheno ${SLURM_JOB_ID}/part3.expression.txt \
          --pheno-name ${ID} \
          --keep ${SLURM_JOB_ID}/part3.individuals.txt \
          --make-bed \
          --out ${SLURM_JOB_ID}/${NAME}/part3

if ( [ ! -f ${SLURM_JOB_ID}/${NAME}/part1.bed ] || [ ! -f ${SLURM_JOB_ID}/${NAME}/part2.bed ] || [ ! -f ${SLURM_JOB_ID}/${NAME}/part3.bed ] ); then
printf "Unable to extract genotype data for ${ID}. Skipping gene.\n"
rm -rf ${SLURM_JOB_ID}/${NAME}
continue
fi

Rscript --vanilla src/drab.R ${SLURM_JOB_ID} ${CONTEXT_A} ${CONTEXT_B} ${GENES} ${BOOT} ${NAME} ${ID}

rm -rf ${SLURM_JOB_ID}/${NAME}

done < annotations/${GENES}.txt

rm -rf ${SLURM_JOB_ID}

