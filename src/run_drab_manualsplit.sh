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

cd ${DRAB}

# Delete (or comment out) the following line if not using SLURM:
module load R/4.2.2-openblas

# Uncomment the following two lines if not using SLURM:
#SLURM_JOB_ID=$((1 + $RANDOM % 1000))
#exec 1> logs/${SLURM_JOB_ID}.out 2> logs/${SLURM_JOB_ID}.err

printf "DRAB runtime parameters:\n"
printf "D_A PATH = ${DA_PATH}\n"
printf "D_B PATH = ${DB_PATH}\n"
printf "D_T PATH = ${DT_PATH}\n"
printf "GENES = ${GENES}\n"
printf "BOOT = ${BOOT}\n\n"

mkdir ${SLURM_JOB_ID}

DA_EXPRESSION=$(basename "${DA_PATH}")
DB_EXPRESSION=$(basename "${DB_PATH}")
DT_EXPRESSION=$(basename "${DT_PATH}")

DA_SIZE=$(wc -l < ${DA_PATH}.expression.txt)
if (( ${DA_SIZE} < 10 )); then
printf "Error: insufficient data. The training sample size for ${DA_EXPRESSION} is ${DA_SIZE}.\n"
rm -rf ${SLURM_JOB_ID}
exit 1
fi

DB_SIZE=$(wc -l < ${DB_PATH}.expression.txt)
if (( ${DB_SIZE} < 10 )); then
printf "Error: insufficient data. The training sample size for ${DB_EXPRESSION} is ${DB_SIZE}.\n"
rm -rf ${SLURM_JOB_ID}
exit 1
fi

DT_SIZE=$(wc -l < ${DT_PATH}.expression.txt)
if (( ${DT_SIZE} < 10 )); then
printf "Error: insufficient data. The testing sample size is ${DT_SIZE}.\n"
rm -rf ${SLURM_JOB_ID}
exit 1
fi

awk 'NR > 1 {print $1 "\t" $2}' ${DA_PATH}.expression.txt > ${SLURM_JOB_ID}/da_individuals.txt
awk 'NR > 1 {print $1 "\t" $2}' ${DB_PATH}.expression.txt > ${SLURM_JOB_ID}/db_individuals.txt
awk 'NR > 1 {print $1 "\t" $2}' ${DT_PATH}.expression.txt > ${SLURM_JOB_ID}/dt_individuals.txt

read -r HEADER_A < ${DA_PATH}.expression.txt
read -r HEADER_B < ${DB_PATH}.expression.txt
read -r HEADER_T < ${DT_PATH}.expression.txt

cp ${DA_PATH}.covariates.txt ${SLURM_JOB_ID}/da.covariates.txt
cp ${DB_PATH}.covariates.txt ${SLURM_JOB_ID}/db.covariates.txt
cp ${DT_PATH}.covariates.txt ${SLURM_JOB_ID}/dt.covariates.txt

while read -r NAME ID CHR START END ETC; do

if ( [[ ${HEADER_A} != *"${ID}"* ]] || [[ ${HEADER_B} != *"${ID}"* ]] || [[ ${HEADER_T} != *"${ID}"* ]] ); then
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
          --pheno ${DA_PATH}.expression.txt \
          --pheno-name ${ID} \
          --keep ${SLURM_JOB_ID}/da_individuals.txt \
          --make-bed \
          --out ${SLURM_JOB_ID}/${NAME}/da

./plink --bfile genotypes/dosages \
          --silent \
          --allow-no-sex \
          --chr ${CHR} \
          --from-bp ${START} \
          --to-bp ${END} \
          --pheno ${DB_PATH}.expression.txt \
          --pheno-name ${ID} \
          --keep ${SLURM_JOB_ID}/db_individuals.txt \
          --make-bed \
          --out ${SLURM_JOB_ID}/${NAME}/db

./plink --bfile genotypes/dosages \
          --silent \
          --allow-no-sex \
          --chr ${CHR} \
          --from-bp ${START} \
          --to-bp ${END} \
          --pheno ${DT_PATH}.expression.txt \
          --pheno-name ${ID} \
          --keep ${SLURM_JOB_ID}/dt_individuals.txt \
          --make-bed \
          --out ${SLURM_JOB_ID}/${NAME}/dt

if ( [ ! -f ${SLURM_JOB_ID}/${NAME}/da.bed ] || [ ! -f ${SLURM_JOB_ID}/${NAME}/db.bed ] || [ ! -f ${SLURM_JOB_ID}/${NAME}/dt.bed ] ); then
printf "Unable to extract genotype data for ${ID}. Skipping gene.\n"
rm -rf ${SLURM_JOB_ID}/${NAME}
continue
fi

Rscript --vanilla src/drab.R ${SLURM_JOB_ID} ${DA_EXPRESSION} ${DB_EXPRESSION} ${DT_EXPRESSION} ${BOOT} ${NAME} ${ID}

rm -rf ${SLURM_JOB_ID}/${NAME}

done < annotations/${GENES}.txt

rm -rf ${SLURM_JOB_ID}

