#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --mem=256g
#SBATCH --time=20:00:00
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

printf "Training full-sample models and simulating expression for:\n"
printf "CONTEXT_A = ${CONTEXT_A}\n"
printf "CONTEXT_B = ${CONTEXT_B}\n"
printf "GENES = ${GENES}\n"

mkdir ${SLURM_JOB_ID}

awk 'NR > 1 {print $1 "\t" $2}' expression/${CONTEXT_A}.expression.txt > ${SLURM_JOB_ID}/A.individuals.txt
awk 'NR > 1 {print $1 "\t" $2}' expression/${CONTEXT_B}.expression.txt > ${SLURM_JOB_ID}/B.individuals.txt

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
          --pheno expression/${CONTEXT_A}.expression.txt \
          --pheno-name ${ID} \
          --keep ${SLURM_JOB_ID}/A.individuals.txt \
          --make-bed \
          --out ${SLURM_JOB_ID}/${NAME}/A.genotypes

./plink --bfile genotypes/dosages \
          --silent \
          --allow-no-sex \
          --chr ${CHR} \
          --from-bp ${START} \
          --to-bp ${END} \
          --pheno expression/${CONTEXT_B}.expression.txt \
          --pheno-name ${ID} \
          --keep ${SLURM_JOB_ID}/B.individuals.txt \
          --make-bed \
          --out ${SLURM_JOB_ID}/${NAME}/B.genotypes

if ( [ ! -f ${SLURM_JOB_ID}/${NAME}/A.genotypes.bed ] || [ ! -f ${SLURM_JOB_ID}/${NAME}/B.genotypes.bed ] ); then
printf "Unable to extract genotype data for ${ID}. Skipping gene.\n"
printf "Unable to extract genotype data for ${ID}. Skipping gene.\n"
rm -rf ${SLURM_JOB_ID}/${NAME}
continue
fi

Rscript --vanilla simulations/simulate_expression.R ${SLURM_JOB_ID} ${CONTEXT_A} ${CONTEXT_B} ${NAME} ${ID}

rm -rf ${SLURM_JOB_ID}/${NAME}

done < annotations/${GENES}.txt

rm -rf ${SLURM_JOB_ID}

