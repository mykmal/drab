#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64g
#SBATCH --time=24:00:00
#SBATCH --tmp=10g
#SBATCH --partition=amd512,amd2tb,ram256g,ram1t
#SBATCH --mail-type=ALL
#SBATCH --mail-user=scientist@university.edu
#SBATCH -o logs/%j.out
#SBATCH -e logs/%j.err

module load R/4.1.0

export PATH=${PATH}:${DREX}/plink

cd ${DREX}

mkdir temp/${SLURM_JOB_ID}

if [ "$(wc -l < expression/${TISSUE_A}.expression_matrix.txt)" -lt "$(wc -l < expression/${TISSUE_B}.expression_matrix.txt)" ]; then
TISSUE_TEMP=${TISSUE_A}
TISSUE_A=${TISSUE_B}
TISSUE_B=${TISSUE_TEMP}
fi

Rscript scripts/SplitData.R ${TISSUE_A} temp/${SLURM_JOB_ID}

read -r HEADER_A < expression/${TISSUE_A}.expression_matrix.txt
read -r HEADER_B < expression/${TISSUE_B}.expression_matrix.txt

while read -r NAME ID CHR START END; do

if ( [[ ${HEADER_A} != *"${ID}"* ]] || [[ ${HEADER_B} != *"${ID}"* ]] ); then
printf "WARNING: expression data not found for ${ID}. Skipping gene.\n\n"
continue
fi

mkdir temp/${SLURM_JOB_ID}/${NAME}

CHR=$(echo ${CHR} | tr -d "chr")

((START=${START}-500000))
if (( ${START} < 0 )); then
START=0
fi

((END=${END}+500000))

plink --bfile genotypes/dosages \
          --silent \
          --allow-no-sex \
          --chr ${CHR} \
          --from-bp ${START} \
          --to-bp ${END} \
          --pheno temp/${SLURM_JOB_ID}/${TISSUE_A}_part1.expression_matrix.txt \
          --pheno-name ${ID} \
          --keep temp/${SLURM_JOB_ID}/${TISSUE_A}_part1.expression_matrix.txt \
          --make-bed \
          --out temp/${SLURM_JOB_ID}/${NAME}/${TISSUE_A}_part1

plink --bfile genotypes/dosages \
          --silent \
          --allow-no-sex \
          --chr ${CHR} \
          --from-bp ${START} \
          --to-bp ${END} \
          --pheno temp/${SLURM_JOB_ID}/${TISSUE_A}_part2.expression_matrix.txt \
          --pheno-name ${ID} \
          --keep temp/${SLURM_JOB_ID}/${TISSUE_A}_part2.expression_matrix.txt \
          --make-bed \
          --out temp/${SLURM_JOB_ID}/${NAME}/${TISSUE_A}_part2

plink --bfile genotypes/dosages \
          --silent \
          --allow-no-sex \
          --chr ${CHR} \
          --from-bp ${START} \
          --to-bp ${END} \
          --pheno expression/${TISSUE_B}.expression_matrix.txt \
          --pheno-name ${ID} \
          --keep expression/${TISSUE_B}.expression_matrix.txt \
          --make-bed \
          --out temp/${SLURM_JOB_ID}/${NAME}/${TISSUE_B}

if ( [ ! -f temp/${SLURM_JOB_ID}/${NAME}/${TISSUE_A}_part1.bed ] || [ ! -f temp/${SLURM_JOB_ID}/${NAME}/${TISSUE_A}_part2.bed ] || [ ! -f temp/${SLURM_JOB_ID}/${NAME}/${TISSUE_B}.bed ] ); then
printf "Unable to extract genotype data for ${ID}. Skipping gene.\n\n"
rm -rf temp/${SLURM_JOB_ID}/${NAME}
continue
fi

Rscript scripts/DREX.R ${NAME} ${ID} ${TISSUE_A} ${TISSUE_B} ${SLURM_JOB_ID} ${GENES} ${REPLICATES}

rm -rf temp/${SLURM_JOB_ID}/${NAME}

done < annotations/${GENES}.txt

rm -rf temp/${SLURM_JOB_ID}

