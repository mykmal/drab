#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64g
#SBATCH --time=10:00:00
#SBATCH --tmp=10g
#SBATCH -p amdlarge,amd512,amd2tb,ram256g,ram1t
#SBATCH --mail-type=ALL
#SBATCH --mail-user=scientist@university.edu
#SBATCH -o logs/%A_%a.out
#SBATCH -e logs/%A_%a.err

DREX="/path/to/drex"

TISSUE_A=Muscle_Skeletal
TISSUE_B=Brain_Cortex
PC_ONLY=TRUE

module load R/4.1.0

export PATH=${PATH}:${DREX}/plink

cd ${DREX}

BATCH=${SLURM_ARRAY_TASK_ID}

mkdir temp/${TISSUE_A}.${BATCH}
mkdir temp/${TISSUE_B}.${BATCH}

if [ ${PC_ONLY} = "FALSE" ] || [ ${PC_ONLY} = "F" ]; then
FILE="gene_annotation.txt"
else
FILE="pc_gene_annotation.txt"
fi

cat reference/${FILE} | awk -v i=${BATCH} 'NR > (i-1)*1000 && NR <= i*1000' | while read GENE; do

NAME=$(echo ${GENE} | awk '{print $1}')
ID=$(echo ${GENE} | awk '{print $2}')
CHR=$(echo ${GENE} | awk '{print $3}' | tr -d "chr")
START=$(echo ${GENE} | awk '{p=$4 - 500e3; if(p<0) p=0; print p;}')
END=$(echo ${GENE} | awk '{print $5 + 500e3}')

EXPRESSION_A=$(zcat expression_matrices/${TISSUE_A}.v8.EUR.normalized_expression.bed.gz | awk -v i=${ID} '$4 == i')
EXPRESSION_B=$(zcat expression_matrices/${TISSUE_B}.v8.EUR.normalized_expression.bed.gz | awk -v i=${ID} '$4 == i')
if [ -z "${EXPRESSION_A}" ] || [ -z "${EXPRESSION_B}" ]; then
continue
fi

OUT_A="temp/${TISSUE_A}.${BATCH}/${TISSUE_A}.${ID}"
OUT_B="temp/${TISSUE_B}.${BATCH}/${TISSUE_B}.${ID}"

echo ${EXPRESSION_A} | tr ' ' '\n' | tail -n+5 | paste expression_covariates/${TISSUE_A}.v8.EUR.covariates.txt.HEADER - | awk '{print 0 "\t" $0}' > ${OUT_A}.pheno
echo ${EXPRESSION_B} | tr ' ' '\n' | tail -n+5 | paste expression_covariates/${TISSUE_B}.v8.EUR.covariates.txt.HEADER - | awk '{print 0 "\t" $0}' > ${OUT_B}.pheno

plink --silent --bfile genotypes/dosages_processed --allow-no-sex --chr ${CHR} --from-bp ${START} --to-bp ${END} --pheno ${OUT_A}.pheno --keep ${OUT_A}.pheno --make-bed --out ${OUT_A}
plink --silent --bfile genotypes/dosages_processed --allow-no-sex --chr ${CHR} --from-bp ${START} --to-bp ${END} --pheno ${OUT_B}.pheno --keep ${OUT_B}.pheno --make-bed --out ${OUT_B}
if [ ! -f ${OUT_A}.bed ] || [ ! -f ${OUT_B}.bed ]; then
rm -f ${OUT_A}*
rm -f ${OUT_B}*
continue
fi

printf "${GENE}\n"

Rscript code/DREX.R ${NAME} ${ID} ${TISSUE_A} ${OUT_A} ${TISSUE_B} ${OUT_B} ${BATCH}

rm -f ${OUT_A}*
rm -f ${OUT_B}*

done

rm -fr temp/${TISSUE_A}.${BATCH}
rm -fr temp/${TISSUE_B}.${BATCH}

