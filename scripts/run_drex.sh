#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64g
#SBATCH --time=24:00:00
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

Rscript code/SplitData.R ${TISSUE_A}
Rscript code/SplitData.R ${TISSUE_B}

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

printf "${GENE}\n"

EXPRESSION_A1=$(cat temp/${TISSUE_A}_half1.v8.EUR.normalized_expression.bed | awk -v i=${ID} '$4 == i')
EXPRESSION_B1=$(cat temp/${TISSUE_B}_half1.v8.EUR.normalized_expression.bed | awk -v i=${ID} '$4 == i')
EXPRESSION_A2=$(cat temp/${TISSUE_A}_half2.v8.EUR.normalized_expression.bed | awk -v i=${ID} '$4 == i')
EXPRESSION_B2=$(cat temp/${TISSUE_B}_half2.v8.EUR.normalized_expression.bed | awk -v i=${ID} '$4 == i')
if ([ -z "${EXPRESSION_A1}" ] || [ -z "${EXPRESSION_B1}" ] || [ -z "${EXPRESSION_A2}" ] || [ -z "${EXPRESSION_B2}" ]); then
printf "Expression data not found. Skipping gene.\n"
continue
fi

OUT_A="temp/${TISSUE_A}.${BATCH}/${TISSUE_A}.${ID}"
OUT_B="temp/${TISSUE_B}.${BATCH}/${TISSUE_B}.${ID}"

echo ${EXPRESSION_A1} | tr ' ' '\n' | tail -n+5 | paste temp/${TISSUE_A}_half1.v8.EUR.covariates.txt.HEADER - | awk '{print 0 "\t" $0}' > ${OUT_A}_half1.pheno
echo ${EXPRESSION_B1} | tr ' ' '\n' | tail -n+5 | paste temp/${TISSUE_B}_half1.v8.EUR.covariates.txt.HEADER - | awk '{print 0 "\t" $0}' > ${OUT_B}_half1.pheno
echo ${EXPRESSION_A2} | tr ' ' '\n' | tail -n+5 | paste temp/${TISSUE_A}_half2.v8.EUR.covariates.txt.HEADER - | awk '{print 0 "\t" $0}' > ${OUT_A}_half2.pheno
echo ${EXPRESSION_B2} | tr ' ' '\n' | tail -n+5 | paste temp/${TISSUE_B}_half2.v8.EUR.covariates.txt.HEADER - | awk '{print 0 "\t" $0}' > ${OUT_B}_half2.pheno

plink --silent --bfile genotypes/dosages_processed --allow-no-sex --chr ${CHR} --from-bp ${START} --to-bp ${END} --pheno ${OUT_A}_half1.pheno --keep ${OUT_A}_half1.pheno --make-bed --out ${OUT_A}_half1
plink --silent --bfile genotypes/dosages_processed --allow-no-sex --chr ${CHR} --from-bp ${START} --to-bp ${END} --pheno ${OUT_B}_half1.pheno --keep ${OUT_B}_half1.pheno --make-bed --out ${OUT_B}_half1
plink --silent --bfile genotypes/dosages_processed --allow-no-sex --chr ${CHR} --from-bp ${START} --to-bp ${END} --pheno ${OUT_A}_half2.pheno --keep ${OUT_A}_half2.pheno --make-bed --out ${OUT_A}_half2
plink --silent --bfile genotypes/dosages_processed --allow-no-sex --chr ${CHR} --from-bp ${START} --to-bp ${END} --pheno ${OUT_B}_half2.pheno --keep ${OUT_B}_half2.pheno --make-bed --out ${OUT_B}_half2
if ([ ! -f ${OUT_A}_half1.bed ] || [ ! -f ${OUT_B}_half1.bed ] || [ ! -f ${OUT_A}_half2.bed ] || [ ! -f ${OUT_B}_half2.bed ]); then
printf "Unable to extract genotype data. Skipping gene.\n"
rm -f ${OUT_A}*
rm -f ${OUT_B}*
continue
fi

Rscript code/DREX.R ${NAME} ${ID} ${TISSUE_A} ${OUT_A} ${TISSUE_B} ${OUT_B} ${BATCH}

rm -f ${OUT_A}*
rm -f ${OUT_B}*

done

rm -fr temp/${TISSUE_A}.${BATCH}
rm -fr temp/${TISSUE_B}.${BATCH}

rm -f temp/${TISSUE_A}*
rm -f temp/${TISSUE_B}*

