#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64g
#SBATCH --time=7:00:00
#SBATCH --tmp=10g
#SBATCH -p amdlarge,amd512,amd2tb,ram256g,ram1t
#SBATCH --mail-type=ALL
#SBATCH --mail-user=scientist@university.edu
#SBATCH -o %A_%a.out
#SBATCH -e %A_%a.err

export DREX="/path/to/drex"

module load R/4.1.0

export PATH=${PATH}:${DREX}/plink

cd ${DREX}

TISSUE=Muscle_Skeletal

mkdir -p weights/${TISSUE}

BATCH=${SLURM_ARRAY_TASK_ID}

mkdir temp/${TISSUE}.${BATCH}

cat reference/pc_gene_annotation.txt | awk -v i=${BATCH} 'NR > (i-1)*1000 && NR <= i*1000' | while read GENE; do

NAME=$(echo ${GENE} | awk '{print $1}')
ID=$(echo ${GENE} | awk '{print $2}')
CHR=$(echo ${GENE} | awk '{print $3}' | tr -d "chr")
START=$(echo ${GENE} | awk '{p=$4 - 500e3; if(p<0) p=0; print p;}')
END=$(echo ${GENE} | awk '{print $5 + 500e3}')

PARAM=$(zcat expression_matrices/${TISSUE}.v8.EUR.normalized_expression.bed.gz | awk -v i=${ID} '$4 == i')
if [ -z "${PARAM}" ]; then
continue
fi

printf "\n${GENE}\n"

OUT="temp/${TISSUE}.${BATCH}/${TISSUE}.${ID}"

echo ${PARAM} | tr ' ' '\n' | tail -n+5 | paste expression_covariates/${TISSUE}.v8.EUR.covariates.txt.HEADER - | awk '{print 0 "\t" $0}' > ${OUT}.pheno

plink --silent --bfile genotypes/dosages_processed --allow-no-sex --chr ${CHR} --from-bp ${START} --to-bp ${END} --pheno $OUT.pheno --keep $OUT.pheno --make-bed --out ${OUT}
if [ ! -f ${OUT}.bed ]; then
continue
fi

FINAL_OUT="weights/${TISSUE}/${TISSUE}.${ID}"

Rscript ./fusion_twas/FUSION.compute_weights.R \
--bfile ${OUT} --tmp ${OUT}.tmp --out ${FINAL_OUT} --PATH_gcta fusion_twas/gcta_nr_robust --models enet --covar expression_covariates/${TISSUE}.v8.EUR.covariates.plink.txt --hsq_p 1.0 --verbose 1

echo ${ID} >> "weights/${TISSUE}/gene_list.txt

rm -f ${OUT}.*

done

rm -fr temp/${TISSUE}.${BATCH}
