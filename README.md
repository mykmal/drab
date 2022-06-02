# DREX: Differential Regulation of EXpression

DREX is an R package for identifying genes with tissue-specific genetic regulation of expression. The documentation below covers DREX installation, required input data, and usage. We also provide an appendix with information on how to obtain the data used in our paper.

**Note:** The DREX workflow is designed to be run on a linux HPC system and all provided commands are for bash (with the exception of a few R commands, which are prefaced by the `>` symbol).

## Setup

* First, clone the DREX repository and create the required folder structure.
```
git clone https://github.com/MykMal/drex.git
cd drex
mkdir annotations covariates expression genotypes logs output plink temp
```
* Launch R and install the packages BEDMatrix and glmnet. We used R v4.1.0 x86_64, BEDMatrix 2.0.3, and glmnet 4.1.4.
```
> install.packages(c("BEDMatrix", "glmnet"))
```
* Download plink to the `drex/plink` folder. We used plink v1.90b6.26 64-bit (2 Apr 2022).
```
cd plink
wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20220402.zip
unzip plink_linux_x86_64_20220402.zip
rm plink_linux_x86_64_20220402.zip
```
* Each of the `*.sh` files in `drex/scripts` is a shell script prefaced by SLURM commands. Modify the `#SBATCH` commands at the beginning of each SLURM script as appropriate for your HPC cluster. In particular, be sure to set the `--mail-user` flag to your own email address and the `--partition` flag to the list of SLURM partitions on your cluster. The remaining commands may be left at their defaults.

## Input data formats

The subsections below explain the files and file formats that DREX expects. All of the files mentioned below are required.

### Gene annotations

DREX requires a gene annotation file listing all of the genes on which it should run. This has to be a plain-text, tab-delimited file without a header line. Each line should contain information for a single gene, with the following five fields:

1. Gene name
2. Gene ID (e.g. from ENSEMBL)
3. Chromosome (with or without the chr prefix)
4. Start position (in base pairs)
5. End position (in base pairs)

For example, the lines for the first three guanylate binding protein genes would be
```
GBP1	ENSG00000117228.9	chr1	89052319	89065360
GBP2	ENSG00000162645.12	chr1	89106132	89126114
GBP3	ENSG00000117226.11	chr1	89006666	89022894
```
Save your gene annotation files with file names of your choice in `drex/annotations`.

### Genotype data

DREX requires individual-level whole genome sequencing data in plink bed/bim/fam format. After performing all desired quality control, save your fully-processed genotype data as `dosages.bed`, `dosages.bim`, and `dosages.fam` in `drex/genotypes`.

### Gene expression data

DREX requires individual-level gene expression data for each tissue of interest, and it is assumed that the RNA-Seq values have already been fully processed and normalized. Expression data should be in tissue-specific, plain-text, tab-delimited files that begin with a header line. Each line should contain information for a single individual with family ID in the first field, within-family ID in the second field, and per-gene expression levels in the remaining fields. For example, the first three lines might be
```
FID	IID	ENSG00000117228.9	ENSG00000162645.12	ENSG00000117226.11
0	indivA	-0.083051769477432	0.808844404113396	1.31169125330214
0	indivB	0.00672624465727554	-1.09866518781071	0.350055616620479
```
Use the naming convention `<tissue>.expression_matrix.txt` and save all of the gene expression files in `drex/expression`.

**Note:** The FIDs and IIDs must be consistent with those used for the genotype data, while the gene IDs in the header line must be consistent with those used in the gene annotation files.

### Expression covariates

The format for expression covariates is analogous to the format for gene expression described above. Covariates should be in tissue-specific, plain-text, tab-delimited files that begin with a header line. Each line should contain information for a single individual with family ID in the first field, within-family ID in the second field, and covariates in the remaining fields. For example, the first three lines might be
```
FID	IID	PC1	PC2	InferredCov	pcr
0	indivA	0.0147	-0.0072	0.0262378174811602	1
0	indivB	0.0161	0.0037	-0.0514548756182194	1
```
Use the naming convention `<tissue>.expression_covariates.txt` and save all of the covariate files in `drex/covariates`.

**Note:** The FIDs and IIDs must be consistent with those used in the gene expression files.

## Running DREX

To run DREX, submit the `scripts/run_drex.sh` shell script as a SLURM job with the appropriate flags as shown below. For example, to test whether the expression of genes listed in the annotation file `all_genes.txt` is differentially regulated in tissues labeled as `Whole_Blood` and `Brain_Cortex`, run the command
```
sbatch --export=TISSUE_A="Whole_Blood",TISSUE_B="Brain_Cortex",GENES="all_genes",PERMUTATIONS="1000",DREX=$(pwd) scripts/run_drex.sh
```
In practice, replace `Whole_Blood` and `Brain_Cortex` with the names of your desired tissues and `all_genes` with the name of your annotation file. If desired, you can also change the default value of 1,000 permutations for approximating the null distribution.

The results will be saved to `output/Whole_Blood-Brain_Cortex-all_genes.txt`. (Here `Whole_Blood`, `Brain_Cortex`, and `all_genes` will be replaced with the tissue names and annotation file name you specified when running DREX.) This is a tab-delimited, plain-text file without a header line. Each line contains information for a single gene, with the following fields:

1. Gene name
2. Gene ID (e.g. from ENSEMBL)
3. The likelihood-ratio test p-value
4. The distribution-free test p-value

If the p-values for a given gene are significant, then we conclude that the genetic regulation of that gene's expression is significantly different between the two tissues.

**Important:** The reported p-values are not corrected for multiple testing. A Bonferroni correction or some other appropriate family-wise error rate method should be applied before drawing conclusions.

## Appendix: download and prepare GTEx data

First, create the folder `drex/raw` to store the unprocessed GTEx data sets. This folder may be safely deleted after completing all of the steps in this section.

From the **GTEx Analysis V8 (dbGaP Accession phs000424.v8.p2)** section of https://www.gtexportal.org/home/datasets, download the following files:

* `GTEx_Analysis_v8_eQTL_EUR.tar` (under the sub-heading "Single-Tissue cis-QTL Data")  
Unpack this archive and move the folders `expression_matrices` and `expression_covariates` to `drex/raw`. (The other folder in the tar is not needed.)
* `GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt.gz` (under the sub-heading "Reference")  
Uncompress this file and move it to `drex/raw`.
* `gencode.v26.GRCh38.genes.gtf` (under the sub-heading "Reference")  
Move this file to `drex/raw`.

After obtaining access to the GTEx data in [dbGaP](https://www.ncbi.nlm.nih.gov/gap/) (accession phs000424.v8.p2), follow the dbGaP documentation to download the following files:

* `phg001219.v1.GTEx_v8_WGS.genotype-calls-vcf.c1.GRU.tar`  
Extract the file `GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.vcf.gz` from the tar and move it to `drex/raw`. (The other files in the tar are not needed.)
* `phs000424.v8.pht002742.v8.p2.c1.GTEx_Subject_Phenotypes.GRU.txt.gz`  
Uncompress this file and move it to `drex/raw`.

To prepare the GTEx data for use with DREX, from your main `drex` folder run
```
sbatch --export=DREX=$(pwd) scripts/prepare_data.sh
```
This script will create an annotation file with all GTEx genes and another one with only protein-coding genes, perform standard quality control steps on the genotype data, and reformat the gene expression data and expression covariates.

