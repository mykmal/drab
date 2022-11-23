# DRAB: Differential Regulation Analysis by Bootstrapping

DRAB is a tool for identifying genes with significantly different patterns of local genetic regulation between two tissues or other biological contexts. The documentation below covers DRAB installation, required input data, and usage. For information on the DRAB methodology, see our paper "Identifying genes with tissue-specific patterns of genetic regulation" by Malakhov et al.

**Note:** The DRAB pipeline is designed to be run on a Linux cluster through SLURM. No other platforms are currently supported.

## Setup

* First, clone the DRAB repository and create the required folder structure.
```
git clone https://github.com/MykMal/drab.git
cd drab
mkdir annotations covariates expression genotypes logs output
```
* Launch R and install the packages BEDMatrix and glmnet. We used R v4.1.0 x86_64, BEDMatrix 2.0.3, and glmnet 4.1.4.
```
install.packages(c("BEDMatrix", "glmnet"))
```
* Download PLINK to the main `drab` folder. We used PLINK v1.90b6.26 64-bit (2 Apr 2022).
```
wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20220402.zip
unzip plink_linux_x86_64_20220402.zip
rm plink_linux_x86_64_20220402.zip prettify toy*
```
* The files `drab/src/run_drex.sh` and `drab/util/prepare_data.sh` are shell scripts prefaced by SLURM commands. Modify the `#SBATCH` commands at the beginning of both SLURM scripts as appropriate for your compute cluster. In particular, be sure to set the `--partition` flag to the list of SLURM partitions on your cluster and the `--mail-user` flag to your own email address. The remaining commands may be left at their defaults.

## Input data formats

The subsections below explain the files and file formats that DRAB expects. All of the files mentioned here are required.

### Gene annotations

DRAB requires a gene annotation file listing all of the genes on which it should run. This has to be a plain-text, tab-delimited file without a header line. Each line should contain information for a single gene, with the following five fields:

1. Gene name
2. Gene ID (e.g. from ENSEMBL)
3. Chromosome (with or without the chr prefix)
4. Start position (in base pairs)
5. End position (in base pairs)

For example, the lines for the first three guanylate binding protein genes might be
```
GBP1	ENSG00000117228.9	chr1	89052319	89065360
GBP2	ENSG00000162645.12	chr1	89106132	89126114
GBP3	ENSG00000117226.11	chr1	89006666	89022894
```
Additional fields are allowed but will be ignored. Save your gene annotation files with file names of your choice in `drab/annotations`.

### Genotype data

DRAB requires individual-level whole genome sequencing data in PLINK bed/bim/fam format. After performing all desired quality control, save your fully processed genotype data as `dosages.bed`, `dosages.bim`, and `dosages.fam` in `drab/genotypes`.

### Gene expression data

DRAB requires individual-level gene expression data for each tissue or context of interest, and it is assumed that the RNA-Seq values have already been fully processed and normalized. Expression data should be in context-specific, plain-text, tab-delimited files that begin with a header line. Each line after the header should contain information for a single individual with family ID in the first field, within-family ID in the second field, and per-gene expression levels in the remaining fields. For example, the first three lines might be
```
FID	IID	ENSG00000117228.9	ENSG00000162645.12	ENSG00000117226.11
0	indivA	-0.083051769477432	0.808844404113396	1.31169125330214
0	indivB	0.00672624465727554	-1.09866518781071	0.350055616620479
```
Use the naming convention `<context>.expression.txt` and save all of the gene expression files in `drab/expression`.

**Note:** The FIDs and IIDs must be consistent with those used for the genotype data, while the gene IDs in the header line must be consistent with those used in the gene annotation file(s).

### Expression covariates

The format for expression covariates is analogous to the format for gene expression described above. Covariates should be in context-specific, plain-text, tab-delimited files that begin with a header line. Each line after the header should contain information for a single individual with family ID in the first field, within-family ID in the second field, and covariates in the remaining fields. For example, the first three lines might be
```
FID	IID	PC1	PC2	InferredCov	pcr
0	indivA	0.0147	-0.0072	0.0262378174811602	1
0	indivB	0.0161	0.0037	-0.0514548756182194	1
```
Use the naming convention `<context>.covariates.txt` and save all of the covariate files in `drab/covariates`.

**Note:** The FIDs and IIDs must be consistent with those used in the gene expression files.

## Running DRAB

To run DRAB, submit the `src/run_drab.sh` shell script as a SLURM job with `CONTEXT_A`, `CONTEXT_B`, `GENES`, and `BOOT` as exported environment variables. For example, to test whether the expression of genes listed in the annotation file `all_genes.txt` is differentially regulated in tissues labeled as `Whole_Blood` and `Brain_Cortex` using 50 bootstrap iterations, run the command
```
sbatch --export=CONTEXT_A="Whole_Blood",CONTEXT_B="Brain_Cortex",GENES="all_genes",BOOT="50",DRAB=$(pwd) src/run_drab.sh
```
In practice, replace `Whole_Blood` and `Brain_Cortex` with the names of your desired tissues/contexts and `all_genes` with the name of your annotation file. The variable `BOOT` determines the number of bootstrap iterations to perform. Although we have found that as few as 10 bootstrap iterations will give reasonable results, we recommend using more iterations (e.g. 50) if computationally feasible.

The results will be saved to `output/Whole_Blood-Brain_Cortex-all_genes.txt`. (Here `Whole_Blood`, `Brain_Cortex`, and `all_genes` will be replaced with the tissue/context names and annotation file name you specified when running DRAB.) This is a tab-delimited, plain-text file without a header line. Each line contains information for a single gene, with the following fields:

1. Gene name
2. Gene ID
3. P-value for testing H_0: the gene is regulated identically in both contexts
4. P-value from a paired t-test for H_0: the trained models have equal squared prediction residuals
5. Number of individuals in each training set
6. Number of individuals in the testing set

If the DRAB P-value (in field 3) for a given gene is sufficiently small, then we conclude that the genetic regulation of that gene's expression is significantly different between the two contexts.

**Note:** DRAB does not account for multiple testing. When drawing conclusions on a set of genes, a multiple testing correction such as the Benjamini-Hochberg method should be used.

## Appendix: download and prepare GTEx data

This appendix describes how to obtain and prepare the data used in our paper "Identifying genes with tissue-specific patterns of genetic regulation."

First, create the folder `drab/raw` to store the unprocessed GTEx data sets. This folder may be safely deleted after completing all of the steps in this appendix.

From the **GTEx Analysis V8 (dbGaP Accession phs000424.v8.p2)** section of https://www.gtexportal.org/home/datasets, download the following files:

* `GTEx_Analysis_v8_eQTL_EUR.tar` (under the heading "Single-Tissue cis-QTL Data")  
Extract the folders `expression_matrices` and `expression_covariates` from the tar, and move them to `drab/raw`. (The other folder in the archive is not needed.)
* `gencode.v26.GRCh38.genes.gtf` (under the heading "Reference")  
Move this file to `drab/raw`.

After obtaining access to the GTEx data in [dbGaP](https://www.ncbi.nlm.nih.gov/gap/) (accession phs000424.v8.p2), follow the dbGaP documentation to download the following files:

* `phg001219.v1.GTEx_v8_WGS.genotype-calls-vcf.c1.GRU.tar`  
Extract the file `GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.vcf.gz` from the tar and move it to `drab/raw`. (The other files in the archive are not needed.)
* `phs000424.v8.pht002742.v8.p2.c1.GTEx_Subject_Phenotypes.GRU.txt.gz`  
Move this file to `drab/raw`.

To prepare the GTEx data for use with DRAB, from your main `drab` folder run
```
sbatch --export=DRAB=$(pwd) util/prepare_data.sh
```
This script will create an annotation file with all GTEx genes and another one with only protein-coding genes, perform standard quality control steps on the genotype data, and reformat the expression matrices and expression covariates.

