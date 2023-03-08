# Overview of DRAB

DRAB (Differential Regulation Analysis by Bootstrapping) is a tool for identifying genes with context-specific (e.g. tissue-specific) patterns of local genetic regulation. DRAB first leverages elastic net regression to learn the effects of genetic variation on gene expression in each context, and then applies a bootstrap model comparison test to determine whether the context-specific models are equivalent. Notably, our approach is able to test population-level models by accounting for the variability of feature selection and model training.

DRAB can be applied to any functional/molecular phenotypes that have a genetic component, including mRNA expression and isoform expression levels.

## Citation

If you use this software, please star the repository and cite the following paper:

```bibtex
@article{malakhov_drab_2023,
     author = {Malakhov, Mykhaylo M. and Dai, Ben and Shen, Xiaotong T. and Pan, Wei},
     title = {A bootstrap model comparison test for identifying genes with context-specific patterns of genetic regulation},
     journal = {bioRxiv},
     publisher = {Cold Spring Harbor Laboratory},
     year = {2023},
     month = {03},
     doi = {10.1101/2023.03.06.531446},
     url = {https://www.biorxiv.org/content/10.1101/2023.03.06.531446}
}
```

# Running DRAB

DRAB is primarily intended to be used on a Linux cluster through SLURM, but it can also be run as a shell script from any bash session.

## Installation

1. Download the DRAB software package and create the required folder structure.
```
wget https://github.com/MykMal/drab/archive/refs/heads/main.zip
unzip main.zip && mv drab-main drab && rm main.zip
cd drab
mkdir annotations covariates expression genotypes logs output
```
2. Install the R packages BEDMatrix and glmnet. We used R 4.2.2 x86_64, BEDMatrix 2.0.3, and glmnet 4.1-6.
```
Rscript -e 'install.packages(c("BEDMatrix", "glmnet"), repos="http://cran.us.r-project.org")'
```
3. Download PLINK to the main `drab` folder. We used PLINK v1.90b7 64-bit (16 Jan 2023).
```
wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20230116.zip
unzip plink_linux_x86_64_20230116.zip plink
rm plink_linux_x86_64_20230116.zip
```
4. The files `src/run_drab.sh` and `util/prepare_data.sh` are shell scripts prefaced by SLURM commands. If you intend to run DRAB through SLURM, modify the `#SBATCH` commands at the beginning of each script as appropriate for your compute cluster. Otherwise, if you will run DRAB without a job scheduling system, delete (or comment out) the `module load R/4.2.2-openblas` line near the top of the scripts and uncomment the two lines below it.

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

### Expression covariates

The format for expression covariates is analogous to the format for gene expression described above. Covariates should be in context-specific, plain-text, tab-delimited files that begin with a header line. Each line after the header should contain information for a single individual with family ID in the first field, within-family ID in the second field, and covariates in the remaining fields. For example, the first three lines might be
```
FID	IID	PC1	PC2	InferredCov	pcr
0	indivA	0.0147	-0.0072	0.0262378174811602	1
0	indivB	0.0161	0.0037	-0.0514548756182194	1
```
Use the naming convention `<context>.covariates.txt` and save all of the covariate files in `drab/covariates`.

## Example usage

The shell script `run_drab.sh` runs the program. It requires five arguments in the form of environment variables, which are described in the table below:

| variable | type | description | example |
|-|-|-|-|
| CONTEXT_A | string | prefix of gene expression filename for the first tissue/context | "Whole_Blood" |
| CONTEXT_B | string | prefix of gene expression filename for the second tissue/context | "Brain_Cortex" |
| GENES | string | basename of gene annotation file | "all_genes" |
| BOOT | integer | number of bootstrap iterations to use | 50 |
| DRAB | string | path to DRAB installation directory | "~/drab" |

Using the example values, the full command to run DRAB through SLURM from within the DRAB installation directory would be
```
sbatch --export=CONTEXT_A="Whole_Blood",CONTEXT_B="Brain_Cortex",GENES="all_genes",BOOT="50",DRAB=$(pwd) src/run_drab.sh
```

To run DRAB without submitting a SLURM job, the commands would instead be
```
export CONTEXT_A="Whole_Blood" CONTEXT_B="Brain_Cortex" GENES="all_genes" BOOT="50" DRAB=$(pwd)
./src/run_drab.sh
```

## Output format

The results will be saved to `output/<CONTEXT_A>-<CONTEXT_B>-<GENES>.txt`. (For the example above, this would be `output/Whole_Blood-Brain_Cortex-all_genes.txt`.) The output file is a tab-delimited, plain-text file without a header line. Each line contains information for a single gene, with the following fields:

1. Gene name
2. Gene ID
3. P-value for testing H_0: the gene is regulated identically in both contexts
4. P-value from a paired t-test for H_0: the trained models have equal mean squared prediction errors
5. Number of individuals in each training set
6. Number of individuals in the testing set

If the DRAB test P-value (in field 3) for a given gene is sufficiently small, then we conclude that the genetic regulation of that gene's expression is significantly different between the two contexts. Note that the reported P-values are from single-gene tests, so a multiple testing correction may be necessary.

# Appendix: download and prepare GTEx data

This appendix describes how to obtain and prepare the data used in our paper.

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

