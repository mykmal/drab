# DRAB: Differential Regulation Analysis by Bootstrapping

This repository provides code for the differential regulation analysis by bootstrapping (DRAB) method. The goal of DRAB is to identify genes with context-specific (e.g., tissue-specific) patterns of local genetic regulation. More generally, our method introduces a novel test for comparing machine learning models at the population level.

GitHub repo: <https://github.com/mykmal/drab>  
Paper: <https://doi.org/10.1214/23-AOAS1859>

DRAB is run for two tissues or other biological contexts at a time. First, elastic net regression models are trained within each context to learn the context-specific effects of genetic variation on gene expression. Next, a bootstrap-based model comparison test is applied to determine whether these context-specific models are equivalent. Our test accounts for the uncertainty of performing feature selection and model fitting on randomly sampled training sets, allowing us to interpret it as a test of equivalency between population-level models of genetic regulation. DRAB can be applied not only to gene expression, but also to any other functional or molecular phenotypes that have a genetic component, such as proteomic concentrations.

## Running DRAB

DRAB is primarily intended to be used on a Linux cluster through SLURM, but it can also be run as a shell script directly in a bash session.

### Installation

1. Download the DRAB repository and create the required folder structure.
```bash
wget https://github.com/mykmal/drab/archive/refs/heads/main.zip
unzip main.zip && mv drab-main drab && rm main.zip
cd drab
mkdir annotations covariates expression genotypes logs output
```
2. Install the R packages BEDMatrix and glmnet. We used R 4.2.2 x86_64, BEDMatrix 2.0.3, and glmnet 4.1-6.
```bash
Rscript -e 'install.packages(c("BEDMatrix", "glmnet"), repos="http://cran.us.r-project.org")'
```
3. Download PLINK to the main `drab` folder. We used PLINK v1.90b7 64-bit (16 Jan 2023).
```bash
wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20230116.zip
unzip plink_linux_x86_64_20230116.zip plink
rm plink_linux_x86_64_20230116.zip
```
4. The files `src/run_drab.sh` and `util/prepare_data.sh` are shell scripts prefaced by SLURM commands. If you intend to run DRAB through SLURM, modify the `#SBATCH` commands at the beginning of each script as appropriate for your high-performance computing system. Otherwise, if you will run DRAB without a job scheduling system, delete (or comment out) the `module load R/4.2.2-openblas` line near the top of the scripts and uncomment the two lines below it.

### Input data formats

The subsections below explain the files and file formats that DRAB expects. All of the files described here are required.

#### Gene annotations

DRAB requires a gene annotation file listing all of the genes on which it should run. This has to be a plain-text, tab-delimited file without a header line. Each line should contain information for a single gene, with the following five fields:

1. Gene name
2. Gene ID (e.g. from ENSEMBL)
3. Chromosome (with or without the chr prefix)
4. Start position (in base pairs)
5. End position (in base pairs)

For example, the lines for the first three guanylate binding protein genes might be:
```
GBP1	ENSG00000117228.9	chr1	89052319	89065360
GBP2	ENSG00000162645.12	chr1	89106132	89126114
GBP3	ENSG00000117226.11	chr1	89006666	89022894
```
Additional fields are allowed but will be ignored. Save your gene annotation files with file names of your choice in `drab/annotations`.

#### Genotype data

DRAB requires individual-level genotype data in PLINK 1.9 bed/bim/fam format. After performing all desired quality control steps, save your fully processed genotype data as `dosages.bed`, `dosages.bim`, and `dosages.fam` in `drab/genotypes`.

#### Gene expression data

DRAB requires individual-level gene expression data for each tissue or context of interest, and it is assumed that the RNA-Seq values have already been fully processed and normalized. Expression data should be in context-specific, plain-text, tab-delimited files that begin with a header line. Each line after the header should contain information for a single individual with family ID in the first field, within-family ID in the second field, and per-gene expression levels in the remaining fields. For example, the first three lines might be:
```
FID	IID	ENSG00000117228.9	ENSG00000162645.12	ENSG00000117226.11
0	indivA	-0.083051769477432	0.808844404113396	1.31169125330214
0	indivB	0.00672624465727554	-1.09866518781071	0.350055616620479
```
Use the naming convention `<context>.expression.txt` and save all of the gene expression files in `drab/expression`.

#### Expression covariates

The format for expression covariates is analogous to the format for gene expression described above. Covariates should be in context-specific, plain-text, tab-delimited files that begin with a header line. Each line after the header should contain information for a single individual with family ID in the first field, within-family ID in the second field, and covariates in the remaining fields. For example, the first three lines might be:
```
FID	IID	PC1	PC2	InferredCov	pcr
0	indivA	0.0147	-0.0072	0.0262378174811602	1
0	indivB	0.0161	0.0037	-0.0514548756182194	1
```
Use the naming convention `<context>.covariates.txt` and save all of the covariate files in `drab/covariates`.

### Example usage

The shell script `run_drab.sh` runs the program. It requires five arguments in the form of environment variables, which are described in the table below:

| Variable | Type | Description | Example |
|-|-|-|-|
| CONTEXT_A | string | prefix of gene expression filename for the first tissue/context | "Whole_Blood" |
| CONTEXT_B | string | prefix of gene expression filename for the second tissue/context | "Brain_Cortex" |
| GENES | string | basename of gene annotation file | "all_genes" |
| BOOT | integer | number of bootstrap iterations to use | 50 |
| DRAB | string | path to DRAB installation folder | "~/drab" |

Using the example values, the full command to run DRAB through SLURM from within the DRAB installation folder would be:
```bash
sbatch --export=CONTEXT_A="Whole_Blood",CONTEXT_B="Brain_Cortex",GENES="all_genes",BOOT="50",DRAB=$(pwd) src/run_drab.sh
```

To run DRAB without submitting a SLURM job, the commands would instead be:
```bash
export CONTEXT_A="Whole_Blood" CONTEXT_B="Brain_Cortex" GENES="all_genes" BOOT="50" DRAB=$(pwd)
./src/run_drab.sh
```

### Output format

The results will be saved to `output/<CONTEXT_A>-<CONTEXT_B>-<GENES>.txt`. (For the example above, this would be `output/Whole_Blood-Brain_Cortex-all_genes.txt`.) The output file is a tab-delimited, plain-text file without a header line. Each line contains information for a single gene, with the following six fields:

1. Gene name
2. Gene ID
3. $P$-value for testing $H_0$: the gene is regulated identically in both contexts
4. $P$-value from a paired $t$-test for $H_0$: the trained models have equal mean squared prediction errors
5. Number of individuals in each training set
6. Number of individuals in the test set

If the DRAB test $P$-value (in field 3) for a given gene is sufficiently small, then we can conclude that the genetic regulation of that gene's expression is significantly different between the two contexts. Note that the reported $P$-values are from single-gene tests, so a multiple testing correction may be necessary.

## Appendix A: Run DRAB with custom data splits

The main implementation of DRAB automatically splits the data you provide into a training set for context A ($D_A$), a training set for context B ($D_B$), and a test set ($D_T$). In some situations, however, it may be useful to manually define the sets $D_A$, $D_B$, and $D_T$. To facilitate such use cases, we provide the shell script `src/run_drab_manualsplit.sh`.

To use DRAB with pre-defined training and test sets, first split your gene expression data and expression covariates into separate files for $D_A$, $D_B$, and $D_T$. The expression and covariate files for each set should have the same prefix and end in `.expression.txt` and `.covariates.txt`, respectively. These files can be saved anywhere, but the expression data and covariate data for each set should stay within a single folder.

For example, suppose you saved your expression data in the files `da_split.expression.txt`, `db_split.expression.txt`, and `dt_split.expression.txt` and your covariate data in the files `da_split.covariates.txt`, `db_split.covariates.txt`, and `dt_split.covariates.txt`, all within a folder named `custom_splits`. Then the command to run DRAB through SLURM using those pre-defined data splits would be:
```bash
sbatch --export=DA_PATH="custom_splits/da_split",DB_PATH="custom_splits/db_split",DT_PATH="custom_splits/dt_split",GENES="all_genes",BOOT="50",DRAB=$(pwd) src/run_drab_manualsplit.sh
```

## Appendix B: Download and prepare GTEx data

This appendix describes how to obtain and prepare the data used in our paper.

First, create the folder `drab/raw` to store the unprocessed GTEx data sets. This folder may be safely deleted after completing all of the steps in this section.

Download the following files from the GTEx Portal website (<https://www.gtexportal.org/home/downloads/adult-gtex>):

* `GTEx_Analysis_v8_eQTL_EUR.tar` (under the ["Single-Tissue cis-QTL Data"](https://www.gtexportal.org/home/downloads/adult-gtex/qtl#qtl-gtex_analysis_v8-single-tissue_cis-qtl_data) heading in the "GTEx Analysis V8" section of the "QTL" tab). Extract the folders `expression_matrices` and `expression_covariates` from the tar, and move them to `drab/raw`. (The other folder in the archive is not needed.)
* `gencode.v26.GRCh38.genes.gtf` (under the ["Reference Tables"](https://www.gtexportal.org/home/downloads/adult-gtex/reference#reference-gtex_analysis_v8-reference_tables) heading in the "GTEx Analysis V8" section of the "Reference" tab). Move this file to `drab/raw`.

After obtaining access to the GTEx data in [dbGaP](https://www.ncbi.nlm.nih.gov/gap/) (accession number phs000424.v8.p2), follow the dbGaP documentation to download the following files:

* `phg001219.v1.GTEx_v8_WGS.genotype-calls-vcf.c1.GRU.tar`. Extract the file `GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.vcf.gz` from the tar and move it to `drab/raw`. (The other files in the archive are not needed.)
* `phs000424.v8.pht002742.v8.p2.c1.GTEx_Subject_Phenotypes.GRU.txt.gz`. Move this file to `drab/raw`.

To prepare the GTEx data for use with DRAB, from your main `drab` folder run:
```bash
sbatch --export=DRAB=$(pwd) util/prepare_data.sh
```
This script will create an annotation file with all GTEx genes and another one with only protein-coding genes, perform standard quality control steps on the genotype data, and reformat the expression matrices and expression covariates.

## Appendix C: Simulate gene expression data

This appendix describes how to simulate context-specific gene expression data using real genotype data, as done for the simulation studies described in our paper.

1. Create the required folder structure:
```bash
mkdir saved_models expression_simulated
```
2. Save an annotation file with the genes for which you wish to simulate gene expression levels, as described under the "Input data formats" section above.
3. Run the `simulations/simulate_expression.sh` shell script. This will train context-specific transcriptome imputation models, extract their weights, and then use those weights to simulate context-specific gene expression levels for all genotyped individuals. The required parameters for this script are the same as for `src/run_drab.sh`, except that the `BOOT` variable is not needed since no bootstrapping is performed. Below is an example run:
```bash
sbatch --export=CONTEXT_A="Whole_Blood",CONTEXT_B="Brain_Cortex",GENES="simulated_genes",DRAB=$(pwd) simulations/simulate_expression.sh
```

The simulated context-specific expression values for each gene will be saved to `expression_simulated/<CONTEXT>_<GENE ID>_expression.simulated.txt`. These files can then be used as inputs to DRAB.

