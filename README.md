# DREX: Differential Regulation of EXpression

This tutorial walks you through performing the drex workflow on European subjects data from the Genotype-Tissue Expression Project (GTEx) v8 release, as described in our paper. If using drex on your own data, make sure that it is formatted analogously to the GTEx data.

**Note:** the drex workflow is designed to be run on a linux HPC system and all provided commands are for bash (with the exception of a few R commands, which are prefaced by the `>` symbol).

## Setup

* First, clone the drex repository and create the required folder structure.
```
git clone https://github.com/MykMal/drex.git
cd drex
mkdir genotypes output plink reference temp weights
```
* Launch R and install the packages BEDMatrix, glmnet, optparse, and RcppEigen. We used R v4.1.0 x86_64, BEDMatrix 2.0.3, glmnet 4.1-3, optparse 1.7.1, and RcppEigen 0.3.3.9.1.
```
> install.packages(c('BEDMatrix','glmnet', 'optparse', 'RcppEigen'))
```
* While still in the main `drex` folder, download and unpack the FUSION software.
```
wget https://github.com/gusevlab/fusion_twas/archive/master.zip
unzip master.zip
rm master.zip
mv fusion_twas-master fusion_twas
```
* Download the plink2R package to the `fusion_twas` folder.
```
cd fusion_twas
wget https://github.com/gabraham/plink2R/archive/master.zip
unzip master.zip
rm master.zip
```
Launch R and install the plink2R package.
```
> install.packages('plink2R-master/plink2R/', repos = NULL)
```
* Download plink to the `plink` folder. We tested this workflow with plink v1.90b6.24 64-bit (6 Jun 2021), but other versions of plink 1.90 should work too.
```
cd ../plink
wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20210606.zip
unzip plink_linux_x86_64_20210606.zip
```
* Modify each of the files in `drex/scripts`, replacing `/path/to/drex` with the full path to your main `drex` folder. Also replace the HPC directives with appropriate ones for your cluster. No other changes are needed, unless you want to adapt the workflow to your custom use case.

## Download data

Some data sets are available from the publicly-accessible GTEx portal at [https://www.gtexportal.org/home/datasets](https://www.gtexportal.org/home/datasets), while others require General Research Use approval through [dbGaP](https://www.ncbi.nlm.nih.gov/gap/).

After obtaining access to the GTEx data in dbGaP (accession phs000424.v8.p2), follow the dbGaP documentation to download the following files:

* `phg001219.v1.GTEx_v8_WGS.genotype-calls-vcf.c1.GRU.tar`. Extract the file `GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.vcf.gz` from the tar and move it to the `drex/genotypes` folder. (The other files in the tar are not needed.)
* `phs000424.v8.pht002742.v8.p2.c1.GTEx_Subject_Phenotypes.GRU.txt.gz`. Uncompress this file and move it to the `drex/reference` folder.

From the **GTEx Analysis V8 (dbGaP Accession phs000424.v8.p2)** section of the GTEx portal, download the following files:

* `GTEx_Analysis_v8_eQTL_EUR.tar` (under the sub-heading "Single-Tissue cis-QTL Data"). Unpack this archive and move the folders `expression_matrices` and `expression_covariates` to the main `drex` folder. (The other folders in the tar are not needed.)
* `GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt.gz` (under the sub-heading "Reference"). Uncompress this file and move it to the `drex/reference` folder.
* `gencode.v26.GRCh38.genes.gtf` (under the sub-heading "Reference"). Move this file to the `drex/reference` folder.

## Prepare data

Run the following shell scripts on your HPC system in any order:

* `process_covariates.sh`. This will reformat the GTEx covariate data into the proper format for drex.
* `prepare_annotation.sh`. This will create files in `drex/reference` with the name, chromosome, and location of each of the genes that drex will run on.
* `process_genotypes.sh`. This will perform quality control on the GTEx genotype data, as described in our paper.

Optionally, run the `split_data.sh` if you want to verify that drex properly controls false positives. This will randomly split the expression and covariate data for any tissue (specified in the first line) into 1/2, 1/2 and 1/3, 2/3 partitions. You can then specify the partitions as separate "tissues" in drex to perform a null test.

## Run drex

As described in our paper, drex consists of three stages:

1. Select tissue-specific cis-eQTLs.
1. Train tissue-specific models that predict gene expression using the cis-eQTLs from the previous step as features.
1. Perform a non-nested model selection test on the trained prediction models.

The first stage is performed by running `select_eqtls.sh`. This shell script runs on a per-tissue basis, so be sure to modify it to specify the desired tissue before running. Note that as written, this script will only run as a batch SLURM job. If using SLURM, submit it with the command
```
sbatch --array=1-13 train_models.sh
```
The second number in the array specifies the upper bound for the number of genes (in thousands) that you want to train models for. For example, to run drex on protein-coding genes in GTEx v8 specify `--array=1-13` because there are 12,438 protein-coding genes. To run drex on all genes in GTEx v8, then instead specify `--array=1-36` because there are 35,036 genes in total. For the latter case you will also need to change `cat reference/pc_gene_annotation.txt` to `cat reference/gene_annotation.txt` in line 28 of the script.  

The second and third stages are performed by running
```
Rscript code/DREX.R <TISSUE_A> <TISSUE_B>
```
where `<TISSUE_A>` and `<TISSUE_B>` should be replaced with the tissues between which you want to test for differential regulation. Note that `select_eqtls.sh` needs to be run on both tissues prior to performing this step. This might take several hours depending on the sample sizes of the tissues and the number of genes being considered, so we recommend running this step as a batch HPC job as well. When the process finishes, it will write the file `drex/output/<TISSUE_A>_<TISSUE_B>.txt`. This is a tab-separated file with one gene per row. The columns specify the gene Ensembl ID, the likelihood-ratio test p-value with `<TISSUE_A>` as the baseline, the likelihood-ratio test p-value with `<TISSUE_B>` as the baseline, the distribution-free test p-value with `<TISSUE_A>` as the baseline, and the distribution-free test p-value with `<TISSUE_B>` as the baseline.

