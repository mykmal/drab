# DREX: Differential Regulation of EXpression

This guide walks you through performing the drex workflow on European subjects data from the Genotype-Tissue Expression (GTEx) Project v8 release, as described in our paper. If using drex on your own data, make sure that it is formatted analogously to the GTEx data.

**Note:** The drex workflow is designed to be run on a linux HPC system and all provided commands are for bash (with the exception of a few R commands, which are prefaced by the `>` symbol).

## Setup

* First, clone the drex repository and create the required folder structure.
```
git clone https://github.com/MykMal/drex.git
cd drex
mkdir genotypes logs output plink reference temp weights
```
* Launch R and install the packages BEDMatrix, glmnet, insight, optparse, and RcppEigen. We used R v4.1.0 x86_64, BEDMatrix 2.0.3, glmnet 4.1.3, insight 0.15.0, optparse 1.7.1, and RcppEigen 0.3.3.9.1.
```
> install.packages(c('BEDMatrix','glmnet', 'insight', 'optparse', 'RcppEigen'))
```
* While still in the main `drex` folder, download and unpack the FUSION software.
```
wget https://github.com/gusevlab/fusion_twas/archive/refs/heads/master.zip
unzip master.zip
rm master.zip
mv fusion_twas-master fusion_twas
```
* Download the plink2R package to the `fusion_twas` folder.
```
cd fusion_twas
wget https://github.com/gabraham/plink2R/archive/refs/heads/master.zip
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
rm plink_linux_x86_64_20210606.zip
```
* Modify each of the files in `drex/scripts`, replacing `/path/to/drex` right below the HPC directives with the full path to your main `drex` folder. Also replace the HPC directives with appropriate ones for your HPC cluster. No other changes are needed, unless you want to adapt the workflow to your custom use case.

## Download data

Some data sets are available from the publicly-accessible [GTEx portal](https://www.gtexportal.org/home/), while others require General Research Use approval through [dbGaP](https://www.ncbi.nlm.nih.gov/gap/).

After obtaining access to the GTEx data in dbGaP (accession phs000424.v8.p2), follow the dbGaP documentation to download the following files:

* `phg001219.v1.GTEx_v8_WGS.genotype-calls-vcf.c1.GRU.tar`  
Extract the file `GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.vcf.gz` from the tar and move it to the `drex/genotypes` folder. (The other files in the tar are not needed.)
* `phs000424.v8.pht002742.v8.p2.c1.GTEx_Subject_Phenotypes.GRU.txt.gz`  
Uncompress this file and move it to the `drex/reference` folder.

From the **GTEx Analysis V8 (dbGaP Accession phs000424.v8.p2)** section of https://www.gtexportal.org/home/datasets, download the following files:

* `GTEx_Analysis_v8_eQTL_EUR.tar` (under the sub-heading "Single-Tissue cis-QTL Data")  
Unpack this archive and move the folders `expression_matrices` and `expression_covariates` to the main `drex` folder. (The other folder in the tar is not needed.)
* `GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt.gz` (under the sub-heading "Reference")  
Uncompress this file and move it to the `drex/reference` folder.
* `gencode.v26.GRCh38.genes.gtf` (under the sub-heading "Reference")  
Move this file to the `drex/reference` folder.

## Prepare data

Run the following shell scripts on your HPC system in any order:

* `process_covariates.sh` -- this will reformat the GTEx covariate data into the proper format for drex.
* `prepare_annotation.sh` -- this will create files in `drex/reference` with the name, chromosome, and location of each of the genes that drex will run on.
* `process_genotypes.sh` -- this will perform quality control on the GTEx genotype data, as described in our paper.

Optionally, run `split_data.sh` to randomly split the expression and covariate data for any tissue (specified in the first line of `code/SplitData.R`) into 1/2, 1/2 and 1/3, 2/3 partitions. You can then specify the partitions as separate "tissues" in drex to perform a null test.

## Run drex

As described in our paper, the drex method consists of three stages:

1. Select tissue-specific cis-eQTLs.
1. Train tissue-specific models that predict gene expression using the cis-eQTLs from the previous step as features.
1. Perform a non-nested model selection test on the trained prediction models.

### Stage 1

The first stage is performed by running `select_eqtls.sh`. This shell script runs on a per-tissue basis, so be sure to specify your desired tissue in the `TISSUE` environmental variable at the beginning of the script file prior to running. Note that this script will only run as a batch array job. If using SLURM, submit it with the command
```
sbatch --array=1-13 scripts/select_eqtls.sh
```
and after it finishes run
```
ls weights/<TISSUE> | tr ' ' '\n' | awk -F. '{print $2 "." $3}' > gene_list.txt && mv gene_list.txt weights/<TISSSUE>
```
where `<TISSUE>` should be replaced with the name of the tissue for which you are selecting eQTLs.  
The second number in the array flag above specifies the upper bound for the number of genes (in thousands) that you want to select eQTLs for. For example, to run drex on protein-coding genes in GTEx v8 you need to specify `--array=1-13` because there are 12,438 protein-coding genes in the GTEx gene model. To run drex on all genes in GTEx v8, you should instead specify `--array=1-36` because the GTEx data has annotations for 35,036 genes in total. For the latter case you will also need to change `cat reference/pc_gene_annotation.txt` to `cat reference/gene_annotation.txt` in line 30 of the script.

### Stage 2 and stage 3

The second and third stages are performed by running
```
Rscript code/DREX.R <TISSUE_A> <TISSUE_B>
```
where `<TISSUE_A>` and `<TISSUE_B>` should be replaced with the tissues between which you want to test for differential regulation of expression. Note that `select_eqtls.sh` needs to be run on both tissues prior to performing this step. The main drex script might take over 10 hours to run, depending on the sample sizes of the tissues and the number of genes being considered, so we recommend running it as an HPC job as well.  
When the process finishes, it will write the file `drex/output/<TISSUE_A>_<TISSUE_B>.txt`. This is a tab-separated file with one gene per row. The columns specify the gene Ensembl ID, the likelihood-ratio test p-value with `<TISSUE_A>` as the baseline, the likelihood-ratio test p-value with `<TISSUE_B>` as the baseline, the distribution-free test p-value with `<TISSUE_A>` as the baseline, and the distribution-free test p-value with `<TISSUE_B>` as the baseline. If the p-values for a given gene are significant, then we conclude that the genetic regulation of that gene's expression is significantly different between `<TISSUE_A>` and `<TISSUE_B>`.  
**Important:** The reported p-values are not corrected for multiple testing. A Bonferroni correction or some other appropriate family-wise error rate method should be applied before inference.
