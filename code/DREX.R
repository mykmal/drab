######################################################################
# INTERNAL FUNCTIONS
######################################################################

# Read plink bed/bim/fam files
ReadPlink <- function(root)
{
  # Set full paths to each of the plink-format files
  bedfile <- paste(root, ".bed", sep = "")
  famfile <- paste(root, ".fam", sep = "")
  bimfile <- paste(root, ".bim", sep = "")
  
  # Import bed/bim/fam files
  # Note that specifying n and p speeds up BEDMatrix, since otherwise it would read the fam and bim files again
  bim <- read.table(bimfile, header = FALSE, stringsAsFactors = FALSE)
  fam <- read.table(famfile, header = FALSE, stringsAsFactors = FALSE)
  bed <- BEDMatrix::BEDMatrix(bedfile, n = nrow(fam), p = nrow(bim))
  
  # add individual and variant names to the imported bed
  rownames(bed) <- fam[, 2]
  colnames(bed) <- bim[, 2]
  
  return(list(
    bed = bed,
    bim = bim,
    fam = fam))
}

# Select eQTLs via elastic net regularization
ElasticNetSelection <- function(genos, pheno, alpha = 0.5)
{
  enet_model <- glmnet::cv.glmnet(x = genos, y = pheno, alpha = alpha, nfold = 5, intercept = TRUE, standardize = FALSE)
  eqtls <- rownames(coef(enet_model, s = "lambda.min"))[coef(enet_model, s = "lambda.min")[, 1] != 0]
  
  return(eqtls[-1])
}

# Compute vectors of individual log-likelihoods for regression models with expression as the response and SNPs from features_1, features_2 as predictors
GetLogLik <- function(features_1, features_2, genotypes, expression)
{
  # Get dosages for the SNPs in features_1, features_2
  features_1 <- intersect(features_1, colnames(genotypes))
  features_2 <- intersect(features_2, colnames(genotypes))
  dosages_1 <- genotypes[, features_1]
  dosages_2 <- genotypes[, features_2]
  
  # If no eQTLs exist, we can't proceed further
  if ((length(features_1) == 0) | (length(features_2) == 0)) {
    cat("WARNING: no eQTLs selected. Skipping gene.\n")
    return(list(
      n = 0,
      coef_1 = 0,
      coef_2 = 0,
      loglik_1 = 0,
      loglik_2 = 0))
  }
  
  # If multiple SNPs are present, remove the highly correlated ones
  if (dim(dosages_1)[2] > 1) {
    tmp <- cor(dosages_1)
    tmp[!lower.tri(tmp)] <- 0
    dosages_1 <- dosages_1[, apply(tmp, 2, function(x) all(abs(x) <= 0.9, na.rm = TRUE))]
  }
  if (dim(dosages_2)[2] > 1) {
    tmp <- cor(dosages_2)
    tmp[!lower.tri(tmp)] <- 0
    dosages_2 <- dosages_2[, apply(tmp, 2, function(x) all(abs(x) <= 0.9, na.rm = TRUE))]
  }
  
  # Fit a regression model with the SNPs in dosages_1 as features and expression as the response
  model_1 <- lm(expression ~ ., data = as.data.frame(dosages_1))
  
  # Fit a regression model with the SNPs in dosages_2 as features and expression as the response
  model_2 <- lm(expression ~ ., data = as.data.frame(dosages_2))
  
  # Get numbers of observations
  n_1 <- insight::n_obs(model_1)
  n_2 <- insight::n_obs(model_2)
  
  # If the two models were trained on different numbers of individuals, something went wrong
  if (n_1 != n_2) {
    cat("WARNING: prediction models have nonequal numbers of individuals. Skipping gene.\n")
    return(list(
      n = 0,
      coef_1 = 0,
      coef_2 = 0,
      loglik_1 = 0,
      loglik_2 = 0))
  }
  
  # Get numbers of parameters
  coef_1 <- insight::n_parameters(model_1)
  coef_2 <- insight::n_parameters(model_2)
  
  # Get individual log-likelihoods
  loglik_1 <- attributes(insight::get_loglikelihood(model_1))$per_obs
  loglik_2 <- attributes(insight::get_loglikelihood(model_2))$per_obs
  
  return(list(
    n = n_1,
    coef_1 = coef_1,
    coef_2 = coef_2,
    loglik_1 = loglik_1,
    loglik_2 = loglik_2))
}

# Calculate the likelihood-ratio test p-value
LRT <- function(x)
{
  # AIC-based correction is applied to the model 1 likelihood
  numerator <- sum(x$loglik_1) - x$coef_1 - sum(x$loglik_2)
  variance <- (1 / x$n) * sum((x$loglik_1 - x$loglik_2)**2) - ((1 / x$n) * sum(x$loglik_1 - x$loglik_2))**2
  statistic <- numerator / sqrt(x$n * variance)
  p <- 2 * pnorm(-abs(stat))
  
  return(p)
}

# Calculate the distribution-free test p-value
DFT <- function(x)
{
  # AIC-based correction is applied to the model 1 likelihood
  d <- x$loglik_1 - (x$coef_1 / x$n) - x$loglik_2
  b <- sum(d > 0)
  statistic <- min(b, x$n - b)
  p <- 2 * pbinom(statistic, x$n, 0.5)
  
  return(p)
}

######################################################################
# MAIN PROGRAM
######################################################################

# Check that all required arguments are supplied
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7) {
  cat(
    "ERROR: missing arguments. Specify (in the following order) a gene name, the gene's Ensembl ID, the first tissue name, the path to plink binary files (minus .bed/.bim/.fam) for the first tissue, the second tissue name, the path to plink binary files (minus .bed/.bim/.fam) for the second tissue, and a batch number.\n",
    file = stderr())
  q()
}

gene <- args[1]
gene_id <- args[2]
tissue_A <- args[3]
plink_path_A <- args[4]
tissue_B <- args[5]
plink_path_B <- args[6]
batch <- args[7]

covar_path_A <- paste("expression_covariates/", tissue_A, ".v8.EUR.covariates.plink.txt", sep = "")
covar_path_B <- paste("expression_covariates/", tissue_B, ".v8.EUR.covariates.plink.txt", sep = "")

# Check that required files exist
files_A <- paste(plink_path_A, c(".bed", ".bim", ".fam"), sep = "")
files_B <- paste(plink_path_B, c(".bed", ".bim", ".fam"), sep = "")
files <- c(files_A, files_B, covar_path_A, covar_path_B)

for (f in files) {
  if (!file.exists(f)) {
    cat("ERROR:", f , "input file does not exist\n", file = stderr())
    q()
  }
}

# Load genotype and expression data
data_A <- ReadPlink(plink_path_A)
data_B <- ReadPlink(plink_path_B)
genotypes_A <- data_A$bed
genotypes_B <- data_B$bed
expression_A <- data_A$fam[, c(2, 6)]
expression_B <- data_B$fam[, c(2, 6)]

# Load covariates
covar_A <- read.table(covar_path_A, header = TRUE, stringsAsFactors = FALSE)[, -1]
covar_B <- read.table(covar_path_B, header = TRUE, stringsAsFactors = FALSE)[, -1]
cat("Loaded", ncol(covar_A)-1, "covariates for", tissue_A, "and", ncol(covar_B)-1, "covariates for", tissue_B, "\n")

# Find individuals for which we have both covariate and genotype/expression data
individuals_A <- intersect(expression_A[, 1], covar_A[, 1])
individuals_B <- intersect(expression_B[, 1], covar_B[, 1])

# Subset genotype, expression, and covariate data to the shared lists of individuals
genotypes_A <- genotypes_A[individuals_A, ]
genotypes_B <- genotypes_B[individuals_B, ]
expression_A <- expression_A[expression_A[, 1] %in% individuals_A, ]
expression_B <- expression_B[expression_B[, 1] %in% individuals_B, ]
covar_A <- covar_A[covar_A[, 1] %in% individuals_A, ]
covar_B <- covar_B[covar_B[, 1] %in% individuals_B, ]

# Adjust expression values for covariates
expression_covar_A <- summary(lm(expression_A[, 2] ~ ., data = as.data.frame(covar_A[, -1])))
cat(expression_covar_A$r.squared, "variance in", tissue_A, "expression explained by covariates\n")
expression_covar_B <- summary(lm(expression_B[, 2] ~ ., data = as.data.frame(covar_B[, -1])))
cat(expression_covar_B$r.squared, "variance in", tissue_B, "expression explained by covariates\n")

# Scale and center the genotypes and adjusted expression values
genotypes_A <- scale(genotypes_A)
genotypes_B <- scale(genotypes_B)
expression_A[, 2] <- scale(expression_covar_A$residuals)
expression_B[, 2] <- scale(expression_covar_B$residuals)

# Check if any genotypes are NA
na_snps_A <- apply(!is.finite(genotypes_A), 2, sum)
if (sum(na_snps_A) != 0) {
  cat("WARNING:", sum(na_snps_A != 0), "SNPs could not be scaled and were zeroed out in", tissue_A, "\n", file = stderr())
  genotypes_A[, na_snps_A != 0] <- 0
}
na_snps_B <- apply(!is.finite(genotypes_B), 2, sum)
if (sum(na_snps_B) != 0) {
  cat("WARNING:", sum(na_snps_B != 0), "SNPs could not be scaled and were zeroed out in", tissue_B, "\n", file = stderr())
  genotypes_B[, na_snps_B != 0] <- 0
}

# Adjust genotypes for covariates
for (i in seq_len(ncol(genotypes_A))) {
  genotypes_A[, i] <- summary(lm(genotypes_A[, i] ~ ., data = as.data.frame(covar_A[, -1])))$residuals
}
for (i in seq_len(ncol(genotypes_B))) {
  genotypes_B[, i] <- summary(lm(genotypes_B[, i] ~ ., data = as.data.frame(covar_B[, -1])))$residuals
}

# Scale and center the genotypes again
genotypes_A <- scale(genotypes_A)
genotypes_B <- scale(genotypes_B)

# Remove monomorphic SNPs
sds_A <- apply(genotypes_A, 2, sd)
keep_A <- (sds_A != 0) & (!is.na(sds_A))
genotypes_A <- genotypes_A[, keep_A]
sds_B <- apply(genotypes_B, 2, sd)
keep_B <- (sds_B != 0) & (!is.na(sds_B))
genotypes_B <- genotypes_B[, keep_B]

cat(tissue_A, ": ", nrow(expression_A), " individuals and ", ncol(genotypes_A), " cis-SNPs\n", sep = "")
cat(tissue_B, ": ", nrow(expression_B), " individuals and ", ncol(genotypes_B), " cis-SNPs\n", sep = "")

# Select eQTLs for each tissue
eqtls_A <- ElasticNetSelection(genotypes_A, expression_A[, 2])
eqtls_B <- ElasticNetSelection(genotypes_B, expression_B[, 2])

cat(tissue_A, ": ", length(eqtls_A), " cis-eQTLs selected\n", sep = "")
cat(tissue_B, ": ", length(eqtls_B), " cis-eQTLs selected\n", sep = "")

# Get individual log-likelihoods with tissue A as the baseline
likelihoods_A <- GetLogLik(eqtls_A, eqtls_B, genotypes_A, expression_A[, 2])

# Get individual log-likelihoods with tissue B as the baseline
likelihoods_B <- GetLogLik(eqtls_B, eqtls_A, genotypes_B, expression_B[, 2])

# Calculate p-values for each pair of baseline tissue and test
if (likelihoods_A$n == 0) {
  lr_pval_A <- "NA"
  df_pval_A <- "NA"
} else {
  lr_pval_A <- LRT(likelihoods_A)
  df_pval_A <- DFT(likelihoods_A)
}
if (likelihoods_B$n == 0) {
  lr_pval_B <- "NA"
  df_pval_B <- "NA"
} else {
  lr_pval_B <- LRT(likelihoods_B)
  df_pval_B <- DFT(likelihoods_B)
}

# Append the test results to the working file
results <- paste(gene, gene_id, lr_pval_A, lr_pval_B, df_pval_A, df_pval_B, sep = "\t")
cat(results, file = paste("temp/", tissue_A, "_", tissue_B, "_", batch, sep = ""), append = TRUE, sep = "\n")

