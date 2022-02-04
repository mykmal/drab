######################################################################
# INTERNAL FUNCTIONS
######################################################################

# Compute vectors of individual log-likelihoods for regression models with SNPs from snps_A, snps_B as features
GetLogLik <- function(features_1, features_2, expression)
{
  # Subset genotype and expression data to the shared list of individuals
  individuals <- intersect(rownames(expression), rownames(genotypes))
  expression <- expression[individuals, ]
  dosages <- genotypes[individuals, ]
  
  # Subset expression data for the current gene
  expression <- subset(expression, select = current_gene)
  colnames(expression) <- "X"
  
  # Get genotype data for the model SNPs
  dosages_1 <- data.frame(dosages[, features_1])
  dosages_2 <- data.frame(dosages[, features_2])
  
  # Remove SNPs for which standard deviation is zero
  dosages_1 <- Filter(sd, dosages_1)
  dosages_2 <- Filter(sd, dosages_2)
  
  # If no SNPs are left by now, we can't proceed further
  if ((dim(dosages_1)[2] == 0) | (dim(dosages_2)[2] == 0)) {
    cat("WARNING: no eQTLs found for", current_gene, ". Skipping gene.\n")
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
  model_1 <- lm(X ~ ., data = cbind(expression, dosages_1))
  
  # Fit a regression model with the SNPs in dosages_2 as features and expression as the response
  model_2 <- lm(X ~ ., data = cbind(expression, dosages_2))
  
  # Get numbers of observations
  n_1 <- insight::n_obs(model_1)
  n_2 <- insight::n_obs(model_2)
  
  # If the two models were trained on different numbers of individuals, something went wrong
  if (n_1 != n_2) {
    cat("WARNING: models for", current_gene, "have nonequal numbers of individuals. Skipping gene.\n")
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
  # AIC is used instead of the likelihood for model 1
  num <- sum(x$loglik_1) - x$coef_1 - sum(x$loglik_2)
  denom <-
    sqrt((1 / x$n) * sum((x$loglik_1 - x$loglik_2)**2) - ((1 / x$n) * sum(x$loglik_1 - x$loglik_2))**2)
  stat <- num / (sqrt(x$n) * denom)
  p <- 2 * pnorm(-abs(stat))
  
  return(p)
}

# Calculate the distribution-free test p-value
DFT <- function(x)
{
  # AIC-based correction is applied to model 1
  d <- x$loglik_1 - (x$coef_1 / x$n) - x$loglik_2
  stat <- sum(d > 0)
  b <- min(stat, x$n - stat)
  p <- 2 * pbinom(b, x$n, 0.5)
  
  return(p)
}

######################################################################
# MAIN PROGRAM
######################################################################

# Check that both arguments are supplied
args = commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Error: please specify two tissue names", call. = FALSE)
}

tissue_A <- args[1]
tissue_B <- args[2]

# Read expression data
expression_A <- read.delim(paste("expression_matrices/", tissue_A, ".v8.EUR.normalized_expression.bed.gz", sep = ""))
expression_B <- read.delim(paste("expression_matrices/", tissue_B, ".v8.EUR.normalized_expression.bed.gz", sep = ""))

# Read genotype data
genotypes <- BEDMatrix::BEDMatrix("genotypes/dosages_processed", simple_names = TRUE)

# Get names of all genes for which we have eQTLs in both tissues
genes_A <- read.table(file = paste("weights/", tissue_A, "/gene_list.txt", sep = ""))
genes_B <- read.table(file = paste("weights/", tissue_B, "/gene_list.txt", sep = ""))
genes <- intersect(genes_A, genes_B)

# Pre-populate a data frame to store p-values for all gene/test pairs
results <- data.frame(matrix(data = NA, nrow = length(genes), ncol = 5))
colnames(results) <- c(
    "gene_id",
    paste(tissue_A, "_lr_pval", sep = ""),
    paste(tissue_B, "_lr_pval", sep = ""),
    paste(tissue_A, "_df_pval", sep = ""),
    paste(tissue_B, "_df_pval", sep = ""))
results$gene_id <- genes

# Loop through all genes
for (current_gene in results$gene_id) {
  
  # Read tissue-specific eQTLs for the gene
  load(paste("weights/", tissue_A, "/", tissue_A, ".", current_gene, ".wgt.RDat", sep = ""))
  snps_A <- rownames(wgt.matrix[which(wgt.matrix[, "enet"] != 0), ])
  load(paste("weights/", tissue_B, "/", tissue_B, ".", current_gene, ".wgt.RDat", sep = ""))
  snps_B <- rownames(wgt.matrix[which(wgt.matrix[, "enet"] != 0), ])
  
  # Get individual log-likelihoods with tissue A as the baseline
  likelihoods_A <- GetLogLik(snps_A, snps_B, expression_A)
  
  # Get individual log-likelihoods with tissue B as the baseline
  likelihoods_B <- GetLogLik(snps_B, snps_A, expression_B)
  
  # Calculate p-values for each gene/test pair
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
  
  results[results$gene_id == current_gene, c(2, 3, 4, 5)] = cbind(lr_pval_A, lr_pval_B, df_pval_A, df_pval_B)
}

# Write all test results to a file
write.table(results,
  file = paste("output/", tissue_A, "_", tissue_B, ".txt", sep = ""),
  quote = FALSE,
  sep = "\t",
  row.names = FALSE)