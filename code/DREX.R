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
  bed <- as.data.frame(as.matrix(BEDMatrix::BEDMatrix(bedfile, n = nrow(fam), p = nrow(bim))))
  
  # Add individual and variant names to the imported bed
  rownames(bed) <- fam$V2
  colnames(bed) <- bim$V2
  
  # Remove unused fam columns and rename the remaining ones
  fam <- subset(fam, select = c(V2, V6))
  colnames(fam) <- c("IID", "expression")
  
  return(list(
    genotypes = bed,
    expression = fam))
}

# Load expression and genotype data, adjust for covariates, and normalize
CovarAdjust <- function(tissue, plink_path)
{
  covar_path_train <- paste("temp/", tissue, ".", batch, "/", tissue, "_half1.v8.EUR.covariates.plink.txt", sep = "")
  covar_path_test <- paste("temp/", tissue, ".", batch, "/", tissue, "_half2.v8.EUR.covariates.plink.txt", sep = "")
  
  # Check that required files exist
  plink_files <- paste(plink_path, c("_half1.bed", "_half2.bed", "_half1.bim", "_half2.bim", "_half1.fam", "_half2.fam"), sep = "")
  files <- c(plink_files, covar_path_train, covar_path_test)
  
  for (f in files) {
    if (!file.exists(f)) {
      cat("ERROR:", f , "input file does not exist\n", file = stderr())
      q()
    }
  }
  
  # Load genotype and expression data
  data_train <- ReadPlink(paste(plink_path, "_half1", sep = ""))
  data_test <- ReadPlink(paste(plink_path, "_half2", sep = ""))
  genotypes_train <- data_train$genotypes
  genotypes_test <- data_test$genotypes
  expression_train <- data_train$expression
  expression_test <- data_test$expression
  
  # Load covariates
  covar_train <- read.table(covar_path_train, header = TRUE, stringsAsFactors = FALSE)
  covar_train <- subset(covar_train, select = -FID)
  covar_test <- read.table(covar_path_test, header = TRUE, stringsAsFactors = FALSE)
  covar_test <- subset(covar_test, select = -FID)
  cat(tissue, ": Loaded ", ncol(covar_train)-1, " covariates for training data and ", ncol(covar_test)-1, " covariates for testing data\n", sep = "")
  
  # Match expression data with covariates by individual IDs
  merged_expression_train <- merge(expression_train, covar_train)
  rownames(merged_expression_train) <- merged_expression_train$IID
  merged_expression_train <- subset(merged_expression_train, select = -IID)
  merged_expression_test <- merge(expression_test, covar_test)
  rownames(merged_expression_test) <- merged_expression_test$IID
  merged_expression_test <- subset(merged_expression_test, select = -IID)
  
  # Adjust expression values for covariates
  expression_train_summary <- summary(lm(expression ~ ., data = merged_expression_train))
  cat(tissue, "training data:", expression_train_summary$r.squared, "variance in expression explained by covariates\n")
  expression_test_summary <- summary(lm(expression ~ ., data = merged_expression_test))
  cat(tissue, "testing data:", expression_test_summary$r.squared, "variance in expression explained by covariates\n")
  
  # Scale and center the genotypes and the adjusted expression values
  genotypes_train <- as.data.frame(scale(genotypes_train))
  genotypes_test <- as.data.frame(scale(genotypes_test))
  expression_train <- as.data.frame(scale(expression_train_summary$residuals))
  colnames(expression_train) <- "expression"
  expression_test <- as.data.frame(scale(expression_test_summary$residuals))
  colnames(expression_test) <- "expression"
  
  # Check if any genotypes are NA
  na_snps_train <- apply(genotypes_train, 2, function(x) !all(is.finite(x)))
  if (sum(na_snps_train) != 0) {
    cat("WARNING:", sum(na_snps_train), "SNPs could not be scaled and were removed in", tissue, "training data\n")
    genotypes_train <- subset(genotypes_train, select = !na_snps_train)
  }
  na_snps_test <- apply(genotypes_test, 2, function(x) !all(is.finite(x)))
  if (sum(na_snps_test) != 0) {
    cat("WARNING:", sum(na_snps_test), "SNPs could not be scaled and were removed in", tissue, "testing data\n")
    genotypes_test <- subset(genotypes_test, select = !na_snps_test)
  }
  
  # Match genotypes with covariates by individual IDs
  merged_genotypes_train <- cbind(rownames(genotypes_train), genotypes_train)
  colnames(merged_genotypes_train) <- c("IID", colnames(genotypes_train))
  merged_genotypes_train <- merge(merged_genotypes_train, covar_train)
  rownames(merged_genotypes_train) <- merged_genotypes_train$IID
  merged_genotypes_train <- subset(merged_genotypes_train, select = -IID)
  merged_genotypes_test <- cbind(rownames(genotypes_test), genotypes_test)
  colnames(merged_genotypes_test) <- c("IID", colnames(genotypes_test))
  merged_genotypes_test <- merge(merged_genotypes_test, covar_test)
  rownames(merged_genotypes_test) <- merged_genotypes_test$IID
  merged_genotypes_test <- subset(merged_genotypes_test, select = -IID)
  
  # Adjust genotypes for covariates
  for (snp in colnames(genotypes_train)) {
    genotypes_train_summary <- summary(lm(reformulate(colnames(covar_train)[-1], response = snp), data = merged_genotypes_train))
    merged_genotypes_train[[snp]] <- genotypes_train_summary$residuals
  }
  for (snp in colnames(genotypes_test)) {
    genotypes_test_summary <- summary(lm(reformulate(colnames(covar_test)[-1], response = snp), data = merged_genotypes_test))
    merged_genotypes_test[[snp]] <- genotypes_test_summary$residuals
  }
  
  genotypes_train <- merged_genotypes_train[, 1:ncol(genotypes_train)]
  genotypes_test <- merged_genotypes_test[, 1:ncol(genotypes_test)]
  
  # Scale and center the genotypes again
  genotypes_train <- as.data.frame(scale(genotypes_train))
  genotypes_test <- as.data.frame(scale(genotypes_test))
  
  # Remove monomorphic SNPs
  sds_train <- apply(genotypes_train, 2, sd)
  keep_train <- (sds_train != 0) & (is.finite(sds_train))
  genotypes_train <- subset(genotypes_train, select = keep_train)
  sds_test <- apply(genotypes_test, 2, sd)
  keep_test <- (sds_test != 0) & (is.finite(sds_test))
  genotypes_test <- subset(genotypes_test, select = keep_test)
  
  cat(tissue, "training data:", nrow(expression_train), "individuals and", ncol(genotypes_train), "cis-SNPs\n")
  cat(tissue, "testing data:", nrow(expression_test), "individuals and", ncol(genotypes_test), "cis-SNPs\n")
  
  return(list(
    genotypes_train = genotypes_train,
    genotypes_test = genotypes_test,
    expression_train = expression_train,
    expression_test = expression_test))
}

# Select eQTLs via elastic net regularization
ElasticNetSelection <- function(genos, pheno, alpha = 0.5)
{
  enet_model <- glmnet::cv.glmnet(x = as.matrix(genos), y = as.matrix(pheno), alpha = alpha, nfold = 5, intercept = TRUE, standardize = FALSE)
  eqtls <- rownames(coef(enet_model, s = "lambda.1se"))[coef(enet_model, s = "lambda.1se")[, 1] != 0]
  
  return(eqtls[-1])
}

# Compute vectors of individual log-likelihoods for regression models with expression as the response and SNPs from eqtls_A, eqtls_B as predictors
GetLogLik <- function(tissue, genotypes, expression)
{
  # Get dosages for the SNPs in eqtls_A, eqtls_B
  features_A <- intersect(eqtls_A, colnames(genotypes))
  features_B <- intersect(eqtls_B, colnames(genotypes))
  dosages_A <- subset(genotypes, select = features_A)
  dosages_B <- subset(genotypes, select = features_B)
  
  # If no eQTLs exist, we can't proceed further
  if ((ncol(dosages_A) == 0) | (ncol(dosages_B) == 0)) {
    cat("WARNING: no eQTLs selected in ", tissue, ". Skipping gene.\n", sep = "")
    return(list(
      n = 0,
      coef_A = 0,
      coef_B = 0,
      loglik_A = 0,
      loglik_B = 0))
  }
  
  # If multiple SNPs are present, remove the highly correlated ones
  if (ncol(dosages_A) > 1) {
    tmp <- cor(dosages_A)
    tmp[!lower.tri(tmp)] <- 0
    dosages_A <- dosages_A[, apply(tmp, 2, function(x) all(abs(x) < 0.8, na.rm = TRUE))]
  }
  if (ncol(dosages_B) > 1) {
    tmp <- cor(dosages_B)
    tmp[!lower.tri(tmp)] <- 0
    dosages_B <- dosages_B[, apply(tmp, 2, function(x) all(abs(x) < 0.8, na.rm = TRUE))]
  }
  
  # Scale dosages
  dosages_A <- as.data.frame(scale(dosages_A))
  dosages_B <- as.data.frame(scale(dosages_B))
  
  # Match dosage and expression data by individual IDs
  merged_data_A <- merge(expression, dosages_A, by = 0)
  rownames(merged_data_A) <- merged_data_A$Row.names
  merged_data_A <- subset(merged_data_A, select = -Row.names)
  merged_data_B <- merge(expression, dosages_B, by = 0)
  rownames(merged_data_B) <- merged_data_B$Row.names
  merged_data_B <- subset(merged_data_B, select = -Row.names)
  
  # Fit a regression model with the SNPs in dosages_A as features and expression as the response
  model_A <- lm(expression ~ ., data = merged_data_A)
  
  # Fit a regression model with the SNPs in dosages_B as features and expression as the response
  model_B <- lm(expression ~ ., data = merged_data_B)
  
  cat(tissue, ": ", summary(model_A)$r.squared, " variance explained by ", tissue_A, "-specific eQTLS and ", summary(model_B)$r.squared, " variance explained by ", tissue_B, "-specific eQTLs\n", sep = "")
  
  # Get numbers of parameters
  coef_A <- insight::n_parameters(model_A)
  coef_B <- insight::n_parameters(model_B)
  
  # Get individual log-likelihoods
  loglik_A <- attributes(insight::get_loglikelihood(model_A))$per_obs
  loglik_B <- attributes(insight::get_loglikelihood(model_B))$per_obs
  
  return(list(
    n = insight::n_obs(model_A),
    coef_A = coef_A,
    coef_B = coef_B,
    loglik_A = loglik_A,
    loglik_B = loglik_B))
}

# Calculate the likelihood-ratio test p-value
LRT <- function(x, correction = "both")
{
  # No AIC-based correction is applied
  if (correction == "none")
    numerator <- sum(x$loglik_A) - sum(x$loglik_B)
  # AIC-based correction is applied to the model A log likelihood
  if (correction == "A")
    numerator <- sum(x$loglik_A) - x$coef_A - sum(x$loglik_B)
  # AIC-based correction is applied to the model B log likelihood
  if (correction == "B")
    numerator <- sum(x$loglik_A) - sum(x$loglik_B) + x$coef_B
  # AIC-based correction is applied to both log likelihoods
  else
    numerator <- sum(x$loglik_A) - x$coef_A - sum(x$loglik_B) + x$coef_B
  
  variance <- (1 / x$n) * sum((x$loglik_A - x$loglik_B)**2) - ((1 / x$n) * sum(x$loglik_A - x$loglik_B))**2
  statistic <- numerator / sqrt(x$n * variance)
  p <- 2 * pnorm(-abs(statistic))
  
  return(p)
}

# Calculate the distribution-free test p-value
DFT <- function(x, correction = "both")
{
  # No AIC-based correction is applied
  if (correction == "none")
    d <- x$loglik_A - x$loglik_B
  # AIC-based correction is applied to the model A log likelihood
  if (correction == "A")
    d <- x$loglik_A - (x$coef_A / x$n) - x$loglik_B
  # AIC-based correction is applied to the model B log likelihood
  if (correction == "B")
    d <- x$loglik_A - x$loglik_B + (x$coef_B / x$n)
  # AIC-based correction is applied to both log likelihoods
  else
    d <- x$loglik_A - (x$coef_A / x$n) - x$loglik_B + (x$coef_B / x$n)
  
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
  cat("ERROR: missing arguments.\n", file = stderr())
  q()
}

gene <- args[1]
gene_id <- args[2]
tissue_A <- args[3]
plink_path_A <- args[4]
tissue_B <- args[5]
plink_path_B <- args[6]
batch <- args[7]

# Import expression and genotype data, adjust it for covariates, and normalize it
data_A <- CovarAdjust(tissue_A, plink_path_A)
data_B <- CovarAdjust(tissue_B, plink_path_B)

# Select eQTLs for each tissue
eqtls_A <- ElasticNetSelection(data_A$genotypes_train, data_A$expression_train)
eqtls_B <- ElasticNetSelection(data_B$genotypes_train, data_B$expression_train)

cat(tissue_A, ": ", length(eqtls_A), " cis-eQTLs selected\n", sep = "")
cat(tissue_B, ": ", length(eqtls_B), " cis-eQTLs selected\n", sep = "")

# Train expression prediction models with tissue-specific eQTLs on data from tissue A, and return their log-likelihoods
likelihoods_A <- GetLogLik(tissue_A, data_A$genotypes_test, data_A$expression_test)

# Train expression prediction models with tissue-specific eQTLs on data from tissue B, and return their log-likelihoods
likelihoods_B <- GetLogLik(tissue_B, data_B$genotypes_test, data_B$expression_test)

cat("\n")

# Calculate p-values for each pair of baseline tissue and model selection test
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
cat(results, file = paste("temp/out_", tissue_A, "_", tissue_B, "_", batch, ".txt", sep = ""), append = TRUE, sep = "\n")

