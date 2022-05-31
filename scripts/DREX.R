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
  
  # Import bed/bim/fam files.
  # Note that specifying n and p speeds up BEDMatrix, since otherwise it would read the fam and bim files again.
  bim <- read.table(bimfile, header = FALSE, stringsAsFactors = FALSE)
  fam <- read.table(famfile, header = FALSE, stringsAsFactors = FALSE)
  bed <- as.data.frame(as.matrix(BEDMatrix::BEDMatrix(bedfile, n = nrow(fam), p = nrow(bim))))
  
  # Add individual and variant names to the imported bed.
  # Individual names are coded as FID_IID;
  # Variant names are coded as rsid_A1.
  rownames(bed) <- paste(fam$V1, fam$V2, sep = "_")
  colnames(bed) <- paste(bim$V2, bim$V5, sep = "_")
  
  # Set row/column names and remove unused columns in the imported fam
  rownames(fam) <- paste(fam$V1, fam$V2, sep = "_")
  fam <- subset(fam, select = V6)
  colnames(fam) <- "value"
  
  return(list(genotypes = bed, expression = fam))
}

# Remove columns with NAs or a standard deviation of zero
RemoveSuperfluous <- function(x)
{
  sds <- apply(x, 2, sd, na.rm = FALSE)
  keep <- (sds != 0) & (is.finite(sds))
  x <- subset(x, select = keep)
  return(x)
}

# Normalize a data frame or matrix, using the 1/n variance formula
PopulationScale <- function(original)
{
  means <- apply(original, 2, mean)
  centered <- sweep(original, 2, means, "-", check.margin = FALSE)
  sds <- apply(centered, 2, function(x) sqrt(sum(x^2) / length(x)))
  scaled <- sweep(centered, 2, sds, "/", check.margin = FALSE)
  return(scaled)
}

# Load expression and genotype data, adjust for covariates, and normalize
CovarAdjust <- function(plink_path, covar_path)
{
  # Load genotype and expression data
  data <- ReadPlink(plink_path)
  genotypes <- data$genotypes
  expression <- data$expression
  
  # Load covariates and concatenate IDs to match with expression data
  covars <- read.table(covar_path, header = TRUE, stringsAsFactors = FALSE)
  rownames(covars) <- paste(covars$FID, covars$IID, sep = "_")
  covars <- subset(covars, select = -c(FID, IID))
  covars <- RemoveSuperfluous(covars)
  
  # Subset to intersection of individuals for whom genotype/expression and covariate data both exist
  expression_and_covars <- merge(expression, covars, by = 0)
  
  # Adjust expression values for covariates
  expression <- summary(lm(value ~ ., data = expression_and_covars))$residuals
  colnames(expression) <- "value"
  rownames(expression) <- rownames(expression_and_covars)
  
  # Normalize the genotypes and the adjusted expression values
  genotypes <- PopulationScale(genotypes)
  expression <- PopulationScale(expression)
  
  # Remove monomorphic SNPs
  genotypes <- RemoveSuperfluous(genotypes)
  
  # Subset to intersection of individuals for whom genotype/expression and covariate data both exist
  genotypes_and_covars <- merge(genotypes, covars, by = 0)
  
  # Adjust normalized genotypes for covariates
  for (snp in colnames(genotypes)) {
    genotypes_and_covars[[snp]] <- summary(lm(reformulate(colnames(covars), response = snp), data = genotypes_and_covars))$residuals
  }
  
  # Remove covars from the merged data frame
  genotypes <- subset(genotypes_and_covars, select = colnames(genotypes))
  
  # Normalize and remove monomorphics again
  genotypes <- PopulationScale(genotypes)
  genotypes <- RemoveSuperfluous(genotypes)
  
  return(list(genotypes = genotypes, expression = expression))
}

# Train an elastic net model to predict gene expression, using BIC to select the complexity parameter lambda
ElasticNetBIC <- function(genos, pheno, alpha = 0.5)
{
  enet_model <- glmnet::glmnet(x = as.matrix(genos), y = as.matrix(pheno),
                                family = "gaussian", alpha = alpha, intercept = FALSE)
  
  # Here we calculate BIC according to a modified formula (see comments below)
  eff_df <- NULL
  for (l in enet_model$lambda) {
    
    nonzero_indices <- predict(enet_model, s = l, type = "nonzero")
    
    # If all coefficients are zero, skip to the next lambda
    if (is.null(nonzero_indices$s1)) {
      eff_df <- c(eff_df, 0)
      next
    }
    
    # Subset the design matrix to those columns that correspond to nonzero coefficients
    x <- as.matrix(subset(genos, select = nonzero_indices$s1))
    
    # Effective degrees of freedom for elastic net, formula from http://dx.doi.org/10.1214/12-AOS1003
    eff_df <- c(eff_df, sum(diag(x %*% tcrossprod(chol2inv(chol(crossprod(x, x) + diag(l*(1-alpha), nrow = ncol(x)))), x))))
  }
  
  # This gives 2*(loglikelihood - loglikelihood(null)) instead of the actual loglikelihood,
  # But that's fine for model selection since the log-likelihood of the null model does not depend on lambda.
  biased_lik <- enet_model$nulldev - deviance(enet_model)
  biased_bic <- eff_df * log(enet_model$nobs) - biased_lik
  
  # Extract the lambda and df values corresponding to lowest BIC
  best_index <- which(biased_bic == min(biased_bic))
  best_lambda <- enet_model$lambda[best_index]
  
  return(list(model = enet_model, lambda = best_lambda))
}

# Compute per-observation log-likelihoods for an elastic net model on the specified data
GetPredictions <- function(fit, genotypes, expression)
{
  expression_predicted <- predict(fit$model, newx = as.matrix(genotypes), s = fit$lambda)
  sigma <- sqrt(sum((expression - expression_predicted)^2) / length(expression))
  log_liks <- log(dnorm(x = expression, mean = expression_predicted, sd = sigma))
  
  return(log_liks)
}

# Compute the likelihood-ratio test statistic.
# Namely, we use the statistic described in Section 5 of Vuong (1989), accessible at https://www.jstor.org/stable/1912557.
LRT <- function(ll_A, ll_B)
{
  # Eq. 3.1 in Vuong (1989)
  LR <- sum(ll_A - ll_B)
  
  # Eq. 4.2 in Vuong (1989)
  variance <- (1 / length(ll_A)) * sum((ll_A - ll_B)^2) - ((1 / length(ll_A)) * sum(ll_A - ll_B))^2
  
  # Eq. 5.6 in Vuong (1989)
  lrt_stat <- LR / sqrt(length(ll_A) * variance)
  
  return(lrt_stat)
}

# Compute the distribution-free test statistic.
# Namely, we use the statistic described in Section 2.2 of Clarke (2007), accessible at https://doi.org/10.1093/pan/mpm004.
DFT <- function(ll_A, ll_B)
{
  # Eq. 7 in Clarke (2007)
  d <- ll_A - ll_B
  B <- sum(d > 0)
  
  return(B)
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

name <- args[1]
id <- args[2]
tissue_A <- args[3]
tissue_B <- args[4]
job <- args[5]
out_name <- args[6]
replicates <- args[7]

# Import expression and genotype data, adjust them for covariates, and normalize them
data_A_train <- CovarAdjust(paste("temp/", job, "/", name, "/", tissue_A, "_part1", sep = ""), paste("temp/", job, "/", tissue_A, "_part1.expression_covariates.txt", sep = ""))
data_A_test <- CovarAdjust(paste("temp/", job, "/", name, "/", tissue_A, "_part2", sep = ""), paste("temp/", job, "/", tissue_A, "_part2.expression_covariates.txt", sep = ""))
data_B <- CovarAdjust(paste("temp/", job, "/", name, "/", tissue_B, sep = ""), paste("covariates/", tissue_B, ".expression_covariates.txt", sep = ""))

# Train elastic net models for each tissue
elnet_A <- ElasticNetBIC(data_A_train$genotypes, data_A_train$expression)
elnet_B <- ElasticNetBIC(data_B$genotypes, data_B$expression)

# Compute prediction log-likelihoods on data_A_test for model trained on data_A_train
likelihoods_A <- GetPredictions(elnet_A, data_A_test$genotypes, data_A_test$expression)

# Compute prediction log-likelihoods on data_A_test for model trained on data_B
likelihoods_B <- GetPredictions(elnet_B, data_A_test$genotypes, data_A_test$expression)

# The prediction log-likelihoods should have the same length (just a sanity check)
if (length(likelihoods_A) != length(likelihoods_B)) {
  cat(paste("ERROR: Testing data is inconsistent. Skipping gene", id, ".\n"), file = stderr())
  q()
}

# Calculate the test statistics
lrt_stat <- LRT(likelihoods_A, likelihoods_B)
dft_stat <- DFT(likelihoods_A, likelihoods_B)

# Approximate the distribution of each statistic under the null of no tissue-specific differences in regulation.
# Namely, we perform a permutation test by independently shuffling the prediction log-likelihoods for each individual.

likelihoods <- cbind(likelihoods_A, likelihoods_B)
lrt_stat_permutations <- numeric(replicates)
dft_stat_permutations <- numeric(replicates)

for (i in seq_len(replicates)) {
  likelihoods_resampled <- t(apply(likelihoods, 1, function(x) { sample(x) }))
  lrt_stat_permutations[i] <- LRT(likelihoods_resampled[, 1], likelihoods_resampled[, 2])
  dft_stat_permutations[i] <- DFT(likelihoods_resampled[, 1], likelihoods_resampled[, 2])
}

lrt_pval <- sum(abs(lrt_stat_permutations) >= abs(lrt_stat)) / replicates

dft_stat_translated <- 2 * mean(dft_stat_permutations) - dft_stat
if (dft_stat >= mean(dft_stat_permutations)) {
  dft_pval <- (sum(dft_stat_permutations >= dft_stat) + sum(dft_stat_permutations <= dft_stat_translated)) / replicates
} else {
  dft_pval <- (sum(dft_stat_permutations <= dft_stat) + sum(dft_stat_permutations >= dft_stat_translated)) / replicates
}

# Append the test results to the output file
results <- paste(name, id, lrt_pval, dft_pval, sep = "\t")
cat(results, file = paste("output/", tissue_A, "-", tissue_B, "-", out_name, ".txt", sep = ""), append = TRUE, sep = "\n")

