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
  # Individual names are coded as FID_IID
  # Variant names are coded as rsid_A1
  rownames(bed) <- paste(fam$V1, fam$V2, sep = "_")
  colnames(bed) <- paste(bim$V2, bim$V5, sep = "_")
  
  # Set row/column names and remove unused columns in the imported fam
  rownames(fam) <- paste(fam$V1, fam$V2, sep = "_")
  fam <- subset(fam, select = V6)
  colnames(fam) <- "value"
  
  return(list(genotypes = bed, expression = fam))
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

# Remove columns with NAs or a standard deviation of zero
RemoveSuperfluous <- function(x)
{
  sds <- apply(x, 2, sd, na.rm = FALSE)
  keep <- (sds != 0) & (is.finite(sds))
  x <- subset(x, select = keep)
  return(x)
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
  
  # This gives 2*(loglikelihood - loglikelihood(null)) instead of the actual loglikelihood
  # But that's fine for model selection since the log-likelihood of the null model does not depend on lambda
  biased_lik <- enet_model$nulldev - deviance(enet_model)
  biased_bic <- eff_df * log(enet_model$nobs) - biased_lik
  
  # Extract the lambda and df values corresponding to lowest BIC
  best_index <- which(biased_bic == min(biased_bic))
  best_lambda <- enet_model$lambda[best_index]
  best_df <- eff_df[best_index]
  
  return(list(models = enet_model, lambda = best_lambda, df = best_df))
}

# Compute per-observation log-likelihoods for an elastic net model on the specified data
GetLogLik <- function(fit, genotypes, expression)
{
  expression_predicted <- predict(fit$models, newx = as.matrix(genotypes), s = fit$lambda)
  sigma <- sqrt(sum((expression - expression_predicted)^2) / length(expression))
  log_liks <- log(dnorm(x = expression, mean = expression_predicted, sd = sigma))
  
  return(list(log_liks = log_liks, df = fit$df, n = length(expression)))
}

# @@ To-Do: finish this function, following the nonnest2 package
# Calculate estimates of A and B, as defined in Eq. 2.1 and 2.2 of Vuong (1989)
calcAB <- function(loglik, n) {
  # Eq. 2.1
  scaling <- summary(object)$sigma
  
  if(is.null(scaling)){
    scaling <- 1
  } else {
    scaling <- scaling^2
  }
  tmpvc <- n * vc(object)
  
  A <- chol2inv(chol(tmpvc))
  
  ## Eq (2.2)
  if(!is.null(scfun)){
    sc <- scfun(object)
  } else if(class(object)[1] == "lavaan"){
    sc <- estfun(object, remove.duplicated=TRUE)
  } else if(class(object)[1] %in% c("SingleGroupClass", "MultipleGroupClass")){
    wts <- mirt::extract.mirt(object, "survey.weights")
    if(length(wts) > 0){
      sc <- mirt::estfun.AllModelClass(object, weights = sqrt(wts))
    } else {
      sc <- mirt::estfun.AllModelClass(object)
    }
  } else if(class(object)[1] %in% c("lm", "glm", "nls")){
    sc <- (1/scaling) * estfun(object)
  } else {
    sc <- estfun(object)
  }
  sc.cp <- crossprod(sc)/n
  B <- matrix(sc.cp, nrow(A), nrow(A))
  
  list(A=A, B=B, sc=sc)
}

# @@ To-Do: finish this function, following the nonnest2 package
# Calculate the eigenvalues of the W matrix defined in Eq. 3.6 of Vuong (1989)
calcLambda <- function(logliks_A, logliks_B, n) {
  AB1 <- CalcAB(logliks_A, n)
  AB2 <- CalcAB(logliks_B, n)
  Bc <- CalcBcross(AB1$sc, AB2$sc, n)
  
  W <- cbind(rbind(-AB1$B %*% chol2inv(chol(AB1$A)),
                   t(Bc) %*% chol2inv(chol(AB1$A))),
             rbind(-Bc %*% chol2inv(chol(AB2$A)),
                   AB2$B %*% chol2inv(chol(AB2$A))))
  
  lambda_star <- eigen(W, only.values = TRUE)$values
  
  # Sometimes the eigenvalues are complex, so we only return the real parts
  Return(Re(lambda_star))
}

# @@ To-Do: finish this function, following the nonnest2 package
# Calculate the Vuong test p-values
# This is the two-step test described in Section 6 of Vuong (1989), accessible at https://www.jstor.org/stable/1912557
# We report both p-values, which users can then use to conduct an overall test at their desired significance level
LRT <- function(A, B, correction = "none")
{
  # Step 1: the variance test -- this p-value should be considered first
  
  # Eq. 4.2 in Vuong (1989)
  variance <- (1 / A$n) * sum((A$log_liks - B$log_liks)^2) - ((1 / A$n) * sum(A$log_liks - B$log_liks))^2
  
  # Get p-value of weighted chi-square distribution
  lambda_star <- CalcLambda(A$log_liks, B$log_liks, A$n)
  var_p <- CompQuadForm::imhof(A$n * variance, lambda_star^2)$Qq
  
  # Step 2: likelihood-ratio test -- if the null hypothesis in the variance test is rejected, then consider this p-value
  
  # Eq. 3.1 in Vuong (1989)
  LR <- sum(A$log_liks - B$log_liks)
  
  # AIC-based correction is applied to both log likelihoods
  if (correction == "both")
    LR <- LR - (A$df - B$df)
  # AIC-based correction is applied to the model A log likelihood
  if (correction == "A")
    LR <- LR - A$df
  # AIC-based correction is applied to the model B log likelihood
  if (correction == "B")
    LR <- LR + B$df
  
  # Eq. 5.6 in Vuong (1989)
  lrt_stat <- LR / sqrt(A$n * variance)
  lrt_p <- 2 * pnorm(-abs(statistic))
  
  return(list(var_p = var_p, lrt_p = lrt_p))
}

# @@ To-Do: rewrite this function (it's still the pre-bootstrap version)
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
if (length(args) != 6) {
  cat("ERROR: missing arguments.\n", file = stderr())
  q()
}

name <- args[1]
id <- args[2]
tissue_A <- args[3]
tissue_B <- args[4]
job <- args[5]
out_name <- args[6]

# Import expression and genotype data, adjust it for covariates, and normalize it
data_A_train <- CovarAdjust(paste("temp/", job, "/", tissue_A, "_part1", sep = ""), paste("temp/", job, "/", tissue_A, "_part1.expression_covariates.txt", sep = ""))
data_A_test <- CovarAdjust(paste("temp/", job, "/", tissue_A, "_part2", sep = ""), paste("temp/", job, "/", tissue_A, "_part2.expression_covariates.txt", sep = ""))
data_B <- CovarAdjust(paste("expression/", tissue_B, sep = ""), paste("covariates/", tissue_B, ".expression_covariates.txt", sep = ""))

# Train elastic net models for each tissue
elnet_A <- ElasticNetBIC(data_A_train$genotypes, data_A_train$expression)
elnet_B <- ElasticNetBIC(data_B$genotypes, data_B$expression)

# Compute prediction log-likelihoods on data_A_test for model trained on data_A_train
likelihoods_A <- GetLogLik(elnet_A, data_A_test$genotypes, data_A_test$expression)

# Compute prediction log-likelihoods on data_A_test for model trained on data_B
likelihoods_B <- GetLogLik(elnet_B, data_A_test$genotypes, data_A_test$expression)

# Calculate p-values for each model selection test
if (likelihoods_A$n != likelihoods_B$n) {
  cat(paste("ERROR: Testing data is inconsistent. Skipping gene", id, ".\n"), file = stderr())
} else {
  lr_pval <- LRT(likelihoods_A, likelihoods_B)
  df_pval <- DFT(likelihoods_A, likelihoods_B)
}
# Append the test results to the output file
results <- paste(name, id, lr_pval, df_pval, sep = "\t")
cat(results, file = paste("output/", tissue_A, "-", tissue_B, "-", out_name, ".txt", sep = ""), append = TRUE, sep = "\n")

