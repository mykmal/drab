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
  
  return(list(genotypes = bed, expression = fam, n = nrow(fam)))
}

# Remove columns with NAs or a standard deviation of zero
RemoveSuperfluous <- function(x)
{
  sds <- apply(x, 2, sd)
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
  n <- data$n
  
  # Load covariates and concatenate IDs to match with expression data
  covars <- read.table(covar_path, header = TRUE, stringsAsFactors = FALSE)
  rownames(covars) <- paste(covars$FID, covars$IID, sep = "_")
  covars <- subset(covars, select = -c(FID, IID))
  covars <- RemoveSuperfluous(covars)
  
  # Adjust expression values for covariates
  adjust_expression_formula <- reformulate(paste("covars$", colnames(covars), sep = ""), response = "expression$value")
  expression <- as.data.frame(resid(lm(adjust_expression_formula)))
  colnames(expression) <- "value"
  
  # Normalize the genotypes and the adjusted expression values
  # Note that we use the 1/n variance formula to match what glmnet expects
  genotypes <- as.data.frame((scale(genotypes) / sqrt(n - 1)) * sqrt(n))
  expression <- as.data.frame((scale(expression) / sqrt(n - 1)) * sqrt(n))
  
  # Remove monomorphic SNPs
  genotypes <- RemoveSuperfluous(genotypes)
  
  # Adjust normalized genotypes for covariates
  for (snp in colnames(genotypes)) {
    adjust_genotypes_formula <- reformulate(paste("covars$", colnames(covars), sep = ""), response = paste("genotypes$", snp, sep = ""))
    genotypes[[snp]] <- resid(lm(adjust_genotypes_formula))
  }
  
  # Normalize and remove monomorphics again
  genotypes <- as.data.frame((scale(genotypes) / sqrt(n - 1)) * sqrt(n))
  genotypes <- RemoveSuperfluous(genotypes)
  
  return(list(genotypes = genotypes, expression = expression))
}

# Train and test an expression imputation model, returning the squared prediction residuals
TrainTest <- function(training_genotypes, training_expression, testing_genotypes, testing_expression)
{
  # Fit an elastic net regression model on the training data
  model <- glmnet::cv.glmnet(x = as.matrix(training_genotypes), y = as.matrix(training_expression),
                                  family = "gaussian", alpha = 0.5, nfolds = 5, type.measure = "mse", standardize = FALSE)
  
  # Get predictions for the model on the testing data
  predictions <- predict(model, newx = as.matrix(testing_genotypes), type = "response", s = "lambda.1se")
  
  # Compute per-individual squared residuals
  loss <- (testing_expression - predictions[,1])^2
  
  return(loss)
}


######################################################################
# MAIN PROGRAM
######################################################################

args <- commandArgs(trailingOnly = TRUE)

job <- args[1]
boot <- as.numeric(args[2])
context_A <- args[3]
context_B <- args[4]
out_name <- args[5]
name <- args[6]
id <- args[7]

# Import data and adjust for covariates
data_train_1 <- CovarAdjust(paste(job, "/", name, "/part1", sep = ""), paste(job, "/part1.covariates.txt", sep = ""))
data_train_2 <- CovarAdjust(paste(job, "/", name, "/part2", sep = ""), paste(job, "/part2.covariates.txt", sep = ""))
data_test <- CovarAdjust(paste(job, "/", name, "/part3", sep = ""), paste(job, "/part3.covariates.txt", sep = ""))

# Subset the three data sets to a common set of SNPs
all_snps <- Reduce(intersect, list(colnames(data_train_1$genotypes), colnames(data_train_2$genotypes), colnames(data_test$genotypes)))
data_train_1$genotypes <- subset(data_train_1$genotypes, select = all_snps)
data_train_2$genotypes <- subset(data_train_2$genotypes, select = all_snps)
data_test$genotypes <- subset(data_test$genotypes, select = all_snps)

# Find the sample mean and variance of the differences between model prediction errors on the full dataset
loss_full_1 <- TrainTest(data_train_1$genotypes, data_train_1$expression$value, data_test$genotypes, data_test$expression$value)
loss_full_2 <- TrainTest(data_train_2$genotypes, data_train_2$expression$value, data_test$genotypes, data_test$expression$value)
var_differences_full <- var(loss_full_1 - loss_full_2)
mean_differences_full <- mean(loss_full_1 - loss_full_2)

# Find the distribution of sample means of the differences between model prediction errors using bootstrapping
boot_means <- numeric(boot)
for (i in 1:boot) {
  # Sample from the rows of each dataset with replacement
  resamples_train_1 <- sample(seq_len(nrow(data_train_1$expression)), replace = TRUE)
  resamples_train_2 <- sample(seq_len(nrow(data_train_2$expression)), replace = TRUE)
  resamples_test <- sample(seq_len(nrow(data_test$expression)), replace = TRUE)
  
  # Find the sample mean of the differences between model prediction errors on the resampled dataset
  loss_boot_1 <- TrainTest(data_train_1$genotypes[resamples_train_1, ], data_train_1$expression[resamples_train_1, ], data_test$genotypes[resamples_test, ], data_test$expression[resamples_test, ])
  loss_boot_2 <- TrainTest(data_train_2$genotypes[resamples_train_2, ], data_train_2$expression[resamples_train_2, ], data_test$genotypes[resamples_test, ], data_test$expression[resamples_test, ])
  boot_means[i] <- mean(loss_boot_1 - loss_boot_2)
}

# Conduct the DRAB test, comparing whether the models have equal predictive performance
n <- nrow(data_test$expression)
t <- mean_differences_full / sqrt((var_differences_full / n) + var(boot_means))
pval <- 2 * pt(abs(t), df = n - 1, lower.tail = FALSE)
pval_conditional <- t.test(x = loss_full_1, y = loss_full_2, paired = TRUE)$p.value

# Append the test results to the output file
results <- paste(name, id, pval, pval_conditional, nrow(data_train_1$expression), nrow(data_test$expression), sep = "\t")
cat(results, file = paste("output/", context_A, "-", context_B, "-", out_name, ".txt", sep = ""), append = TRUE, sep = "\n")
