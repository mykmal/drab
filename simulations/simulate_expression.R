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
  bim <- read.table(bimfile, header = FALSE, stringsAsFactors = FALSE)
  fam <- read.table(famfile, header = FALSE, stringsAsFactors = FALSE)
  # Note that specifying n and p speeds up BEDMatrix, since otherwise it would read the fam and bim files again
  bed <- as.data.frame(as.matrix(BEDMatrix::BEDMatrix(bedfile, n = nrow(fam), p = nrow(bim))))
  
  # Add individual and variant names to the imported bed:
  # Individual names are coded as FID_IID
  # Variant names are coded the same as in the original plink file
  rownames(bed) <- paste(fam$V1, fam$V2, sep = "_")
  colnames(bed) <- bim$V2
  
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

# Load expression and genotype data, normalize both, and adjust expression for covariates
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
  
  return(list(genotypes = genotypes, expression = expression))
}

# Train a transcriptome imputation model, saving its weights and returning the model object
TrainSave <- function(training_genotypes, training_expression, context)
{
  # Fit an elastic net regression model on the training data
  model <- cv.glmnet(x = as.matrix(training_genotypes), y = as.matrix(training_expression),
                                  family = "gaussian", alpha = 0.5, nfolds = 5, type.measure = "mse", standardize = FALSE)
  
  # Extract the fitted model weights
  wgt_matrix <- as.matrix(coef(model, s = "lambda.1se"))
  colnames(wgt_matrix) <- "weight"
  
  # Save the model weights
  write.table(wgt_matrix, file = paste0("saved_models/", context, "_", name, "_", id, ".weights"),
                                  sep = "\t", quote = FALSE)
  
  return(model)
}


######################################################################
# MAIN PROGRAM
######################################################################

suppressMessages(library(glmnet))

args <- commandArgs(trailingOnly = TRUE)

job <- args[1]
context_A <- args[2]
context_B <- args[3]
name <- args[4]
id <- args[5]

# Import context-specific data and adjust expression for covariates
data_train_A <- CovarAdjust(paste0(job, "/", name, "/A.genotypes"), paste0("covariates/", context_A, ".covariates.txt"))
data_train_B <- CovarAdjust(paste0(job, "/", name, "/B.genotypes"), paste0("covariates/", context_B, ".covariates.txt"))

# Load genotype data for all samples
full_genotypes <- ReadPlink(paste0(job, "/", name, "/full.genotypes"))$genotypes

# Scale using the 1/n variance formula and remove monomorphic SNPs
n <- nrow(full_genotypes)
full_genotypes <- as.data.frame((scale(full_genotypes) / sqrt(n - 1)) * sqrt(n))
full_genotypes <- RemoveSuperfluous(full_genotypes)

# Subset to a common set of SNPs
all_snps <- Reduce(intersect, list(colnames(data_train_A$genotypes),
                                   colnames(data_train_B$genotypes),
                                   colnames(full_genotypes)))
data_train_A$genotypes <- subset(data_train_A$genotypes, select = all_snps)
data_train_B$genotypes <- subset(data_train_B$genotypes, select = all_snps)
full_genotypes <- subset(full_genotypes, select = all_snps)
full_genotypes <- as.matrix(full_genotypes)

# Train and save context-specific models using all available data for each context
trained_model_A <- TrainSave(data_train_A$genotypes, data_train_A$expression$value, context_A)
trained_model_B <- TrainSave(data_train_B$genotypes, data_train_B$expression$value, context_B)

# Estimate the amount of variance unexplained by genetics
prediction_A <- predict(trained_model_A,
                        newx = as.matrix(data_train_A$genotypes),
                        s = "lambda.1se",
                        type = "response")
unexplained_sd_A <- sd(data_train_A$expression$value - prediction_A[, 1])

prediction_B <- predict(trained_model_B,
                        newx = as.matrix(data_train_B$genotypes),
                        s = "lambda.1se",
                        type = "response")
unexplained_sd_B <- sd(data_train_B$expression$value - prediction_B[, 1])

# Using weights from each of the trained models, simulate expression data for all GTEx samples
expression_imputed_A <- predict(trained_model_A,
                                newx = full_genotypes,
                                s = "lambda.1se",
                                type = "response")
expression_imputed_A <- expression_imputed_A + rnorm(length(expression_imputed_A), mean = 0, sd = unexplained_sd_A)
expression_imputed_B <- predict(trained_model_B,
                                newx = full_genotypes,
                                s = "lambda.1se",
                                type = "response")
expression_imputed_B <- expression_imputed_B + rnorm(length(expression_imputed_B), mean = 0, sd = unexplained_sd_B)

# Save the simulated expression values
write.table(expression_imputed_A, file = paste0("expression_simulated/", context_A, "_", name, "_", id, "_expression.simulated.txt"), sep = "\t", quote = FALSE)
write.table(expression_imputed_B, file = paste0("expression_simulated/", context_B, "_", name, "_", id, "_expression.simulated.txt"), sep = "\t", quote = FALSE)

