tissue <- commandArgs(trailingOnly = TRUE)[1]
if (commandArgs(trailingOnly = TRUE)[2] == "SAVE") {
  save <- TRUE
} else {
  save <- FALSE
}

# Read list of individuals for the given tissue
samples <- read.table(paste("expression_covariates/", tissue, ".v8.EUR.covariates.txt.HEADER", sep = ""))
samples <- samples[, 1]

# Create 1/2, 1/2 splits of the individuals
samples_half1 <- sample(samples, size = floor(0.5 * length(samples)))
samples_half2 <- setdiff(samples, samples_half1)

# If creating splits for null testing, then also create 1/3, 2/3 splits of the individuals
if (save) {
  samples_third <- sample(samples, size = floor((1/3) * length(samples)))
  samples_twothirds <- setdiff(samples, samples_third)
}

# Read expression and covariate matrices
expression_matrix <- read.table(paste("expression_matrices/", tissue, ".v8.EUR.normalized_expression.bed.gz", sep = ""), header = TRUE, check.names = FALSE, comment.char = "")
expression_covariates <- read.table(paste("expression_covariates/", tissue, ".v8.EUR.covariates.plink.txt", sep = ""), header = TRUE, check.names = FALSE, comment.char = "")

# Split the expression matrix into the block that contains annotations and the block that contains expression data
annotations <- expression_matrix[c(1,2,3,4)]
expression <- expression_matrix[-c(1,2,3,4)]

# Split the expression data into halves
expression_half1 <- expression[colnames(expression) %in% samples_half1]
expression_half2 <- expression[colnames(expression) %in% samples_half2]
# If creating splits for null testing, then also split the expression data into 1/3, 2/3 partitions
if (save) {
  expression_third <- expression[colnames(expression) %in% samples_third]
  expression_twothirds <- expression[colnames(expression) %in% samples_twothirds]
}

# Split the covariate data into halves
covariates_half1 <- expression_covariates[expression_covariates$IID %in% samples_half1, ]
covariates_half2 <- expression_covariates[expression_covariates$IID %in% samples_half2, ]
# If creating splits for null testing, then also split the covariate data into 1/3, 2/3 partitions
if (save) {
  covariates_third <- expression_covariates[expression_covariates$IID %in% samples_third, ]
  covariates_twothirds <- expression_covariates[expression_covariates$IID %in% samples_twothirds, ]
}

if (!save) {
  # Write the split expression data to new files in a temporary location
  write.table(cbind(annotations, expression_half1), file = paste("temp/", tissue, "_half1.v8.EUR.normalized_expression.bed", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(cbind(annotations, expression_half2), file = paste("temp/", tissue, "_half2.v8.EUR.normalized_expression.bed", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE)
} else {
  # Write the split expression data to new files in a permanent location
  write.table(cbind(annotations, expression_half1), file = paste("expression_matrices/", tissue, "_half1.v8.EUR.normalized_expression.bed", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(cbind(annotations, expression_half2), file = paste("expression_matrices/", tissue, "_half2.v8.EUR.normalized_expression.bed", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(cbind(annotations, expression_third), file = paste("expression_matrices/", tissue, "_third.v8.EUR.normalized_expression.bed", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(cbind(annotations, expression_twothirds), file = paste("expression_matrices/", tissue, "_twothirds.v8.EUR.normalized_expression.bed", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE)
}

if (!save) {
  # Write the split covariate data to new files in a temporary location
  write.table(covariates_half1, file = paste("temp/", tissue, "_half1.v8.EUR.covariates.plink.txt", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(covariates_half2, file = paste("temp/", tissue, "_half2.v8.EUR.covariates.plink.txt", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE)
} else {
  # Write the split covariate data to new files in a permanent location
  write.table(covariates_half1, file = paste("expression_covariates/", tissue, "_half1.v8.EUR.covariates.plink.txt", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(covariates_half2, file = paste("expression_covariates/", tissue, "_half2.v8.EUR.covariates.plink.txt", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(covariates_third, file = paste("expression_covariates/", tissue, "_third.v8.EUR.covariates.plink.txt", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(covariates_twothirds, file = paste("expression_covariates/", tissue, "_twothirds.v8.EUR.covariates.plink.txt", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE)
}

if (!save) {
  # Write the lists of individuals in each split to new header files in a temporary location
  write.table(samples_half1, file = paste("temp/", tissue, "_half1.v8.EUR.covariates.txt.HEADER", sep = ""), sep = "\n", quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(samples_half2, file = paste("temp/", tissue, "_half2.v8.EUR.covariates.txt.HEADER", sep = ""), sep = "\n", quote = FALSE, row.names = FALSE, col.names = FALSE)
} else {
  # Write the lists of individuals in each split to new header files in a permanent location
  write.table(samples_half1, file = paste("expression_covariates/", tissue, "_half1.v8.EUR.covariates.txt.HEADER", sep = ""), sep = "\n", quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(samples_half2, file = paste("expression_covariates/", tissue, "_half2.v8.EUR.covariates.txt.HEADER", sep = ""), sep = "\n", quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(samples_third, file = paste("expression_covariates/", tissue, "_third.v8.EUR.covariates.txt.HEADER", sep = ""), sep = "\n", quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(samples_twothirds, file = paste("expression_covariates/", tissue, "_twothirds.v8.EUR.covariates.txt.HEADER", sep = ""), sep = "\n", quote = FALSE, row.names = FALSE, col.names = FALSE)
}

