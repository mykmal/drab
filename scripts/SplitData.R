args <- commandArgs(trailingOnly = TRUE)

tissue <- args[1]
out_path <- args[2]

split_size <- 0.5

# Read expression and covariate data
expression_matrix <- read.table(paste("expression/", tissue, ".expression_matrix.txt", sep = ""), header = TRUE)
expression_covariates <- read.table(paste("covariates/", tissue, ".expression_covariates.txt", sep = ""), header = TRUE)

# Merge FID and IID, since both are necessary to uniquely identify an individual
expression_matrix$ID <- paste(expression_matrix$FID, expression_matrix$IID, sep = "_")
expression_covariates$ID <- paste(expression_covariates$FID, expression_covariates$IID, sep = "_")
individuals <- intersect(expression_matrix$ID, expression_covariates$ID)

# Randomly partition the individuals
indivs_part1 <- sample(individuals, size = floor(split_size * length(individuals)), replace = FALSE)
indivs_part2 <- setdiff(individuals, indivs_part1)

# Split the expression and covariate data
expression_part1 <- expression_matrix[expression_matrix$ID %in% indivs_part1, ]
expression_part2 <- expression_matrix[expression_matrix$ID %in% indivs_part2, ]
covariates_part1 <- expression_covariates[expression_covariates$ID %in% indivs_part1, ]
covariates_part2 <- expression_covariates[expression_covariates$ID %in% indivs_part2, ]

# Remove the ID column to retain the original data format
expression_part1 <- subset(expression_part1, select = -ID)
expression_part2 <- subset(expression_part2, select = -ID)
covariates_part1 <- subset(covariates_part1, select = -ID)
covariates_part2 <- subset(covariates_part2, select = -ID)

# Write the split data to new files in the specified location
write.table(expression_part1, file = paste(out_path, "/", tissue, "_part1.expression_matrix.txt", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(expression_part2, file = paste(out_path, "/", tissue, "_part2.expression_matrix.txt", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(covariates_part1, file = paste(out_path, "/", tissue, "_part1.expression_covariates.txt", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(covariates_part2, file = paste(out_path, "/", tissue, "_part2.expression_covariates.txt", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE)

