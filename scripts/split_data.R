args <- commandArgs(trailingOnly = TRUE)

tissue_A <- args[1]
tissue_B <- args[2]
out_path <- args[3]

if (tissue_A == tissue_B) {
  
  tissue <- tissue_A
  
  split_size <- 1 / 3
  
  expression_matrix <- read.table(paste("expression/", tissue, ".expression_matrix.txt", sep = ""), header = TRUE)
  expression_covariates <- read.table(paste("covariates/", tissue, ".expression_covariates.txt", sep = ""), header = TRUE)
  
  expression_matrix$ID <- paste(expression_matrix$FID, expression_matrix$IID, sep = "_")
  expression_covariates$ID <- paste(expression_covariates$FID, expression_covariates$IID, sep = "_")
  individuals <- intersect(expression_matrix$ID, expression_covariates$ID)
  
  indivs_part1 <- sample(individuals, size = ceiling(split_size * length(individuals)), replace = FALSE)
  indivs_part2 <- sample(setdiff(individuals, indivs_part1), size = length(indivs_part1), replace = FALSE)
  indivs_part3 <- setdiff(individuals, union(indivs_part1, indivs_part2))
  
  expression_part1 <- expression_matrix[expression_matrix$ID %in% indivs_part1, ]
  expression_part2 <- expression_matrix[expression_matrix$ID %in% indivs_part2, ]
  expression_part3 <- expression_matrix[expression_matrix$ID %in% indivs_part3, ]
  covariates_part1 <- expression_covariates[expression_covariates$ID %in% indivs_part1, ]
  covariates_part2 <- expression_covariates[expression_covariates$ID %in% indivs_part2, ]
  covariates_part3 <- expression_covariates[expression_covariates$ID %in% indivs_part3, ]
  
} else {
  
  expression_matrix_A <- read.table(paste("expression/", tissue_A, ".expression_matrix.txt", sep = ""), header = TRUE)
  expression_covariates_A <- read.table(paste("covariates/", tissue_A, ".expression_covariates.txt", sep = ""), header = TRUE)
  expression_matrix_B <- read.table(paste("expression/", tissue_B, ".expression_matrix.txt", sep = ""), header = TRUE)
  expression_covariates_B <- read.table(paste("covariates/", tissue_B, ".expression_covariates.txt", sep = ""), header = TRUE)
  
  expression_matrix_A$ID <- paste(expression_matrix_A$FID, expression_matrix_A$IID, sep = "_")
  expression_covariates_A$ID <- paste(expression_covariates_A$FID, expression_covariates_A$IID, sep = "_")
  individuals_A <- intersect(expression_matrix_A$ID, expression_covariates_A$ID)
  
  expression_matrix_B$ID <- paste(expression_matrix_B$FID, expression_matrix_B$IID, sep = "_")
  expression_covariates_B$ID <- paste(expression_covariates_B$FID, expression_covariates_B$IID, sep = "_")
  individuals_B <- intersect(expression_matrix_B$ID, expression_covariates_B$ID)
  
  if (length(individuals_A) < length(individuals_B)) {
    
    indivs_part1 <- individuals_A
    mutual <- intersect(individuals_B, indivs_part1)
    indivs_part2 <- union(mutual, sample(setdiff(individuals_B, indivs_part1), size = (length(indivs_part1) - length(mutual)), replace = FALSE))
    indivs_part3 <- setdiff(individuals_B, indivs_part2)
    
    expression_part3 <- expression_matrix_B[expression_matrix_B$ID %in% indivs_part3, ]
    covariates_part3 <- expression_covariates_B[expression_covariates_B$ID %in% indivs_part3, ]
    
  } else {
    
    indivs_part2 <- individuals_B
    mutual <- intersect(individuals_A, indivs_part2)
    indivs_part1 <- union(mutual, sample(setdiff(individuals_A, indivs_part2), size = (length(indivs_part2) - length(mutual)), replace = FALSE))
    indivs_part3 <- setdiff(individuals_A, indivs_part1)
    
    expression_part3 <- expression_matrix_A[expression_matrix_A$ID %in% indivs_part3, ]
    covariates_part3 <- expression_covariates_A[expression_covariates_A$ID %in% indivs_part3, ]
  }
  
  expression_part1 <- expression_matrix_A[expression_matrix_A$ID %in% indivs_part1, ]
  expression_part2 <- expression_matrix_B[expression_matrix_B$ID %in% indivs_part2, ]
  covariates_part1 <- expression_covariates_A[expression_covariates_A$ID %in% indivs_part1, ]
  covariates_part2 <- expression_covariates_B[expression_covariates_B$ID %in% indivs_part2, ]
}

expression_part1 <- subset(expression_part1, select = -ID)
expression_part2 <- subset(expression_part2, select = -ID)
expression_part3 <- subset(expression_part3, select = -ID)
covariates_part1 <- subset(covariates_part1, select = -ID)
covariates_part2 <- subset(covariates_part2, select = -ID)
covariates_part3 <- subset(covariates_part3, select = -ID)

write.table(expression_part1, file = paste(out_path, "/part1.expression_matrix.txt", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(expression_part2, file = paste(out_path, "/part2.expression_matrix.txt", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(expression_part3, file = paste(out_path, "/part3.expression_matrix.txt", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(covariates_part1, file = paste(out_path, "/part1.expression_covariates.txt", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(covariates_part2, file = paste(out_path, "/part2.expression_covariates.txt", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(covariates_part3, file = paste(out_path, "/part3.expression_covariates.txt", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE)

