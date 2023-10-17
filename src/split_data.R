args <- commandArgs(trailingOnly = TRUE)

context_A <- args[1]
context_B <- args[2]
out_path <- args[3]

if (context_A == context_B) {
  
  context <- context_A
  
  split_size <- 1 / 3
  
  expression <- read.table(paste("expression/", context, ".expression.txt", sep = ""), header = TRUE)
  covariates <- read.table(paste("covariates/", context, ".covariates.txt", sep = ""), header = TRUE)
  
  expression$ID <- paste(expression$FID, expression$IID, sep = "_")
  covariates$ID <- paste(covariates$FID, covariates$IID, sep = "_")
  individuals <- intersect(expression$ID, covariates$ID)
  
  individuals_part1 <- sample(individuals, size = ceiling(split_size * length(individuals)), replace = FALSE)
  individuals_part2 <- sample(setdiff(individuals, individuals_part1), size = length(individuals_part1), replace = FALSE)
  individuals_part3 <- setdiff(individuals, union(individuals_part1, individuals_part2))
  
  expression_part1 <- expression[expression$ID %in% individuals_part1, ]
  expression_part2 <- expression[expression$ID %in% individuals_part2, ]
  expression_part3 <- expression[expression$ID %in% individuals_part3, ]
  covariates_part1 <- covariates[covariates$ID %in% individuals_part1, ]
  covariates_part2 <- covariates[covariates$ID %in% individuals_part2, ]
  covariates_part3 <- covariates[covariates$ID %in% individuals_part3, ]
  
} else {
  
  expression_A <- read.table(paste("expression/", context_A, ".expression.txt", sep = ""), header = TRUE)
  covariates_A <- read.table(paste("covariates/", context_A, ".covariates.txt", sep = ""), header = TRUE)
  expression_B <- read.table(paste("expression/", context_B, ".expression.txt", sep = ""), header = TRUE)
  covariates_B <- read.table(paste("covariates/", context_B, ".covariates.txt", sep = ""), header = TRUE)
  
  expression_A$ID <- paste(expression_A$FID, expression_A$IID, sep = "_")
  covariates_A$ID <- paste(covariates_A$FID, covariates_A$IID, sep = "_")
  individuals_A <- intersect(expression_A$ID, covariates_A$ID)
  
  expression_B$ID <- paste(expression_B$FID, expression_B$IID, sep = "_")
  covariates_B$ID <- paste(covariates_B$FID, covariates_B$IID, sep = "_")
  individuals_B <- intersect(expression_B$ID, covariates_B$ID)
  
  if (length(individuals_A) < length(individuals_B)) {
    
    individuals_part1 <- individuals_A
    mutual <- intersect(individuals_B, individuals_part1)
    individuals_part2 <- union(mutual, sample(setdiff(individuals_B, individuals_part1), size = (length(individuals_part1) - length(mutual)), replace = FALSE))
    individuals_part3 <- setdiff(individuals_B, individuals_part2)
    
    expression_part3 <- expression_B[expression_B$ID %in% individuals_part3, ]
    covariates_part3 <- covariates_B[covariates_B$ID %in% individuals_part3, ]
    
  } else {
    
    individuals_part2 <- individuals_B
    mutual <- intersect(individuals_A, individuals_part2)
    individuals_part1 <- union(mutual, sample(setdiff(individuals_A, individuals_part2), size = (length(individuals_part2) - length(mutual)), replace = FALSE))
    individuals_part3 <- setdiff(individuals_A, individuals_part1)
    
    expression_part3 <- expression_A[expression_A$ID %in% individuals_part3, ]
    covariates_part3 <- covariates_A[covariates_A$ID %in% individuals_part3, ]
  }
  
  expression_part1 <- expression_A[expression_A$ID %in% individuals_part1, ]
  expression_part2 <- expression_B[expression_B$ID %in% individuals_part2, ]
  covariates_part1 <- covariates_A[covariates_A$ID %in% individuals_part1, ]
  covariates_part2 <- covariates_B[covariates_B$ID %in% individuals_part2, ]
}

expression_part1 <- subset(expression_part1, select = -ID)
expression_part2 <- subset(expression_part2, select = -ID)
expression_part3 <- subset(expression_part3, select = -ID)
covariates_part1 <- subset(covariates_part1, select = -ID)
covariates_part2 <- subset(covariates_part2, select = -ID)
covariates_part3 <- subset(covariates_part3, select = -ID)

write.table(expression_part1, file = paste(out_path, "/da.expression.txt", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(expression_part2, file = paste(out_path, "/db.expression.txt", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(expression_part3, file = paste(out_path, "/dt.expression.txt", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(covariates_part1, file = paste(out_path, "/da.covariates.txt", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(covariates_part2, file = paste(out_path, "/db.covariates.txt", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(covariates_part3, file = paste(out_path, "/dt.covariates.txt", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE)

