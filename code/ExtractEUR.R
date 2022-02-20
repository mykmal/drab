# Read subject phenotypes file
phenotypes <- read.csv(file = "reference/phs000424.v8.pht002742.v8.p2.c1.GTEx_Subject_Phenotypes.GRU.txt", skip = 10, sep = "\t", header = TRUE)

# Drop unneeded phenotype data
EUR <- subset(phenotypes, RACE == 3, select = c(SUBJID))

# Create column of family IDs
zeros <- matrix(data = 0L, nrow = length(EUR), ncol = 1)

# Write plink format samples list
write.table(cbind(zeros, EUR), file = "reference/EUR_samples.txt", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)

