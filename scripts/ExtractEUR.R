# Read subject phenotypes file
phenotypes <- read.csv(file = "raw/phs000424.v8.pht002742.v8.p2.c1.GTEx_Subject_Phenotypes.GRU.txt", skip = 10, sep = "\t", header = TRUE)

# Drop unneeded phenotype data
europeans <- subset(phenotypes, RACE == 3, select = c(SUBJID))

# Create column of family IDs
zeros <- matrix(data = 0L, nrow = length(europeans), ncol = 1)

# Write plink format samples list
write.table(cbind(zeros, europeans), file = "EUR_samples.txt", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)

