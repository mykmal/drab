data <- read.delim(file = "raw/phs000424.v8.pht002742.v8.p2.c1.GTEx_Subject_Phenotypes.GRU.txt", skip = 10, stringsAsFactors = FALSE)

phenotypes <- subset(data, RACE == 3, select = c(SUBJID, SEX, AGE))

FID <- matrix(data = 0L, nrow = nrow(phenotypes), ncol = 1)

plink_phenotypes <- cbind(FID, phenotypes)
colnames(plink_phenotypes) <- c("FID", "IID", "SEX", "AGE")

write.table(plink_phenotypes, file = "phenotypes.txt", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

