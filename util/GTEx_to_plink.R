tissues <-
  c(
    "Adipose_Subcutaneous",
    "Adipose_Visceral_Omentum",
    "Adrenal_Gland",
    "Artery_Aorta",
    "Artery_Coronary",
    "Artery_Tibial",
    "Brain_Amygdala",
    "Brain_Anterior_cingulate_cortex_BA24",
    "Brain_Caudate_basal_ganglia",
    "Brain_Cerebellar_Hemisphere",
    "Brain_Cerebellum",
    "Brain_Cortex",
    "Brain_Frontal_Cortex_BA9",
    "Brain_Hippocampus",
    "Brain_Hypothalamus",
    "Brain_Nucleus_accumbens_basal_ganglia",
    "Brain_Putamen_basal_ganglia",
    "Brain_Spinal_cord_cervical_c-1",
    "Brain_Substantia_nigra",
    "Breast_Mammary_Tissue",
    "Cells_Cultured_fibroblasts",
    "Cells_EBV-transformed_lymphocytes",
    "Colon_Sigmoid",
    "Colon_Transverse",
    "Esophagus_Gastroesophageal_Junction",
    "Esophagus_Mucosa",
    "Esophagus_Muscularis",
    "Heart_Atrial_Appendage",
    "Heart_Left_Ventricle",
    "Kidney_Cortex",
    "Liver",
    "Lung",
    "Minor_Salivary_Gland",
    "Muscle_Skeletal",
    "Nerve_Tibial",
    "Ovary",
    "Pancreas",
    "Pituitary",
    "Prostate",
    "Skin_Not_Sun_Exposed_Suprapubic",
    "Skin_Sun_Exposed_Lower_leg",
    "Small_Intestine_Terminal_Ileum",
    "Spleen",
    "Stomach",
    "Testis",
    "Thyroid",
    "Uterus",
    "Vagina",
    "Whole_Blood"
  )

EUR_individuals <- read.table(file = "phenotypes.txt")$IID

for (tissue in tissues) {
  
  expression <- read.table(file = paste("raw/GTEx_Analysis_v8_eQTL_expression_matrices/", tissue, ".v8.normalized_expression.bed.gz", sep = ""),
                           comment.char = "", check.names = FALSE)
  expression <- subset(expression, select = c("gene_id", EUR_individuals))
  expression <- as.data.frame(t(expression))
  FID <- matrix(0L, nrow = nrow(expression), ncol = 1)
  expression <- cbind(FID, expression)
  expression[1, 1] <- "FID"
  expression[1, 2] <- "IID"
  write.table(expression, file = paste("expression/", tissue, ".expression.txt", sep = ""),
              quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  
  covariates <- read.table(file = paste("raw/GTEx_Analysis_v8_eQTL_covariates/", tissue, ".v8.covariates.txt", sep = ""))
  covariates <- subset(covariates, select = c("ID", EUR_individuals))
  covariates <- as.data.frame(t(covariates))
  FID <- matrix(0L, nrow = nrow(covariates), ncol = 1)
  covariates <- cbind(FID, covariates)
  covariates[1, 1] <- "FID"
  covariates[1, 2] <- "IID"
  write.table(covariates, file = paste("covariates/", tissue, ".covariates.txt", sep = ""),
              quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}

