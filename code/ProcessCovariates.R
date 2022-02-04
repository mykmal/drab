# Set tissue names
tissues <- c("Adipose_Subcutaneous", "Adipose_Visceral_Omentum", "Adrenal_Gland", "Artery_Aorta", "Artery_Coronary", "Artery_Tibial", "Brain_Amygdala", "Brain_Anterior_cingulate_cortex_BA24", "Brain_Caudate_basal_ganglia", "Brain_Cerebellar_Hemisphere", "Brain_Cerebellum", "Brain_Cortex", "Brain_Frontal_Cortex_BA9", "Brain_Hippocampus", "Brain_Hypothalamus", "Brain_Nucleus_accumbens_basal_ganglia", "Brain_Putamen_basal_ganglia", "Brain_Spinal_cord_cervical_c-1", "Brain_Substantia_nigra", "Breast_Mammary_Tissue", "Cells_Cultured_fibroblasts", "Cells_EBV-transformed_lymphocytes", "Colon_Sigmoid", "Colon_Transverse", "Esophagus_Gastroesophageal_Junction", "Esophagus_Mucosa", "Esophagus_Muscularis", "Heart_Atrial_Appendage", "Heart_Left_Ventricle", "Kidney_Cortex", "Liver", "Lung", "Minor_Salivary_Gland", "Muscle_Skeletal", "Nerve_Tibial", "Ovary", "Pancreas", "Pituitary", "Prostate", "Skin_Not_Sun_Exposed_Suprapubic", "Skin_Sun_Exposed_Lower_leg", "Small_Intestine_Terminal_Ileum", "Spleen", "Stomach", "Testis", "Thyroid", "Uterus", "Vagina", "Whole_Blood")

for (tissue in tissues) {
        # Load data
        covar_path <- paste("expression_covariates/", tissue, ".v8.EUR.covariates.txt", sep = "")
        covariates <- read.table(file = covar_path, row.names = 1)

        # Rename the ID row
        rownames(covariates)[1] <- "IID"
        
        # Transpose the covariate matrix
        covariates <- as.data.frame(t(covariates))
        
        # Add column of zeros for the family IDs
        FID <- matrix(0L, nrow = length(covariates$IID), ncol = 1)
        covariates <- cbind(FID, covariates)

        # Save to a new file
        processed_covar_path <- paste("expression_covariates/", tissue, ".v8.EUR.covariates.plink.txt", sep = "")
        write.table(covariates, file = processed_covar_path, quote = FALSE, sep = "\t", row.names = FALSE)
}
