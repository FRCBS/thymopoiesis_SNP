library(data.table)
library(tidyverse)

# get sample numbers for the manuscript

###################################################################################################################################################

# samples all in all with the thymus SNP genotyped


# donors
donors <- fread("./results/patient_information/ID_lists/IDs_donors.txt", data.table=F, header=F)
# genotype file
geno <- fread("./results/thymic_SNP/extracted_snp/all_datasets_SNP.fam", data.table=F, header=F)
# pheno
pheno <- fread("./results/thymic_SNP/pheno_plink.txt", data.table=F, header=T)
# survival
survival <- fread("./results/thymic_SNP/survival/survival_data_thymus.txt")
survival <- survival[match(pheno$IID, survival$DonorGenotypingID),]

pheno_surv <- cbind(pheno, survival[,3:6])

###################################################################################################################################################

# how many after filtering with ID files

# association genotypic log file says
# allDonors 1491
test_all <- geno[geno$V2 %in% donors$V2,] # 1491
# fin 765
# kat 272
# nc 352
# pol 102


# how may who have at least one phenotype or survival available of these
pheno_surv_no_missing <- pheno_surv[!(is.na(pheno_surv$pheno_aGvHD_2_4)) | !(is.na(pheno_surv$pheno_cGvHD)) | !(is.na(pheno_surv$pheno_relapse)) | !(is.na(pheno_surv$overall_survival)) | !(is.na(pheno_surv$relapse_free_survival)),]


test_all_pheno <- test_all[test_all$V2 %in% pheno_surv_no_missing$IID,] # 1473


# get the numbers for each population

fin <- fread("./results/thymic_SNP/extracted_snp/finland_SNP.fam", data.table=F, header=F)
spain <- fread("./results/thymic_SNP/extracted_snp/spain_SNP.fam", data.table=F, header=F)
uk <- fread("./results/thymic_SNP/extracted_snp/uk_SNP.fam", data.table=F, header=F)
poland <- fread("./results/thymic_SNP/extracted_snp/poland_SNP.fam", data.table=F, header=F)

test_all_pheno_fin <- test_all_pheno[test_all_pheno$V2 %in% fin$V2,] # 765
test_all_pheno_spain <- test_all_pheno[test_all_pheno$V2 %in% spain$V2,] # 272
test_all_pheno_uk <- test_all_pheno[test_all_pheno$V2 %in% uk$V2,] # 334
test_all_pheno_poland <- test_all_pheno[test_all_pheno$V2 %in% poland$V2,] # 102

#----------------------------------

# how many (of the 1473) were 9/10 HLA-matched?

matched_9_10 <- fread("./results/thymic_SNP/IDs_donors_9_10.txt", data.table=F, header=F)
test_all_pheno_9_10 <- test_all_pheno[test_all_pheno$V2 %in% matched_9_10$V2,] # 765


