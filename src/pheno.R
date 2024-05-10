library(data.table)
library(tidyverse)

# check pheno file & make it be in the correct format

# read in pheno
# pheno <- fread("./results/patient_information/pheno/pheno_first_tx.txt", data.table=F, header=T) # does not have aGvHD 2-4 vs rest
pheno <- fread("./results/patient_information/pheno/pheno_first_tx_original_coding.txt", data.table=F, header=T) 

colnames(pheno)
# [1] "ID"                   
# [2] "GenotypingID"         
# [3] "DonorGenotypingID"    
# [4] "DonorType"            
# [5] "aGvHD_status_original"
# [6] "cGvHD_status_original"
# [7] "Relapse_original" 

# add levels needed

# aGvHD 2-4 vs 0-1
unique(pheno$aGvHD_status_original)
# "7"   NA    "2"   "4"   "3"   "1"   "Yes"
pheno_aGvHD <- rep(NA, length(pheno$ID))
pheno_aGvHD[pheno$aGvHD_status_original == 7] <- 0
pheno_aGvHD[pheno$aGvHD_status_original == 1] <- 0
pheno_aGvHD[pheno$aGvHD_status_original == 2] <- 1
pheno_aGvHD[pheno$aGvHD_status_original == 3] <- 1
pheno_aGvHD[pheno$aGvHD_status_original == 4] <- 1
pheno_aGvHD[pheno$aGvHD_status_original == "Yes"] <- NA # only 3 finns
pheno$pheno_aGvHD_2_4 <- pheno_aGvHD

# cGvHD yes,limited,extensive vs no
unique(pheno$cGvHD_status_original)
# "7"   NA    "2"   "1"   "Yes"
pheno_cGvHD <- rep(NA, length(pheno$ID))
pheno_cGvHD[pheno$cGvHD_status_original == 7] <- 0
pheno_cGvHD[pheno$cGvHD_status_original == 1] <- 1
pheno_cGvHD[pheno$cGvHD_status_original == 2] <- 1
pheno_cGvHD[pheno$cGvHD_status_original == "Yes"] <- 1 
pheno$pheno_cGvHD <- pheno_cGvHD


# relapse
unique(pheno$Relapse_original)
#  0  1 NA
pheno_relapse <- rep(NA, length(pheno$ID))
pheno_relapse[pheno$Relapse_original == 0] <- 0
pheno_relapse[pheno$Relapse_original == 1] <- 1
pheno$pheno_relapse <- pheno_relapse

# leave in only donor IDs
pheno_donors <- pheno[,c(3,3,8:10)]
colnames(pheno_donors)[1] <- "#FID"
colnames(pheno_donors)[2] <- "IID"
# remove rows where genotyping ID is NA
pheno_donors <- pheno_donors[!(is.na(pheno_donors$IID)),]

write.table(pheno_donors, "./results/thymic_SNP/pheno_plink.txt", sep = "\t", quote = F, row.names = F, col.names = T)

############################################################################

# association test run and getting errors for it - need to modiphy cGvHD and relapse phenos to be separare files 
# modify pheno for poland: remove aGvHD

pheno <- fread("./results/thymic_SNP/pheno_plink.txt")
pheno <- pheno[,-"pheno_aGvHD_2_4"]
write.table(pheno, "./results/thymic_SNP/pheno_plink_mod.txt", sep = "\t", quote = F, row.names = F, col.names = T)

pheno <- pheno[,-"pheno_cGvHD"]
write.table(pheno, "./results/thymic_SNP/pheno_plink_mod_2.txt", sep = "\t", quote = F, row.names = F, col.names = T)











