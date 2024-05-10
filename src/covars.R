library(data.table)
library(tidyverse)

# covariates for the association analysis

# reusing the covars from the NK receptor project here
covars <- fread("/home/nithiju/work/NK_killing_activity/results/pheno_and_covars/covars_donors_aGvHD_relapse_all_binary.txt")

colnames(covars)
# [1] "FID"                  "IID"                 
# [3] "donor_sibling"        "donor_haplo"         
# [5] "graft_PB"             "graft_PB_BM"         
# [7] "graft_NA"             "sex_combo_risk"      
# [9] "sex_combo_NA"         "preconditioning_MAC" 
# [11] "preconditioning_sekv" "preconditioning_NMA" 
# [13] "population_fin"       "population_spain"    
# [15] "population_poland"    "AgeAll"              
# [17] "DonorAgeAll"          "AML_not"             
# [19] "txyear"            

########################################################################################################

# modify covariates
# they are causing issues in the association testing

# all populations: donortype_haplo & population_poland correlate too much
# -> merge donortype_haplo to be in the donortype_register category
# already: sibling 0/1, haplo 0/1, register = sibling 0, haplo 0
# -> just remove the donortype_haplo column -> register & haplo both 0 in the donortype_sibling column

covars <- covars[,-4]
write.table(covars, "./results/thymic_SNP/covars_mod.txt", quote = F, col.names = T, row.names = F)

# add HLA-matching status for when running the analysis in all unrelated donors
matching <- fread("./results/HLA_imputation/imputed_HLA/HLA_matching_status.txt", data.table=F, header=T)
matching <- matching[match(covars$IID, matching$donor),]
# make the match be a number (now X/10)
split <- str_split(matching$matches, "/")
for (i in 1:nrow(matching)) {
  matching$matches[i] <- split[[i]][1]
}

covars_match <- cbind(covars, matching$matches)
colnames(covars_match)[ncol(covars_match)] <- "match_score"
write.table(covars_match, "./results/thymic_SNP/covars_mod_match.txt", quote = F, col.names = T, row.names = F)

########################################################################################################

# modify covariates for poland
# they are causing issues in the association testing

# 'graft_PB' and 'graft_PB_BM' correlate too much -> are the opposite of each other -> leave only one of the two

covars <- fread("./results/thymic_SNP/covars_mod.txt")
pol <- fread("./results/thymic_SNP/extracted_snp/poland_SNP_ref_A.fam")

covars_pol <- covars[covars$IID %in% pol$V2,]
# all but one have PB, one has PB+BM
covars <- covars[,-"graft_PB_BM"]
write.table(covars, "./results/thymic_SNP/covars_mod_pol.txt", quote = F, col.names = T, row.names = F)



matching <- fread("./results/HLA_imputation/imputed_HLA/HLA_matching_status.txt", data.table=F, header=T)
matching <- matching[match(covars$IID, matching$donor),]
# make the match be a number (now X/10)
split <- str_split(matching$matches, "/")
for (i in 1:nrow(matching)) {
  matching$matches[i] <- split[[i]][1]
}

covars_match <- cbind(covars, matching$matches)
colnames(covars_match)[ncol(covars_match)] <- "match_score"

write.table(covars_match, "./results/thymic_SNP/covars_mod_pol_match.txt", quote = F, col.names = T, row.names = F)


