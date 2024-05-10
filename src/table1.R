library(tidyverse)
library(data.table)
library(flextable)

#####################################################################################################################################

# table 1 for the thymus SNP

#####################################################################################################################################

# clinical and laboratory information
# read in the processed clinical information
harmonized_processed <- fread("/home/nithiju/work/HSCT_predictor/results/patient_information/harmonized_all_first_tx_no_duplicates.txt")

# add hla-matching status
matching <- fread("./results/HLA_imputation/imputed_HLA/HLA_matching_status.txt", data.table=F, header=T)
# "./results/thymic_SNP/IDs_donors_10_10.txt"
# "./results/thymic_SNP/IDs_donors_9_10.txt"
matching <- matching[match(harmonized_processed$DonorGenotypingID, matching$donor),]
harmonized_processed$HLA_match <- matching$matches


# check that all have at least one phenotype tested 
# pheno
pheno <- fread("./results/thymic_SNP/pheno_plink.txt", data.table=F, header=T)
has_pheno <- !(is.na(pheno$pheno_aGvHD_2_4)) | !(is.na(pheno$pheno_cGvHD)) | !(is.na(pheno$pheno_relapse))
pheno <- pheno[has_pheno,]
# check that all have a genotyping ID in the genotype files
geno <- fread("./results/thymic_SNP/extracted_snp/all_datasets_SNP_ref_A.fam", data.table=F, header=F)
has_pheno_geno <- pheno[pheno$IID %in% geno$V2,] # table with IDs (only donors)



# survival
survival <- fread("./results/thymic_SNP/survival/survival_data_thymus.txt")
survival <- survival[!(is.na(survival$chr14_22448641_G_A_A)),] # in survival analysis, division into groups based on the donor genotype -> also leaves in "only donors"
has_OS <- !(is.na(survival$overall_survival) & is.na(survival$OS_status))
has_RFS <- !(is.na(survival$relapse_free_survival) & is.na(survival$RFS_status))
has_some_survival <- has_OS | has_RFS
all_with_survival <- survival[has_some_survival,] # table with IDs 



leave_pheno <- harmonized_processed$DonorGenotypingID %in% has_pheno_geno$IID
leave_survival <- harmonized_processed$DonorGenotypingID %in% all_with_survival$DonorGenotypingID
leave <- leave_pheno | leave_survival

harmonized_processed <- harmonized_processed[leave,] # 1473


# ID list for filtering all donors
# all donors
donors <- fread("./results/patient_information/ID_lists/IDs_donors.txt", data.table=F, header=F) # 1519
# unrelated donors
unrelated <- fread("./results/thymic_SNP/IDs_donors_unrelated.txt", data.table=F, header=F) # 450

# use these to extract information below
harmonized_all <- harmonized_processed[harmonized_processed$DonorGenotypingID %in% donors$V2,]  # 1473
harmonized_unrel<- harmonized_processed[harmonized_processed$DonorGenotypingID %in% unrelated$V2,] # 427


#####################################################################################################################################

# gather all results into this table
table_all <- data.table(subset=c("All donors", "Unrelated donors"))

#-------------------------------------------------------------------------------------------------

# N

table_all$N <- c(nrow(harmonized_all), nrow(harmonized_unrel))

#-------------------------------------------------------------------------------------------------

# Population

table_all$Finland <- c(sum(harmonized_all$population %in% c("ic1", "ic3", "mcgill", "tyks", "hus")), sum(harmonized_unrel$population %in% c("ic1", "ic3", "mcgill", "tyks", "hus")))

table_all$UK <- c(sum(harmonized_all$population == "newcastle"), sum(harmonized_unrel$population == "newcastle"))

table_all$Spain <- c(sum(harmonized_all$population == "katalonia"), sum(harmonized_unrel$population == "katalonia"))

table_all$Poland <- c(sum(harmonized_all$population == "poland"), sum(harmonized_unrel$population == "poland"))

#-------------------------------------------------------------------------------------------------

# HLA-match

table_all$HLA_match_10 <- c(sum(harmonized_all$HLA_match == "10/10", na.rm = T), sum(harmonized_unrel$HLA_match == "10/10", na.rm = T))

table_all$HLA_match_9 <- c(sum(harmonized_all$HLA_match == "9/10", na.rm = T), sum(harmonized_unrel$HLA_match == "9/10", na.rm = T))

table_all$HLA_match_other <- c(sum(!(harmonized_all$HLA_match %in% c("10/10", "9/10"))), sum(!(harmonized_unrel$HLA_match %in% c("10/10", "9/10")))) # this fixed now

#-------------------------------------------------------------------------------------------------

# Donortype

table_all$Donortype_unrelated <- c(sum(harmonized_all$DonorGenotypingID %in% harmonized_unrel$DonorGenotypingID), nrow(harmonized_unrel))

table_all$Donortype_related <- c(sum(!(harmonized_all$DonorGenotypingID %in% harmonized_unrel$DonorGenotypingID)), 0)


#-------------------------------------------------------------------------------------------------

# tXYears

get_tx_years <- function(clin_table){
  
  return(paste0(min(clin_table$TxYear, na.rm = T), "–", max(clin_table$TxYear, na.rm = T)))
  
}


table_all$tx_years <- c(get_tx_years(harmonized_all),
                        get_tx_years(harmonized_unrel))

table_all$tx_years_NA <- c(sum(is.na(harmonized_all$TxYear)),
                           sum(is.na(harmonized_unrel$TxYear)))

#-------------------------------------------------------------------------------------------------

# donor and recipient ages

R_age <- function(clin_table){
  
  return(paste0(signif(median(clin_table$AgeAll, na.rm = T), digits = 2), " (", signif(min(clin_table$AgeAll, na.rm = T), digits = 2), "–", signif(max(clin_table$AgeAll, na.rm = T), digits = 2), ")"))
  
}



table_all$recipient_age_median_range <- c(R_age(harmonized_all),
                                          R_age(harmonized_unrel))

table_all$recipient_age_NA <- c(sum(is.na(harmonized_all$AgeAll)),
                                sum(is.na(harmonized_unrel$AgeAll)))



D_age <- function(clin_table){
  
  return(paste0(signif(median(clin_table$DonorAgeAll, na.rm = T), digits = 2), " (", signif(min(clin_table$DonorAgeAll, na.rm = T), digits = 2), "–", signif(max(clin_table$DonorAgeAll, na.rm = T), digits = 2), ")"))
  
}


table_all$donor_age_median_range <- c(D_age(harmonized_all),
                                      D_age(harmonized_unrel)) 


table_all$donor_age_NA <- c(sum(is.na(harmonized_all$DonorAgeAll)),
                            sum(is.na(harmonized_unrel$DonorAgeAll)))

#---------------------------------------------------------------------------------------------------------

# donor and recipient sex match

sex_match <- function(table){
  
  table$sex_match <- paste0(table$GenderAll, table$DonorGenderAll)
  # coded as M or F 
  # many have NA as well -> remove them
  table$sex_match[table$sex_match %like% "NA"] <- NA
  
  # what is returned in which order
  # D_M+R_M   D_M+R_F   D_F+R_M   D_F+R_F   one or both sexes are NA
  
  return(c(sum(table$sex_match == "MM", na.rm = T), sum(table$sex_match == "FM", na.rm = T), sum(table$sex_match == "MF", na.rm = T), sum(table$sex_match == "FF", na.rm = T), sum(is.na(table$sex_match))))
  
  
}

# get all sex matches in each population into one table
sex_matches <- data.table()
sex_matches$all <- sex_match(harmonized_all)
sex_matches$unrel <- sex_match(harmonized_unrel)
sex_matches <- t(sex_matches)
colnames(sex_matches) <- c("donor M, recipient M", "donor M, recipient F", "donor F, recipient M", "donor F, recipient F", "either or both sexes NA")

table_all <- cbind(table_all, sex_matches)

#---------------------------------------------------------------------------------------------------------

# graft

get_grafts <- function(table){
  
  # what is returned in which order
  # PB, BM, PB+BM, NA
  
  return(c(sum(table$Graft == "PB", na.rm = T), sum(table$Graft == "BM", na.rm = T), sum(table$Graft == "PB+BM", na.rm = T), sum(is.na(table$Graft))))
  
  
}

# get all in each population into one table
graft <- data.table()
graft$fin <- get_grafts(harmonized_all)
graft$spain <- get_grafts(harmonized_unrel)
graft <- t(graft)
colnames(graft) <- c("graft PB", "graft BM", "graft PB+BM", "graft NA")

table_all <- cbind(table_all, graft)

#---------------------------------------------------------------------------------------------------------

# preconditioning

unique(harmonized_all$PreconditioningType)
# [1] "MAC"  "RIC"  "sekv" NA 
# does not contain other types in all data (NMA)

get_preconditioning <- function(table){
  
  # what is returned in which order
  # MAC, RIC, sekv, NA
  
  return(c(sum(table$PreconditioningType == "MAC", na.rm = T), sum(table$PreconditioningType %in% c("RIC", "sekv", "NMA"), na.rm = T), sum(is.na(table$PreconditioningType))))
  
  
}

# get all in each population into one table
preconditioning <- data.table()
preconditioning$fin <- get_preconditioning(harmonized_all)
preconditioning$spain <- get_preconditioning(harmonized_unrel)
preconditioning <- t(preconditioning)
colnames(preconditioning) <- c("preconditioning MAC", "preconditioning RIC", "preconditioning NA")

table_all <- cbind(table_all, preconditioning)

#---------------------------------------------------------------------------------------------------------

# endpoints: aGvHD
unique(harmonized_all$aGvHD_status)
# [1] "7" NA  "2" "4" "3" "1"

get_aGvHD <- function(table){
  
  return(c(sum(table$aGvHD_status %in% c("2", "3", "4"), na.rm = T), sum(is.na(table$aGvHD_status))))
  
  
}

# get all in each population into one table
aGvHD <- data.table()
aGvHD$fin <- get_aGvHD(harmonized_all)
aGvHD$spain <- get_aGvHD(harmonized_unrel)
aGvHD <- t(aGvHD)
colnames(aGvHD) <- c("aGvHD grade 2-4", "aGvHD NA")

table_all <- cbind(table_all, aGvHD)

#---------------------------------------------------------------------------------------------------------

# endpoints: cGvHD
unique(harmonized_all$cGvHD_status)
# [1] "7"   NA    "2"   "1"   "Yes"

get_cGvHD <- function(table){
  
  return(c(sum(table$cGvHD_status %in% c("1", "2", "Yes"), na.rm = T), sum(is.na(table$cGvHD_status))))
  
  
}

# get all in each population into one table
cGvHD <- data.table()
cGvHD$fin <- get_cGvHD(harmonized_all)
cGvHD$spain <- get_cGvHD(harmonized_unrel)
cGvHD <- t(cGvHD)
colnames(cGvHD) <- c("cGvHD grade 1-2", "cGvHD NA")

table_all <- cbind(table_all, cGvHD)

#---------------------------------------------------------------------------------------------------------

# endpoints: relapse

get_relapse <- function(table){
  
  return(c(sum(table$Relapse == "Y", na.rm = T), sum(table$Relapse == "N", na.rm = T), sum(is.na(table$Relapse))))
  
  
}

# get all in each population into one table
relapse <- data.table()
relapse$fin <- get_relapse(harmonized_all)
relapse$spain <- get_relapse(harmonized_unrel)
relapse <- t(relapse)
colnames(relapse) <- c("relapse yes", "relapse no", "relapse NA")

table_all <- cbind(table_all, relapse)

#---------------------------------------------------------------------------------------------------------

# endpoints: survival

survival <- fread("./results/thymic_SNP/survival/survival_data_thymus.txt")
survival_all <- survival[match(harmonized_all$DonorGenotypingID, survival$DonorGenotypingID),] 
survival_unrel <- survival[match(harmonized_unrel$DonorGenotypingID, survival$DonorGenotypingID),] 



all_OS <- !(is.na(survival_all$overall_survival) & is.na(survival_all$OS_status))
unrel_OS <- !(is.na(survival_unrel$overall_survival) & is.na(survival_unrel$OS_status))
all_RFS <- !(is.na(survival_all$relapse_free_survival) & is.na(survival_all$RFS_status))
unrel_RFS <- !(is.na(survival_unrel$relapse_free_survival) & is.na(survival_unrel$RFS_status))

table_all$OS <- c(sum(all_OS), sum(unrel_OS))
table_all$OS_missing <- c(sum(!all_OS), sum(!unrel_OS))

table_all$RFS <- c(sum(all_RFS), sum(unrel_RFS))
table_all$RFS_missing <- c(sum(!all_RFS), sum(!unrel_RFS))

##################################################################################################################################

# transpose table to have populations as columns

table_all_T <- t(table_all)

write.table(table_all_T, file="./results/thymic_SNP/table_1.txt", sep="\t", quote=F, row.names=T, col.names=T)

table_all_T <- fread("./results/thymic_SNP/table_1.txt", header = T)

first <- c("Number of HSCT donors, n", "Country of origin, n (%)", "Country of origin, n (%)", "Country of origin, n (%)", "Country of origin, n (%)", "HLA-match level, n (%)", "HLA-match level, n (%)", "HLA-match level, n (%)", "Donor type, n (%)", "Donor type, n (%)", "HSCT time, years", "HSCT time, years", "Recipient age in years, median (range)", "Recipient age in years, median (range)", "Donor age in years, median (range)", "Donor age in years, median (range)","Donor-recipient gender, n (%)", "Donor-recipient gender, n (%)", "Donor-recipient gender, n (%)", "Donor-recipient gender, n (%)", "Donor-recipient gender, n (%)", "Stem cell source, n (%)", "Stem cell source, n (%)", "Stem cell source, n (%)", "Stem cell source, n (%)", "Conditioning regimen, n (%)", "Conditioning regimen, n (%)", "Conditioning regimen, n (%)", "aGvHD, n (%)", "aGvHD, n (%)", "cGvHD, n (%)", "cGvHD, n (%)", "Relapse, n (%)", "Relapse, n (%)", "Relapse, n (%)", "Overall survival, n (%)", "Overall survival, n (%)", "Relapse-free survival, n (%)", "")
second <- c("", "Finland", "UK", "Spain", "Poland", "10/10", "9/10", "Other", "Unrelated", "Related", "", "Missing, n (%)", "", "Missing, n (%)", "", "Missing, n (%)", "Male-male", "Male-female", "Female-male", "Female-female", "Missing", "Peripheral blood", "Bone marrow", "Both", "Missing", "Myeloablative", "Reduced intensity", "Missing", paste0("grade ", as.roman(2), "-", as.roman(4)), "Missing", "Limited or extensive", "Missing", "Yes", "No", "Missing", "Found", "Missing", "Found", "Missing")


table <- cbind(NA, NA, table_all_T[,2:3])
table <- as.data.frame(table)
table[,1] <- first
table[,2] <- second

colnames(table) <- c(" ", "  ", "All donors", "Unrelated donors")

# add %
rows <- 1:nrow(table)
rows <- rows[-c(1,11,13,15)] 

for (j in 3:4) {
  
  # i = row
  # j = col
  
  for (i in rows) {
    
    value <- (as.numeric(table[i,j]) / as.numeric(table[1,j])) * 100
    
    pasted <- paste0(table[i,j], " (", round(value), ")")
    
    table[i,j] <- pasted
    
  }
  
}



# headers
table_flex <- flextable(table)
table_flex <- autofit(table_flex, add_w = 0, add_h = 0)
# table_flex <- align(table_flex, part = "header", align = "left")

# merge vertical duplicated names
table_flex <- merge_v(table_flex, j = c(" ", "  ")) # which cols
table_flex <- valign(table_flex, j = 1, valign = "top")

# add footnote
table_flex <- add_footer_lines(table_flex, "GvHD, graft-versus-host disease; aGvHD, acute GvHD; cGvHD, chronic GvHD")
table_flex <- footnote(table_flex, i = c(13,15), j = 1,
                       value = as_paragraph(
                         c("Missing ages were imputed")
                       ),
                       ref_symbols = c("*"),
                       part = "body")


set_flextable_defaults(background.color = "white")

save_as_image(table_flex, path = "./results/thymic_SNP/table_1_flextable.png", bg = "white")
save_as_docx(table_flex, path = "./results/thymic_SNP/table_1_flextable.docx", align = "left")

