library(data.table)
library(tidyverse)

# modify HSCT survival data for studying the thymus SNP

##################################################################################################################################

# read in data
data <- fread("./results/patient_information/survival/HSCT_survival_data.txt")

# add covariates

# covars <- fread("./results/thymic_SNP/covars_mod.txt") # modified for thymus SNP
covars <- fread("./results/thymic_SNP/covars_mod_match.txt") # includes HLA-match score
covars <- covars[match(data$DonorGenotypingID, covars$IID),]

# add thymus SNP dosage
dosage <- fread("./results/thymic_SNP/survival/dosage.raw")
# make dosage be two classes: GG/AG vs. AA
# now the dosage calculated for the allele A (colname chr14_22448641_G_A_A)
dosage$chr14_22448641_G_A_A[dosage$chr14_22448641_G_A_A == 1] <- "GG or GA"
dosage$chr14_22448641_G_A_A[dosage$chr14_22448641_G_A_A == 0] <- "GG or GA" # already ok
dosage$chr14_22448641_G_A_A[dosage$chr14_22448641_G_A_A == 2] <- "AA"
dosage$chr14_22448641_G_A_A <- as.factor(dosage$chr14_22448641_G_A_A)

dosage <- dosage[match(data$DonorGenotypingID, dosage$IID),]

data <- cbind(data, covars[,c(3:19)], dosage[,7])

write.table(data, "./results/thymic_SNP/survival/survival_data_thymus.txt", sep = "\t", col.names = T, row.names = F, quote = F)

##################################################################################################################################

# add covariates where all levels are in one column

# read in data
data <- fread("./results/patient_information/survival/HSCT_survival_data.txt")

# add covariates
covars <- fread("./results/patient_information/covars/covars_first_tx_all_imputed_ages.txt")
covars <- covars[match(data$DonorGenotypingID, covars$DonorGenotypingID),]
# remove columns aGvHD, diagnosis
covars <- covars[,c(4:9,11,13:14)]
# make classes be factors (this works nicer when saving results tables in the survival analysis)
# also rename the levels so that their names are automatically ok in the results tables
covars$population <- factor(covars$population)
levels(covars$population) <- c("Finland", "Spain", "UK", "Poland") # was [1] "finland"   "katalonia" "newcastle" "poland"  before
covars$DonorType <- factor(covars$DonorType)
levels(covars$DonorType) <- c("Haplotype", "Unrelated", "Sibling") # was [1] "haplo"    "register" "sibling"   before
covars$Graft <- factor(covars$Graft)
# levels ok for graft already
covars$sex_combo_binary <- factor(covars$sex_combo_binary)
levels(covars$sex_combo_binary) <- c("F donor, M recipient", "Other combination") # was [1] "0" "1"   before
covars$preconditioning <- factor(covars$preconditioning)
# levels ok for preconditioning already
covars$AML_not <- factor(covars$AML_not)
levels(covars$AML_not) <- c("Other", "AML") # was [1] "0" "1"   before

colnames(covars) <- c("Population", "Donor type", "Graft", "Recipient age", "Donor age", "D-R sex match", "Preconditioning", "Diagnosis", "Transplantation year")


# add thymus SNP dosage
dosage <- fread("./results/thymic_SNP/survival/dosage.raw")
# make dosage be two classes: GG/AG vs. AA
# now the dosage calculated for the allele A (colname chr14_22448641_G_A_A)
dosage$chr14_22448641_G_A_A[dosage$chr14_22448641_G_A_A == 1] <- "GG or GA"
dosage$chr14_22448641_G_A_A[dosage$chr14_22448641_G_A_A == 0] <- "GG or GA" # already ok
dosage$chr14_22448641_G_A_A[dosage$chr14_22448641_G_A_A == 2] <- "AA"
dosage$chr14_22448641_G_A_A <- as.factor(dosage$chr14_22448641_G_A_A)

dosage <- dosage[match(data$DonorGenotypingID, dosage$IID),]

# also add HLA-match score
covars_2 <- fread("./results/thymic_SNP/covars_mod_match.txt") # includes HLA-match score
covars_2 <- covars_2[match(data$DonorGenotypingID, covars_2$IID),]


data <- cbind(data, covars, covars_2[,19], dosage[,7])

write.table(data, "./results/thymic_SNP/survival/survival_data_thymus_factors.txt", sep = "\t", col.names = T, row.names = F, quote = F)



