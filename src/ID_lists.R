library(data.table)
library(tidyverse)

# ID lists to keep for the thymic SNP analysis

# donors
donors <- fread("./results/patient_information/ID_lists/IDs_donors.txt", data.table=F, header=F)
# already ok

# donors unrelated
donors <- fread("./results/patient_information/ID_lists/IDs_donors_register_tx.txt", data.table=F, header=F) 
write.table(donors[,1:2], "./results/thymic_SNP/IDs_donors_unrelated.txt", quote = F, col.names = F, row.names = F)

# get hla matching 10/10 or 9/10 or other
matching <- fread("./results/HLA_imputation/imputed_HLA/HLA_matching_status.txt", data.table=F, header=T)
# a few have NAs -> remove them
matching <- matching[!(is.na(matching$matches)),]
# a list of 10/10  and 9/10 matched
matched_10 <- matching[matching$matches == "10/10",]
matched_9 <- matching[matching$matches == "9/10",]

# combine these with register tx ID lists
matched_10_unrelated <- matched_10[matched_10$donor %in% donors$V2,]
matched_9_unrelated <- matched_9[matched_9$donor %in% donors$V2,]


write.table(matched_10[,c(2,2)], "./results/thymic_SNP/IDs_donors_10_10.txt", quote = F, col.names = F, row.names = F)
write.table(matched_9[,c(2,2)], "./results/thymic_SNP/IDs_donors_9_10.txt", quote = F, col.names = F, row.names = F)

write.table(matched_10_unrelated[,c(2,2)], "./results/thymic_SNP/IDs_donors_10_10_unrelated.txt", quote = F, col.names = F, row.names = F)
write.table(matched_9_unrelated[,c(2,2)], "./results/thymic_SNP/IDs_donors_9_10_unrelated.txt", quote = F, col.names = F, row.names = F)

