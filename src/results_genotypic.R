library(data.table)
library(tidyverse)

# make results tables out of the association results of the thymic SNP

###############################################################################################################################

# get folders

# one table for each subset

all_donors <- list.files(path = "./results/thymic_SNP/association_results_genotypic/plink1/allDonors", pattern = "glm.logistic.hybrid", full.names = T)
all_donors_10 <- list.files(path = "./results/thymic_SNP/association_results_genotypic/plink1/allDonors_10_10", pattern = "glm.logistic.hybrid", full.names = T)
all_donors_9 <- list.files(path = "./results/thymic_SNP/association_results_genotypic/plink1/allDonors_9_10", pattern = "glm.logistic.hybrid", full.names = T)

unrelated_donors <- list.files(path = "./results/thymic_SNP/association_results_genotypic/plink1/unrelatedDonors", pattern = "glm.logistic.hybrid", full.names = T)
unrelated_donors_10 <- list.files(path = "./results/thymic_SNP/association_results_genotypic/plink1/unrelatedDonors_10_10", pattern = "glm.logistic.hybrid", full.names = T)
unrelated_donors_9 <- list.files(path = "./results/thymic_SNP/association_results_genotypic/plink1/unrelatedDonors_9_10", pattern = "glm.logistic.hybrid", full.names = T)


make_results_table <- function(file_list){
  
  results <- data.frame(population=character(), pheno=character(), tested=character(), OR=numeric(), CI_low=numeric(), CI_up=numeric(), p=numeric())
  
  for (i in 1:length(file_list)) {
    
    file <- fread(file_list[i])
    
    filename <- str_split(file_list[i], "/")[[1]][7]
    
    # print(filename)
    
    pop <- str_split(filename, "_")[[1]][1]
    tested <- str_sub(str_split(filename, "_")[[1]][2], 1,1)
    pheno <- str_sub(str_split(filename, "\\.")[[1]][2], 7)
    
    results[nrow(results) + 1,] <- c(pop, pheno, tested, file[1,"OR"], file[1,"L95"], file[1,"U95"], file[1,"P"])
    
  }
  
  # modify test names
  results$tested[results$tested == "A"] <- "AA"
  results$tested[results$tested == "G"] <- "GG"
  results$tested[results$tested == "h"] <- "AG"
  
  # order by pop, then pheno, then genotype
  results <- results[order(results[,2], results[,1], results[,3]),]
  
  return(results)
  
}


res_all_donors <- make_results_table(all_donors)
res_all_donors_10 <- make_results_table(all_donors_10)
res_all_donors_9 <- make_results_table(all_donors_9)

res_unrelated_donors <- make_results_table(unrelated_donors)
res_unrelated_donors_10 <- make_results_table(unrelated_donors_10)
res_unrelated_donors_9 <- make_results_table(unrelated_donors_9)



# save results
write.table(res_all_donors, "./results/thymic_SNP/association_results_genotypic/plink1/res_all_donors.txt", quote = F, col.names = T, row.names = F, sep = "\t")
write.table(res_all_donors_10, "./results/thymic_SNP/association_results_genotypic/plink1/res_all_donors_10.txt", quote = F, col.names = T, row.names = F, sep = "\t")
write.table(res_all_donors_9, "./results/thymic_SNP/association_results_genotypic/plink1/res_all_donors_9.txt", quote = F, col.names = T, row.names = F, sep = "\t")

write.table(res_unrelated_donors, "./results/thymic_SNP/association_results_genotypic/plink1/res_unrelated_donors.txt", quote = F, col.names = T, row.names = F, sep = "\t")
write.table(res_unrelated_donors_10, "./results/thymic_SNP/association_results_genotypic/plink1/res_unrelated_donors_10.txt", quote = F, col.names = T, row.names = F, sep = "\t")
write.table(res_unrelated_donors_9, "./results/thymic_SNP/association_results_genotypic/plink1/res_unrelated_donors_9.txt", quote = F, col.names = T, row.names = F, sep = "\t")





