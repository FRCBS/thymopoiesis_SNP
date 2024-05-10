library(data.table)
library(tidyverse)

# make results plots out of the association results of the thymic SNP

###############################################################################################################################

# read in results tables

res_all_donors <- fread("./results/thymic_SNP/association_results_genotypic/plink1/res_all_donors.txt")
res_all_donors_10 <- fread("./results/thymic_SNP/association_results_genotypic/plink1/res_all_donors_10.txt")
res_all_donors_9 <- fread("./results/thymic_SNP/association_results_genotypic/plink1/res_all_donors_9.txt")

res_unrelated_donors <- fread("./results/thymic_SNP/association_results_genotypic/plink1/res_unrelated_donors.txt")
res_unrelated_donors_10 <- fread("./results/thymic_SNP/association_results_genotypic/plink1/res_unrelated_donors_10.txt")
res_unrelated_donors_9 <- fread("./results/thymic_SNP/association_results_genotypic/plink1/res_unrelated_donors_9.txt")

#-------------------------------------------------

# modify data to be in the correct format for plotting
# some columns in res_unrelated_donors_9 are characters in stead of numeric for some reason & this causes errors in plotting
# -> change them to be numeric separately here

res_unrelated_donors_9$CI_low <- as.numeric(res_unrelated_donors_9$CI_low)
res_unrelated_donors_9$CI_up <- as.numeric(res_unrelated_donors_9$CI_up)

res_unrelated_donors$CI_low <- as.numeric(res_unrelated_donors$CI_low)
res_unrelated_donors$CI_up <- as.numeric(res_unrelated_donors$CI_up)

#-------------------------------------------------

# change population names to be better

change_names <- function(results_table){
  
  print(nrow(results_table))
  print(unique(results_table$population))
  
  results_table$population[results_table$population == "all"] <- "Combination"
  results_table$population[results_table$population == "fin"] <- "Finland"
  results_table$population[results_table$population == "kat"] <- "Spain"
  results_table$population[results_table$population == "nc"] <- "UK"
  results_table$population[results_table$population == "pol"] <- "Poland"
  
  return(results_table)
  
}

res_all_donors <- change_names(res_all_donors)
res_all_donors_10 <- change_names(res_all_donors_10)
res_all_donors_9 <- change_names(res_all_donors_9)

res_unrelated_donors <- change_names(res_unrelated_donors)
res_unrelated_donors_10 <- change_names(res_unrelated_donors_10)
res_unrelated_donors_9 <- change_names(res_unrelated_donors_9)

##########################################33

# plots

# common theme elements
th <- theme(panel.spacing = unit(.5, "lines"), 
            panel.border = element_rect(color = "black", fill = NA, linewidth = 1), 
            strip.background = element_rect(color = "black", linewidth = 1), 
            text = element_text(size = 15), 
            legend.position = "bottom", 
            legend.text = element_text(size = 14),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            strip.text = element_text(size=15), 
            axis.title.x = element_text(size=15), 
            axis.text.x = element_text(size=12))

plot_results <- function(results_table, save_to){
  
  plot <- ggplot(results_table, aes(OR, tested))  +
    geom_point(size = 4) +
    geom_linerange(aes(xmin = CI_low, xmax = CI_up), linewidth = 1) + 
    geom_vline(xintercept = 1, color = "grey30", linetype = "dashed") +
    ylab("") + xlab("OR") +
    facet_grid(pheno ~ population, scales = 'free') +
    theme_minimal() + th
  
  ggsave(paste0(save_to, ".png"), bg = "white", width = 14, height = 8, dpi = 600)
  
}

plot_results(res_all_donors, "./results/thymic_SNP/association_results_genotypic/plink1/res_all_donors")
plot_results(res_all_donors_10, "./results/thymic_SNP/association_results_genotypic/plink1/res_all_donors_10")
plot_results(res_all_donors_9, "./results/thymic_SNP/association_results_genotypic/plink1/res_all_donors_9")

# plot_results(res_unrelated_donors, "./results/thymic_SNP/association_results_genotypic/plink1/res_unrelated_donors") # error, one massive CI upper limit causing it
plot_results(res_unrelated_donors_10, "./results/thymic_SNP/association_results_genotypic/plink1/res_unrelated_donors_10")
# plot_results(res_unrelated_donors_9, "./results/thymic_SNP/association_results_genotypic/plink1/res_unrelated_donors_9") # error, one massive CI upper limit causing it

#---------------------------------------------------------------------------------------

# plot not working well for res_unrelated_donors due to massive CI_uppers
# add big ones as dotted lines or open circles

big <- res_unrelated_donors$CI_up > 50
res_unrelated_donors$CI_up[res_unrelated_donors$CI_up > 50] <- NA
dotted <- res_unrelated_donors[big,]
dotted$CI_low <- 0
dotted$CI_up <- 25
dotted$OR <- 1 # having the OR column is causing weird/too big x axis

ors <- res_unrelated_donors$OR > 15
res_unrelated_donors$OR[res_unrelated_donors$OR > 15] <- NA
circles <- res_unrelated_donors[ors,]
circles$OR <- 23

plot <- ggplot(res_unrelated_donors, aes(OR, tested))  +
  geom_point(size = 4) +
  geom_linerange(aes(xmin = CI_low, xmax = CI_up), linewidth = 1) + 
  geom_point(data=circles, size = 4, shape=1) +
  geom_linerange(data = dotted, aes(xmin = CI_low, xmax = CI_up), linewidth = 1, linetype = "dotted") +
  geom_vline(xintercept = 1, color = "grey30", linetype = "dashed") +
  ylab("") + xlab("OR") +
  facet_grid(pheno ~ population, scales = 'free') +
  theme_minimal() + th

ggsave("./results/thymic_SNP/association_results_genotypic/plink1/res_unrelated_donors.png", bg = "white", width = 14, height = 8, dpi = 600)



#---------------------------------------------------------------------------------------

# plot not working well for res_unrelated_donors_9 due to massive CI_uppers
# add big ones as dotted lines or open circles

big <- res_unrelated_donors_9$CI_up > 50
res_unrelated_donors_9$CI_up[res_unrelated_donors_9$CI_up > 50] <- NA
dotted <- res_unrelated_donors_9[big,]
dotted$CI_low <- 0
dotted$CI_up <- 25
dotted$OR <- 1 # having the OR column is causing weird/too big x axis

ors <- res_unrelated_donors_9$OR > 15
res_unrelated_donors_9$OR[res_unrelated_donors_9$OR > 15] <- NA
circles <- res_unrelated_donors_9[ors,]
circles$OR <- 23

plot <- ggplot(res_unrelated_donors_9, aes(OR, tested))  +
  geom_point(size = 4) +
  geom_linerange(aes(xmin = CI_low, xmax = CI_up), linewidth = 1) + 
  geom_point(data=circles, size = 4, shape=1) +
  geom_linerange(data = dotted, aes(xmin = CI_low, xmax = CI_up), linewidth = 1, linetype = "dotted") +
  geom_vline(xintercept = 1, color = "grey30", linetype = "dashed") +
  ylab("") + xlab("OR") +
  facet_grid(pheno ~ population, scales = 'free') +
  theme_minimal() + th

ggsave("./results/thymic_SNP/association_results_genotypic/plink1/res_unrelated_donors_9.png", bg = "white", width = 14, height = 8, dpi = 600)

