library(data.table)
library(tidyverse)
library(survival)
library(ggsurvfit)
library(gtsummary)
library(flextable)
library(cowplot)

# survival analysis for the thymus SNP

###############################################################################################################################################

# univariate analysis

# read in data

data <- fread("./results/thymic_SNP/survival/survival_data_thymus.txt")
# rename the genotype levels for the plots
data$chr14_22448641_G_A_A[data$chr14_22448641_G_A_A == "AA"] <- "Donor AA"
data$chr14_22448641_G_A_A[data$chr14_22448641_G_A_A == "GG or GA"] <- "Donor AG/GG"
data$chr14_22448641_G_A_A <- factor(data$chr14_22448641_G_A_A, levels=c( "Donor AA", "Donor AG/GG"))

# read in files to filter the data into subgroups
files <- c("./results/patient_information/ID_lists/IDs_donors.txt", "./results/thymic_SNP/IDs_donors_10_10.txt", "./results/thymic_SNP/IDs_donors_9_10.txt", "./results/thymic_SNP/IDs_donors_unrelated.txt", "./results/thymic_SNP/IDs_donors_10_10_unrelated.txt", "./results/thymic_SNP/IDs_donors_9_10_unrelated.txt")
names <- c("all_donors", "all_10_10", "all_9_10", "unrelated", "unrelated_10_10", "unrelated_9_10")
titles <- c("all donors", "all 10/10 HLA-matched donors", "all 9/10 HLA-matched donors", "unrelated donors", "unrelated 10/10 HLA-matched donors", "unrelated  9/10 HLA-matched donors")


for (i in 1:length(files)) {
  
  donors <- fread(files[i], header = F)
  data_sub <- data[data$DonorGenotypingID %in% donors$V2,]
  
  print(names[i])
  
  # OS
  
  # modified plot
  OS_plot <- survfit2(Surv(overall_survival, OS_status) ~ chr14_22448641_G_A_A, data = data_sub) %>% 
    ggsurvfit(linetype_aes = TRUE, linewidth = 1) +
    labs(
      x = "Days after transplantation",
      y = "Overall survival"
    ) + 
    add_confidence_interval() +
    scale_color_grey(start = 0, end = 0.4) +
    scale_fill_grey(start = 0, end = 0.4) +
    add_censor_mark(size = 5, shape=124) +
    add_risktable(risktable_stats = c("n.risk", "cum.censor", "cum.event")) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_x_continuous(expand = c(0.02, 0.02)) +
    add_pvalue("annotation", size = 4) +
    labs(title = paste0("OS ", titles[i])) +
    theme_classic() +
    theme(legend.position="bottom")
  
  ggsave(file = paste0("./results/thymic_SNP/survival/OS_univariate_", names[i], ".png"), print(OS_plot), dpi = 600, width = 9, height = 8)
  
  logrank <- survdiff(Surv(overall_survival, OS_status) ~ chr14_22448641_G_A_A, data=data_sub)
  saveRDS(logrank, file=paste0("./results/thymic_SNP/survival/OS_univariate_", names[i], "_logrank"))
  
  #--------------------------------------------
  
  # RFS/DFS
  
  RFS_plot <- survfit2(Surv(relapse_free_survival, RFS_status) ~ chr14_22448641_G_A_A, data = data_sub) %>% 
    ggsurvfit(linetype_aes = TRUE, linewidth = 1) +
    labs(
      x = "Days after transplantation",
      y = "Disease-free survival"
    ) + 
    add_confidence_interval() +
    scale_color_grey(start = 0, end = 0.4) +
    scale_fill_grey(start = 0, end = 0.4) +
    add_censor_mark(size = 5, shape=124) +
    add_risktable(risktable_stats = c("n.risk", "cum.censor", "cum.event")) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_x_continuous(expand = c(0.02, 0.02)) +
    add_pvalue("annotation", size = 4) +
    labs(title = paste0("DFS ", titles[i])) +
    theme_classic() +
    theme(legend.position="bottom")
  ggsave(file = paste0("./results/thymic_SNP/survival/DFS_univariate_", names[i], ".png"), print(RFS_plot), dpi = 600, width = 9, height = 8)
  
  logrank <- survdiff(Surv(relapse_free_survival, RFS_status) ~ chr14_22448641_G_A_A, data=data_sub)
  saveRDS(logrank, file=paste0("./results/thymic_SNP/survival/DFS_univariate_", names[i], "_logrank"))
  
  
}

#-----------------------------------------------------------------------------------------------------------

# make the plots for OS and RFS in unrelated donors to be in one figure, side by side

files <- "./results/thymic_SNP/IDs_donors_unrelated.txt"
names <- "unrelated"
titles <- "unrelated donors"

donors <- fread(files, header = F)
data_sub <- data[data$DonorGenotypingID %in% donors$V2,]

OS_plot <- survfit2(Surv(overall_survival, OS_status) ~ chr14_22448641_G_A_A, data = data_sub) %>% 
  ggsurvfit(linetype_aes = TRUE, linewidth = 1) +
  labs(
    x = "Days after transplantation",
    y = "Overall survival"
  ) + 
  add_confidence_interval() +
  scale_color_grey(start = 0, end = 0.4) +
  scale_fill_grey(start = 0, end = 0.4) +
  add_censor_mark(size = 5, shape=124) +
  add_risktable(risktable_stats = c("n.risk", "cum.censor", "cum.event"), size = 7, theme = theme_risktable_default(axis.text.y.size = 20, plot.title.size = 20)) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(expand = c(0.04, 0.1)) +
  add_pvalue("annotation", size = 7) +
  labs(title = paste0("OS ", titles)) +
  theme_classic() +
  theme(legend.position="bottom", text = element_text(size = 20), plot.title = element_text(size = 20))

RFS_plot <- survfit2(Surv(relapse_free_survival, RFS_status) ~ chr14_22448641_G_A_A, data = data_sub) %>% 
  ggsurvfit(linetype_aes = TRUE, linewidth = 1) +
  labs(
    x = "Days after transplantation",
    y = "Relapse-free survival"
  ) + 
  add_confidence_interval() +
  scale_color_grey(start = 0, end = 0.4) +
  scale_fill_grey(start = 0, end = 0.4) +
  add_censor_mark(size = 5, shape=124) +
  add_risktable(risktable_stats = c("n.risk", "cum.censor", "cum.event"), size = 7, theme = theme_risktable_default(axis.text.y.size = 20, plot.title.size = 20)) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(expand = c(0.04, 0.1)) +
  add_pvalue("annotation", size = 7) +
  labs(title = paste0("DFS ", titles)) +
  theme_classic() +
  theme(legend.position="bottom", text = element_text(size = 20), plot.title = element_text(size = 20))

built_os <- ggsurvfit_build(OS_plot)
built_rfs <- ggsurvfit_build(RFS_plot)

combined <- cowplot::plot_grid(built_os, built_rfs, ncol = 2, labels = c("A", "B"), label_size = 20)

ggsave(file = "./results/thymic_SNP/survival/combined_OS_RFS_univariate_unrelated.png", dpi = 600, width = 20, height = 11)
ggsave(file = "./results/thymic_SNP/survival/combined_OS_RFS_univariate_unrelated.jpg", dpi = 600, width = 20, height = 11)

#-----------------------------------------------------------------------------------------------------------

# get the plots for OS and RFS in unrelated to be in one figure, side by side
# years on the x axis

files <- "./results/thymic_SNP/IDs_donors_unrelated.txt"
names <- "unrelated"
titles <- "unrelated donors"

donors <- fread(files, header = F)
data_sub <- data[data$DonorGenotypingID %in% donors$V2,]

data_sub$overall_survival_years <- data_sub$overall_survival / 365.25
data_sub$relapse_free_survival_years <- data_sub$relapse_free_survival / 365.25

OS_plot <- survfit2(Surv(overall_survival_years, OS_status) ~ chr14_22448641_G_A_A, data = data_sub) %>% 
  ggsurvfit(linetype_aes = TRUE, linewidth = 1) +
  labs(
    x = "Years after transplantation",
    y = "Overall survival"
  ) + 
  add_confidence_interval() +
  scale_color_grey(start = 0, end = 0.4) +
  scale_fill_grey(start = 0, end = 0.4) +
  add_censor_mark(size = 5, shape=124) +
  add_risktable(risktable_stats = c("n.risk", "cum.censor", "cum.event"), size = 7, theme = theme_risktable_default(axis.text.y.size = 20, plot.title.size = 20)) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(expand = c(0.04, 0.1), breaks = c(0,5,10,15,20,25)) +
  add_pvalue("annotation", size = 7) +
  labs(title = paste0("OS ", titles)) +
  theme_classic() +
  theme(legend.position="bottom", text = element_text(size = 20), plot.title = element_text(size = 20))

RFS_plot <- survfit2(Surv(relapse_free_survival_years, RFS_status) ~ chr14_22448641_G_A_A, data = data_sub) %>% 
  ggsurvfit(linetype_aes = TRUE, linewidth = 1) +
  labs(
    x = "Years after transplantation",
    y = "Relapse-free survival"
  ) + 
  add_confidence_interval() +
  scale_color_grey(start = 0, end = 0.4) +
  scale_fill_grey(start = 0, end = 0.4) +
  add_censor_mark(size = 5, shape=124) +
  add_risktable(risktable_stats = c("n.risk", "cum.censor", "cum.event"), size = 7, theme = theme_risktable_default(axis.text.y.size = 20, plot.title.size = 20)) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(expand = c(0.04, 0.1), breaks = c(0,5,10,15,20,25)) +
  add_pvalue("annotation", size = 7) +
  labs(title = paste0("DFS ", titles)) +
  theme_classic() +
  theme(legend.position="bottom", text = element_text(size = 20), plot.title = element_text(size = 20))

built_os <- ggsurvfit_build(OS_plot)
built_rfs <- ggsurvfit_build(RFS_plot)

combined <- cowplot::plot_grid(built_os, built_rfs, ncol = 2, labels = c("A", "B"), label_size = 20)

ggsave(file = "./results/thymic_SNP/survival/combined_OS_DFS_univariate_unrelated_years.jpg", dpi = 600, width = 20, height = 11)


#######################################################################################################################################################

# cox / multivariate
# with classes as factors & all levels in one column -> one set as reference automatically

data <- fread("./results/thymic_SNP/survival/survival_data_thymus_factors.txt")
data$chr14_22448641_G_A_A[data$chr14_22448641_G_A_A == "AA"] <- "Donor AA"
data$chr14_22448641_G_A_A[data$chr14_22448641_G_A_A == "GG or GA"] <- "Donor AG/GG"
data$chr14_22448641_G_A_A <- factor(data$chr14_22448641_G_A_A, levels=c("Donor AG/GG", "Donor AA")) # so that AG/GG is ref in the tables

colnames(data)[16] <- "HLA-match level"  

files <- c("./results/patient_information/ID_lists/IDs_donors.txt", "./results/thymic_SNP/IDs_donors_10_10.txt", "./results/thymic_SNP/IDs_donors_9_10.txt", "./results/thymic_SNP/IDs_donors_unrelated.txt", "./results/thymic_SNP/IDs_donors_10_10_unrelated.txt", "./results/thymic_SNP/IDs_donors_9_10_unrelated.txt")
names <- c("all_donors", "all_10_10", "all_9_10", "unrelated", "unrelated_10_10", "unrelated_9_10")
titles <- c("all donors", "all 10/10 HLA-matched donors", "all 9/10 HLA-matched donors", "unrelated donors", "unrelated 10/10 HLA-matched donors", "unrelated  9/10 HLA-matched donors")

for (i in 1:length(files)) {
  
  donors <- fread(files[i], header = F)
  data_sub <- data[data$DonorGenotypingID %in% donors$V2,]
  
  if (files[i] %like% "10") {
    
    # if a subset like 9/10 or 10/10 matched, leave the matching column out
    
    data_sub <- data_sub[,-"HLA-match level"]
    
  }
  
  print(names[i])
  # print(colnames(data_sub))
  
  # results tables
  
  # OS
  
  res_OS <- coxph(Surv(overall_survival, OS_status) ~ ., data = data_sub[,-c(1,2,4,5)]) %>% 
    tbl_regression(exp = TRUE)  %>%
    as_flex_table() 
  res_OS <- bg(res_OS, bg = "white", part = "all")
  save_as_image(res_OS, path = paste0("./results/thymic_SNP/survival/OS_multivariate_", names[i], "_factor.png"))
  
  
  
  # RFS
  
  res_RFS <- coxph(Surv(relapse_free_survival, RFS_status) ~ ., data = data_sub[,-c(1,2,3,6)]) %>% 
    tbl_regression(exp = TRUE)  %>%
    as_flex_table() 
  res_RFS <- bg(res_RFS, bg = "white", part = "all")
  save_as_image(res_RFS, path = paste0("./results/thymic_SNP/survival/RFS_multivariate_", names[i], "_factor.png"))
  
}

