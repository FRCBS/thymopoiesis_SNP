library(meta)
library(tidyverse)
library(data.table)

#########################################################################################################################

# meta-analysis on the thymus SNP

# tsamadou
# all unrelated

# OS multivariate 9/10 HLA match
# Donor rs2204985: AA, HR 1.48 (1.15, 1.92), p-value 0.003

# DFS multivariate 9/10 HLA match
# Donor rs2204985: AA, HR 1.50 (1.18, 1.91), p-value 0.001


# supple:
# OS multivariate 10/10 HLA match
# Donor rs2204985: AA, HR 0.97 (0.79, 1.19), p-value 0.771

# DFS multivariate 10/10 HLA match
# Donor rs2204985: AA, HR 1.03 (0.86, 1.24), p-value 0.749


#-------------------------------------------------------------------------------------------------

# our analysis

# unrelated, 9/10 HLA match, multivariate, OS: AA, HR 1.24 (0.34, 4.50)
# unrelated, 9/10 HLA match, multivariate, RFS: AA, HR 0.83 (0.18, 3.85)

# unrelated, 10/10 HLA match, multivariate, OS: AA, HR 1.36 (0.95, 1.96)
# unrelated, 10/10 HLA match, multivariate, RFS: AA, HR 1.19 (0.76, 1.86)

#-------------------------------------------------------------------------------------------

# unrelated, 9/10 match
# OS

study <- c("Tsamadou", "Our analysis")
HR <- c(1.48, 1.24)
lower.HR <- c(1.15, 0.34)
upper.HR <- c(1.92, 4.50)

m <- metagen(
  HR = log(HR), lower = log(lower.HR), upper = log(upper.HR),
  sm = "HR", fixed=F, random=T,
  studlab = study,
  method.tau = "DL", ## method to calculate Tau
  method.random.ci = "classic", ## method to calculate estimator's CI
)

summary(m)
#                   HR           95%-CI %W(random)
# Tsamadou     1.4859 [1.1500; 1.9200]       96.2
# Our analysis 1.2369 [0.3400; 4.5000]        3.8
# 
# Number of studies: k = 2
# 
# HR           95%-CI    z p-value
# Random effects model 1.4756 [1.1476; 1.8974] 3.03  0.0024
# 
# Quantifying heterogeneity:
#   tau^2 = 0; tau = 0; I^2 = 0.0%; H = 1.00
# 
# Test of heterogeneity:
#   Q d.f. p-value
# 0.07    1  0.7848
# 
# Details on meta-analytical method:
#   - Inverse variance method
# - DerSimonian-Laird estimator for tau^2

png(file = "./results/thymic_SNP/meta_analysis/unrelated_9_10_OS_forestplot.png", width = 2800, height = 1000, res = 300)
forest(m, xlim = c(0.3,4.5))
dev.off()

# funnel.meta(m, studlab = TRUE)

#--------------------------------------------------------------------------------

# DFS

study <- c("Tsamadou", "Our analysis")
HR <- c(1.50, 0.83)
lower.HR <- c(1.18, 0.18)
upper.HR <- c(1.91, 3.85)
 

m <- metagen(
  HR = log(HR), lower = log(lower.HR), upper = log(upper.HR),
  sm = "HR", fixed=F, random=T,
  studlab = study,
  method.tau = "DL", ## method to calculate Tau
  method.random.ci = "classic", ## method to calculate estimator's CI
)

summary(m)
#                   HR           95%-CI %W(random)
# Tsamadou     1.5013 [1.1800; 1.9100]       97.6
# Our analysis 0.8325 [0.1800; 3.8500]        2.4
# 
# Number of studies: k = 2
# 
# HR           95%-CI    z p-value
# Random effects model 1.4801 [1.1667; 1.8775] 3.23  0.0012
# 
# Quantifying heterogeneity:
#   tau^2 = 0; tau = 0; I^2 = 0.0%; H = 1.00
# 
# Test of heterogeneity:
#   Q d.f. p-value
# 0.56    1  0.4560
# 
# Details on meta-analytical method:
#   - Inverse variance method
# - DerSimonian-Laird estimator for tau^2


png(file = "./results/thymic_SNP/meta_analysis/unrelated_9_10_DFS_forestplot.png", width = 2800, height = 1000, res = 300)
forest(m, xlim = c(0.3,4.5))
dev.off()


#############################################################################################################

# unrelated, 10/10 match
# OS

study <- c("Tsamadou", "Our analysis")
HR <- c(0.97, 1.36)
lower.HR <- c(0.79, 0.95)
upper.HR <- c(1.19, 1.96)

m <- metagen(
  HR = log(HR), lower = log(lower.HR), upper = log(upper.HR),
  sm = "HR", fixed=F, random=T,
  studlab = study,
  method.tau = "DL", ## method to calculate Tau
  method.random.ci = "classic", ## method to calculate estimator's CI
)

summary(m)

png(file = "./results/thymic_SNP/meta_analysis/unrelated_10_10_OS_forestplot.png", width = 2800, height = 1000, res = 300)
forest(m, xlim = c(0.3,4.5))
dev.off()



# DFS

study <- c("Tsamadou", "Our analysis")
HR <- c(1.03, 1.19)
lower.HR <- c(0.86, 0.76)
upper.HR <- c(1.24, 1.86)

m <- metagen(
  HR = log(HR), lower = log(lower.HR), upper = log(upper.HR),
  sm = "HR", fixed=F, random=T,
  studlab = study,
  method.tau = "DL", ## method to calculate Tau
  method.random.ci = "classic", ## method to calculate estimator's CI
)

summary(m)

png(file = "./results/thymic_SNP/meta_analysis/unrelated_10_10_DFS_forestplot.png", width = 2800, height = 1000, res = 300)
forest(m, xlim = c(0.3,4.5))
dev.off()

#######################################################################################################################

# our result form unrelated all (=9/10+10/10+ others), but this not in tsamadou
# -> try combining their 9/10, their 10/10, our 9/10, our 10/10
# (even though this is not exactly the set in which our results were in)

# OS

study <- c("Tsamadou 9/10", "Our analysis 9/10", "Tsamadou 10/10", "Our analysis 10/10")
HR <- c(1.48, 1.24, 0.97, 1.36)
lower.HR <- c(1.15, 0.34, 0.79, 0.95)
upper.HR <- c(1.92, 4.50, 1.19, 1.96)

m <- metagen(
  HR = log(HR), lower = log(lower.HR), upper = log(upper.HR),
  sm = "HR", fixed=F, random=T,
  studlab = study,
  method.tau = "DL", ## method to calculate Tau
  method.random.ci = "classic", ## method to calculate estimator's CI
)

summary(m)

png(file = "./results/thymic_SNP/meta_analysis/unrelated_910_10_OS_forestplot.png", width = 2800, height = 1000, res = 300)
forest(m, xlim = c(0.3,4.5))
dev.off()

#----------------------------------------------------------------------------------------------------------------------------
# -> try combining their 9/10, their 10/10, our all unrelated (= where our result was)

# OS

study <- c("Tsamadou 9/10", "Tsamadou 10/10", "Our analysis all unrelated")
HR <- c(1.48, 0.97, 1.34)
lower.HR <- c(1.15, 0.79, 0.96)
upper.HR <- c(1.92, 1.19, 1.86)

m <- metagen(
  HR = log(HR), lower = log(lower.HR), upper = log(upper.HR),
  sm = "HR", fixed=F, random=T,
  studlab = study,
  method.tau = "DL", ## method to calculate Tau
  method.random.ci = "classic", ## method to calculate estimator's CI
)

summary(m)


# summary(m)
png(file = "./results/thymic_SNP/meta_analysis/unrelated_all_OS_forestplot.png", width = 2800, height = 1000, res = 300)
forest(m, xlim = c(0.3,4.5))
dev.off()






