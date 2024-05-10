library(data.table)
library(tidyverse)

# look for SNPs 

# rs2204985
#----------------------------------------------------------

# positions for the SNPs
# from https://www.ncbi.nlm.nih.gov/snp/

# rs2204985
# 14:22448641 (GRCh38)
# 14:22917633 (GRCh37)

####################################################################################################

# check imputed files
# before info filtering (are imputed or not)

fin <- fread("./results/imputation/ic3_imputed_not_filtered/ic3_b38_imputed_info_all.bim")
eur <- fread("./results/imputation/newcastle_donors/newcastle_donors_hg38_nonreference_flipped_ambiguous_flipped_vcf_cleaned_vcf_INFO_group_chr14.txt")

fin <- fin[fin$V1 == 14 & fin$V4 == 22448641,]
#     V1                 V2 V3       V4 V5 V6
# 1: 14 chr14_22448641_G_A  0 22448641  G  A

eur <- eur[eur$SNP %like% "chr14_22448641",]
#       CHR                SNP REF ALT       AF INFO AF_GROUP
# 1: chr14 chr14_22448641_G_A   G   A 0.513736    1        3

# save the SNP name for extracting it
write.table(fin$V2, "./results/thymic_SNP/snp_name.txt", quote = F, col.names = F, row.names = F)

