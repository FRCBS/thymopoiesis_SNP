# a test for association testing in spanish data

# A allele recessive
plink2 \
	--bfile ./results/thymic_SNP/extracted_snp/spain_SNP_ref_G \
	--keep ./results/thymic_SNP/IDs_donors_9_10.txt \
	--glm omit-ref recessive \
	--1 \
	--ci 0.95 \
	--pheno ./results/thymic_SNP/pheno_plink.txt \
	--out ./results/thymic_SNP/association_results_genotypic/plink1/kat_test_A

# G allele recessive
plink2 \
	--bfile ./results/thymic_SNP/extracted_snp/spain_SNP_ref_A \
	--keep ./results/thymic_SNP/IDs_donors_9_10.txt \
	--glm omit-ref recessive \
	--1 \
	--ci 0.95 \
	--pheno ./results/thymic_SNP/pheno_plink.txt \
	--out ./results/thymic_SNP/association_results_genotypic/plink1/kat_test_G
	
# heterozygote
# use a newer version of plink which has this modifier (the old "default" one on this computer doesn't have it)
/home/nithiju/Programs/plink2_linux_avx2_20221024/plink2 \
	--bfile ./results/thymic_SNP/extracted_snp/spain_SNP \
	--keep ./results/thymic_SNP/IDs_donors_9_10.txt \
	--glm hetonly allow-no-covars \
	--1 \
	--ci 0.95 \
	--pheno ./results/thymic_SNP/pheno_plink.txt \
	--out ./results/thymic_SNP/association_results_genotypic/plink1/kat_test_het


