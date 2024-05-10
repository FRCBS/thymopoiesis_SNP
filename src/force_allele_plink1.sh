# force both alleles

# first allele

all_datasets=("./results/thymic_SNP/extracted_snp/all_datasets_SNP" "./results/thymic_SNP/extracted_snp/finland_SNP" "./results/thymic_SNP/extracted_snp/spain_SNP" "./results/thymic_SNP/extracted_snp/uk_SNP" "./results/thymic_SNP/extracted_snp/poland_SNP" )

for i in 0 1 2 3 4
do
	
	DATASET=${all_datasets[i]}
	
	plink2 \
		--bfile ${DATASET} \
		--ref-allele force ./results/thymic_SNP/change_ref_allele.txt 1 2 \
		--make-bfile \
		--out ${DATASET}_ref_G

		
	
done 

# second allele

all_datasets=("./results/thymic_SNP/extracted_snp/all_datasets_SNP" "./results/thymic_SNP/extracted_snp/finland_SNP" "./results/thymic_SNP/extracted_snp/spain_SNP" "./results/thymic_SNP/extracted_snp/uk_SNP" "./results/thymic_SNP/extracted_snp/poland_SNP" )

for i in 0 1 2 3 4
do
	
	DATASET=${all_datasets[i]}
	
	plink2 \
		--bfile ${DATASET} \
		--ref-allele force ./results/thymic_SNP/change_ref_allele_A.txt 1 2 \
		--make-bfile \
		--out ${DATASET}_ref_A

		
	
done 


