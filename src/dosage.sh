# make dosage for the thymus SNP

plink2 \
	--bfile ./results/thymic_SNP/extracted_snp/all_datasets_SNP \
	--recode A \
	--out ./results/thymic_SNP/survival/dosage


