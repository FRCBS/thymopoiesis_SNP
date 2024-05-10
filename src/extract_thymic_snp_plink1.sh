# extract the thuýmus snp from all populations

all_datasets=("/media/volume/datavolume/harmonizing/all_datasets/merged_all_datasets_all_individuals" "/media/volume/datavolume/harmonizing/all_datasets/merged_finns_all_individuals" "/media/volume/datavolume/imputation/katalonia/katalonia_all_chrs_info_filtered_QC_variantsRenamed" "/media/volume/datavolume/harmonizing/all_datasets/merged_newcastle_all_individuals" "/media/volume/datavolume/imputation/poland/poland_all_chrs_info_filtered_QC_QC2_variantsRenamed")

save_here=("all_datasets" "finland" "spain" "uk" "poland")

for i in 0 1 2 3 4
do
	
	DATASET=${all_datasets[i]}
	save=${save_here[i]}
	
	plink \
		--bfile ${DATASET} \
		--extract /media/volume/datavolume/thymic_SNP/snp_name.txt \
		--make-bed \
		--out /media/volume/datavolume/thymic_SNP/${save}_SNP
	
done
