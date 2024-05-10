# genotypic association for thymic SNP
# in plink1.9 files

function association() {

	genotype_file=$1
	keep_individuals=$2
	save_to=$3
	covars=${4:-./results/thymic_SNP/covars_mod.txt} # the default covar file is this, a modified file given only for poland & when all HLA-match categories present
	pheno=${5:-./results/thymic_SNP/pheno_plink.txt} # the default covar file is this, a modified file given only for poland
	
    	# association for all phenos
    	
    	# A allele recessive
    	plink2 \
		--bfile ${genotype_file}_ref_G \
		--keep ${keep_individuals} \
		--glm omit-ref recessive \
		--1 \
		--ci 0.95 \
		--covar-variance-standardize \
		--covar ${covars} \
		--pheno ${pheno} \
		--out ${save_to}_A
	
    	# G allele recessive
    	plink2 \
		--bfile ${genotype_file}_ref_A \
		--keep ${keep_individuals} \
		--glm omit-ref recessive \
		--1 \
		--ci 0.95 \
		--covar-variance-standardize \
		--covar ${covars} \
		--pheno ${pheno} \
		--out ${save_to}_G
		
	# heterozygote
	# use a newer version of plink which has this modifier (the old "default" one on this computer doesn't have it)
    	/home/nithiju/Programs/plink2_linux_avx2_20221024/plink2 \
		--bfile ${genotype_file} \
		--keep ${keep_individuals} \
		--glm hetonly \
		--1 \
		--ci 0.95 \
		--covar-variance-standardize \
		--covar ${covars} \
		--pheno ${pheno} \
		--out ${save_to}_het
}


# all donors
# all HLA match numbers

# all populations: 
association "./results/thymic_SNP/extracted_snp/all_datasets_SNP" "./results/patient_information/ID_lists/IDs_donors.txt" "./results/thymic_SNP/association_results_genotypic/plink1/allDonors/all" "./results/thymic_SNP/covars_mod_match.txt"

# fin:
association "./results/thymic_SNP/extracted_snp/finland_SNP" "./results/patient_information/ID_lists/IDs_donors.txt" "./results/thymic_SNP/association_results_genotypic/plink1/allDonors/fin" "./results/thymic_SNP/covars_mod_match.txt"

# kat :
association "./results/thymic_SNP/extracted_snp/spain_SNP" "./results/patient_information/ID_lists/IDs_donors.txt" "./results/thymic_SNP/association_results_genotypic/plink1/allDonors/kat" "./results/thymic_SNP/covars_mod_match.txt"

# nc:
association "./results/thymic_SNP/extracted_snp/uk_SNP" "./results/patient_information/ID_lists/IDs_donors.txt" "./results/thymic_SNP/association_results_genotypic/plink1/allDonors/nc" "./results/thymic_SNP/covars_mod_match.txt"

# poland
association "./results/thymic_SNP/extracted_snp/poland_SNP" "./results/patient_information/ID_lists/IDs_donors.txt" "./results/thymic_SNP/association_results_genotypic/plink1/allDonors/pol" "./results/thymic_SNP/covars_mod_pol_match.txt"

#------------------------------------
# all donors 10/10

# all populations: 
association "./results/thymic_SNP/extracted_snp/all_datasets_SNP" "./results/thymic_SNP/IDs_donors_10_10.txt" "./results/thymic_SNP/association_results_genotypic/plink1/allDonors_10_10/all"

# fin:
association "./results/thymic_SNP/extracted_snp/finland_SNP" "./results/thymic_SNP/IDs_donors_10_10.txt" "./results/thymic_SNP/association_results_genotypic/plink1/allDonors_10_10/fin"

# kat :
association "./results/thymic_SNP/extracted_snp/spain_SNP" "./results/thymic_SNP/IDs_donors_10_10.txt" "./results/thymic_SNP/association_results_genotypic/plink1/allDonors_10_10/kat"

# nc:
association "./results/thymic_SNP/extracted_snp/uk_SNP" "./results/thymic_SNP/IDs_donors_10_10.txt" "./results/thymic_SNP/association_results_genotypic/plink1/allDonors_10_10/nc"

# poland
association "./results/thymic_SNP/extracted_snp/poland_SNP" "./results/thymic_SNP/IDs_donors_10_10.txt" "./results/thymic_SNP/association_results_genotypic/plink1/allDonors_10_10/pol" "./results/thymic_SNP/covars_mod_pol.txt"

#------------------------------------
# all donors 9/10


# all populations: 
association "./results/thymic_SNP/extracted_snp/all_datasets_SNP" "./results/thymic_SNP/IDs_donors_9_10.txt" "./results/thymic_SNP/association_results_genotypic/plink1/allDonors_9_10/all"

# fin:
association "./results/thymic_SNP/extracted_snp/finland_SNP" "./results/thymic_SNP/IDs_donors_9_10.txt" "./results/thymic_SNP/association_results_genotypic/plink1/allDonors_9_10/fin"

# kat :
association "./results/thymic_SNP/extracted_snp/spain_SNP" "./results/thymic_SNP/IDs_donors_9_10.txt" "./results/thymic_SNP/association_results_genotypic/plink1/allDonors_9_10/kat"

# nc:
association "./results/thymic_SNP/extracted_snp/uk_SNP" "./results/thymic_SNP/IDs_donors_9_10.txt" "./results/thymic_SNP/association_results_genotypic/plink1/allDonors_9_10/nc"

# poland
association "./results/thymic_SNP/extracted_snp/poland_SNP" "./results/thymic_SNP/IDs_donors_9_10.txt" "./results/thymic_SNP/association_results_genotypic/plink1/allDonors_9_10/pol" "./results/thymic_SNP/covars_mod_pol.txt" "./results/thymic_SNP/pheno_plink_mod_2.txt"

#------------------------------------
# unrelated donors
# all HLA match numbers

# all populations: 
association "./results/thymic_SNP/extracted_snp/all_datasets_SNP" "./results/thymic_SNP/IDs_donors_unrelated.txt" "./results/thymic_SNP/association_results_genotypic/plink1/unrelatedDonors/all" "./results/thymic_SNP/covars_mod_match.txt"

# fin:
association "./results/thymic_SNP/extracted_snp/finland_SNP" "./results/thymic_SNP/IDs_donors_unrelated.txt" "./results/thymic_SNP/association_results_genotypic/plink1/unrelatedDonors/fin" "./results/thymic_SNP/covars_mod_match.txt"

# kat :
association "./results/thymic_SNP/extracted_snp/spain_SNP" "./results/thymic_SNP/IDs_donors_unrelated.txt" "./results/thymic_SNP/association_results_genotypic/plink1/unrelatedDonors/kat" "./results/thymic_SNP/covars_mod_match.txt"

# nc:
association "./results/thymic_SNP/extracted_snp/uk_SNP" "./results/thymic_SNP/IDs_donors_unrelated.txt" "./results/thymic_SNP/association_results_genotypic/plink1/unrelatedDonors/nc" "./results/thymic_SNP/covars_mod_match.txt"

# poland
association "./results/thymic_SNP/extracted_snp/poland_SNP" "./results/thymic_SNP/IDs_donors_unrelated.txt" "./results/thymic_SNP/association_results_genotypic/plink1/unrelatedDonors/pol" "./results/thymic_SNP/covars_mod_pol_match.txt"

#------------------------------------
# unrelated donors 10/10

# all populations: 
association "./results/thymic_SNP/extracted_snp/all_datasets_SNP" "./results/thymic_SNP/IDs_donors_10_10_unrelated.txt" "./results/thymic_SNP/association_results_genotypic/plink1/unrelatedDonors_10_10/all"

# fin:
association "./results/thymic_SNP/extracted_snp/finland_SNP" "./results/thymic_SNP/IDs_donors_10_10_unrelated.txt" "./results/thymic_SNP/association_results_genotypic/plink1/unrelatedDonors_10_10/fin"

# kat :
association "./results/thymic_SNP/extracted_snp/spain_SNP" "./results/thymic_SNP/IDs_donors_10_10_unrelated.txt" "./results/thymic_SNP/association_results_genotypic/plink1/unrelatedDonors_10_10/kat"

# nc:
association "./results/thymic_SNP/extracted_snp/uk_SNP" "./results/thymic_SNP/IDs_donors_10_10_unrelated.txt" "./results/thymic_SNP/association_results_genotypic/plink1/unrelatedDonors_10_10/nc"

# poland
association "./results/thymic_SNP/extracted_snp/poland_SNP" "./results/thymic_SNP/IDs_donors_10_10_unrelated.txt" "./results/thymic_SNP/association_results_genotypic/plink1/unrelatedDonors_10_10/pol" "./results/thymic_SNP/covars_mod_pol.txt"

#------------------------------------
# unrelated donors 9/10

# all populations: 
association "./results/thymic_SNP/extracted_snp/all_datasets_SNP" "./results/thymic_SNP/IDs_donors_9_10_unrelated.txt" "./results/thymic_SNP/association_results_genotypic/plink1/unrelatedDonors_9_10/all"

# fin:
association "./results/thymic_SNP/extracted_snp/finland_SNP" "./results/thymic_SNP/IDs_donors_9_10_unrelated.txt" "./results/thymic_SNP/association_results_genotypic/plink1/unrelatedDonors_9_10/fin"

# kat :
association "./results/thymic_SNP/extracted_snp/spain_SNP" "./results/thymic_SNP/IDs_donors_9_10_unrelated.txt" "./results/thymic_SNP/association_results_genotypic/plink1/unrelatedDonors_9_10/kat"

# nc:
association "./results/thymic_SNP/extracted_snp/uk_SNP" "./results/thymic_SNP/IDs_donors_9_10_unrelated.txt" "./results/thymic_SNP/association_results_genotypic/plink1/unrelatedDonors_9_10/nc"

# poland
association "./results/thymic_SNP/extracted_snp/poland_SNP" "./results/thymic_SNP/IDs_donors_9_10_unrelated.txt" "./results/thymic_SNP/association_results_genotypic/plink1/unrelatedDonors_9_10/pol" "./results/thymic_SNP/covars_mod_pol.txt" "./results/thymic_SNP/pheno_plink_mod_2.txt"





