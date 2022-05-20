# This file takes in an index and runs a batch of the specified settings
# It is assumed that plink, gcta, and LDAK are installed. 

source("gcta.R")
source("genotypes_cor.R")
source("estimators.R")

library(mvtnorm)

# Optionally set a seed before running
#set.seed(23)

# Load SNP information from 1000 genomes data, where the order of the markers is random
# The AFR population's alternate allele frequencies will be taken
# Alternate allele frequencies are filtered to be > 0.05 
sampled_markers_subpop <- read.csv("markers1shuffle.txt", sep="", stringsAsFactors=FALSE)

# Number of iterations
n_iter = 500

# Set the proportion of markers to be repeated for the repeat structure
repeat_proportion = .1

# Set the true simulation heritability
hsq = .8

# If set to true, remove intermediate files
clean_up = T

# Create a list of parameters for analyses

n_vec = c()
p_vec = c()
structure_vec = c()
param_vec = c()
iteration_vec = c()
for(np in list(c(1000, 200), c(200, 1000), c(2000, 1000), c(200, 3000))) {
	n = np[1]
	p = np[2]
	for(structure in c("autocorrelation", "repeat", "block")) {
		if(structure == "autocorrelation") {
			params = c(0, .2, .4, .6, .8)
		} else if (structure == "repeat") {
			params = c(0, 2, 4, 6, 8)
		} else if (structure == "block") {
			params = c(0, .2, .4, .6, .8)
		} else {
			print("Uhm something went wrong!")
			quit()
		}
		for(param in params) {
			for(iteration in 1:n_iter) {
				n_vec = c(n_vec, n)
				p_vec = c(p_vec, p)
				structure_vec = c(structure_vec, structure)
				param_vec = c(param_vec, param)
				iteration_vec = c(iteration_vec, iteration)
			}
		}
		
	}
}
settings = data.frame(n = n_vec,
					  p = p_vec,
					  structure = structure_vec,
					  param = param_vec,
					  iteration = iteration_vec)

# Set the indices of the parameters to test for this batch
# This is currently hard coded to run 50 parameters per batch
array_index = as.integer(commandArgs(trailingOnly = T))[1]
indexes = 1:50 + (array_index - 1) * 50 

# For each chosen index 
for(i in indexes) {
	
	# Get the parameters for this batch
	n = settings$n[i]
	p = settings$p[i]
	structure = settings$structure[i]
	param = settings$param[i]
	iteration = settings$iteration[i]
	
	# Compute different genotypes depending on the structure
	if(structure == "autocorrelation") {
		
		# Autocorrelated markers are simulated with n individuals, p markers
		# Level of autocorrelation is set using param
		# Markers are set to be binarized
		# Minor allele frequency set to be minor allele frequencies from the AFR population
		g = get_genotypes_autocor(n, p, param, T, sampled_markers_subpop$afr_alt[1:p])
		
		# Phenotypes are simulated by taking every other simulated genotype 
		# True heritability is set up to be hsq 
		pheno = get_phenotypes(g[,seq(1,p,by = 2)], hsq)
		
	} else if (structure == "repeat") {
		
		# Repeated markers are simulated with n individuals, p markers
		# p * repeat_proportion markers are repeated for param times
		# Minor allele frequency set to be minor allele frequencies from the AFR population 
		g_list = get_genotypes(n, p, round(p * repeat_proportion), param, sampled_markers_subpop$afr_alt[1:p])
		
		# Phenotypes are set using the markers not being repeated
		pheno = get_phenotypes(g_list$genotypes, hsq)
		
		# Genotypes are set to be the non-repeated markers plus repeated markers
		g = cbind(g_list$genotypes, g_list$genotypes_rep)
	} else if (structure == "block") {
		
		# Block-correlated markers are simulated with n individuals, p markers
		# Level of correlation is set using param
		# Markers are set to be binarized
		# Minor allele frequency set to be minor allele frequencies from the AFR population
		g = get_genotypes_block(n, p, param, T, sampled_markers_subpop$afr_alt[1:p])
		pheno = get_phenotypes(g[,seq(1,p,by = 2)], hsq)
	} else {
		print("Invalid structure option")
		quit()
	}
	
	if(!is.null(g)) {
		# Extract the phenotypes with noise 
		pheno_with_noise = pheno$phenotypes
		
		# Get a folder name to do analyses and store results
		folder = get_folder(n, p, structure, param, iteration)
		
		# Make the folder 
		system(paste("mkdir -p", folder))
		
		# Save the genotypes and phenotypes in the folder
		save(g, file = get_file_name(n, p, structure, param, iteration, "genotypes.Rda"))
		save(pheno, file = get_file_name(n, p, structure, param, iteration, "phenotypes.Rda"))
		
		# GCTA requires files in a particular format, so we take some steps to make these files
		
		# Make up marker names for GCTA
		map_file = get_samp_markers(ncol(g))
		
		# Write the MAP file 
		write.table(map_file, file = get_file_name(n, p, structure, param, iteration, "test.map"), row.names = F, col.names = F, quote = F)
		
		# Create a pedigree file for GCTA. This makes a separate family for each individual
		ped_file = get_ped_file(g, p, map_file)
		write.table(ped_file, file = get_file_name(n, p, structure, param, iteration, "test.ped"), row.names = F, col.names = F, quote = F)
		
		# Create the phenotype file for GCTA 
		phen_file = get_phen_file(pheno_with_noise)
		write.table(phen_file, file = get_file_name(n, p, structure, param, iteration, "test.phen"), row.names = F, col.names = F, quote = F)
		
		# Create a .bed file using plink for GCTA
		system(paste("plink --file", get_file_name(n, p, structure, param, iteration, "test"), "--make-bed --out", get_file_name(n, p, structure, param, iteration, "test")))
		print(paste("Just finished ", get_file_name(n, p, structure, param, iteration, "test")))
		
		
		# Do analysis
		
		# Get file name
		folder = get_folder(n, p, structure, param, iteration)
		
		# Load genotypes and phenotypes that we created above
		load(get_file_name(n, p, structure, param, iteration, "genotypes.Rda"))
		load(get_file_name(n, p, structure, param, iteration, "phenotypes.Rda"))
		
		# The gold standard is the varnace of the true phenotypes
		# divided by the variance of  phenotypes with noise
		gold_standard = var(pheno$true_phenotypes)/var(pheno$phenotypes)
		
		# Calculate the dicker 2 estimator 
		dicker2 = dicker_2014_estimator_2_sd(g, pheno$phenotypes)
		
		# Calculate the dicker 1 estimator 
		# dicker1 = Re(dicker_2014_estimator_1_sd(g, pheno$phenotypes))
		dicker1 = dicker_2014_estimator_1_sd(g, pheno$phenotypes)
		
		# Calculate the schwartzman estimator
		schwartzman = schwartzman_2019_estimator(g, pheno$phenotypes)
		
		# Calculate the haseman elston MoM estimator without including
		# the diagonal
		MoM_no_diag_result = MoM_no_diag(g, pheno$phenotypes)
		
		# Calculate the haseman elston MoM estimator with the diagonal
		MoM_diag_result = MoM_diag(g, pheno$phenotypes)
		
		# Calculate the GCTA estimator 
		gcta_maxlike = gcta(folder)
		
		# Calculate the LDAK estimator with two options for the power parameter.
		ldak_res_1 = ldak(folder, power=-.25)
		ldak_res_2 = ldak(folder, power=-1)
		
		# Create a data frame with each of the estimators
		# Where applicable, keep estimates of sigma_g, sigma_e, and hsq
		names = c("gold_standard", "dicker2", "dicker1", "schwartzman", "MoM_no_diag", "MoM_diag", "gcta", "ldak1", "ldak2")
		result_hsq = c(gold_standard[1], dicker2[1], dicker1[1], schwartzman[1], MoM_no_diag_result[1], MoM_diag_result[1], gcta_maxlike[1], ldak_res_1[1], ldak_res_2[1])
		sigmag = c(gold_standard[2], dicker2[2], dicker1[2], schwartzman[2], MoM_no_diag_result[2], MoM_diag_result[2], gcta_maxlike[2], ldak_res_1[2], ldak_res_2[2])
		sigmae = c(gold_standard[3], dicker2[3], dicker1[3], schwartzman[3], MoM_no_diag_result[3], MoM_diag_result[3], gcta_maxlike[3], ldak_res_1[3], ldak_res_2[3])
		df = data.frame(names, result_hsq, sigmae, sigmag, n, p, structure, param, iteration)
		if(clean_up) {
			system(paste("rm", get_file_name(n, p, structure, param, iteration, "*")))
		}
		
		# Save file 
		save(df,file = get_file_name(n, p, structure, param, iteration, "heritability_results.Rda")) 
	}
}


