source("gcta.R")
source("genotypes_cor.R")
source("estimators.R")
source("genotypes_related.R")
library(rres)

library(mvtnorm)
# Optionally set a seed 
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

#Create a list of different parameters to enter into analysis
n_vec = rep(0, 80000)
p_vec = rep(0, 80000)
iteration_vec = rep(0, 80000)
structure_vec = rep(0, 80000)
param = 0
counter = 1

# For the number of markers between 400 and 8000
for(p in seq(400, 8000, by = 400)) {
	
	# The number of individuals is set to 400
    n = 400
    
    # Set different levels of cousins: unrelated, 1st, 2nd, and 3rd cousins
    for(cousinship in c(0,1,2,3)) {
        for(iteration in 1:n_iter) {
            n_vec[counter] = n
            p_vec[counter] = p
            structure_vec[counter] = cousinship
            iteration_vec[counter] = iteration
            counter = counter + 1
        }
    }
}
settings = data.frame(n = n_vec,
                      p = p_vec,
                      structure = paste("cousin", structure_vec, sep = ""),
                      iteration = iteration_vec)

# Select 100 parameter sets to run for this batch
array_index = as.integer(commandArgs(trailingOnly = T))[1]
indexes = 1:100 + (array_index - 1) * 100
for(i in indexes) {
	n = settings$n[i]
	p = settings$p[i]
	structure = settings$structure[i]
	iteration = settings$iteration[i]
	
	# Get genotypes and phenotypes for unrelated individuals 
	if (structure == "cousin0") {
	    g = get_genotypes(n, p, 0, 0, mafs = sampled_markers_subpop$afr_alt[1:p])$genotypes
	    pheno = get_phenotypes(g, hsq)
	 
	# Get genotypes and phenotypes for first cousins
	} else if(structure == "cousin1") {
	    g = get_related_genotypes_1st_cous(n, p, sampled_markers_subpop$afr_alt[1:p])
	    pheno = get_phenotypes(g, hsq)
	
	# Get genotypes and phenotypes for second cousins
	} else if (structure == "cousin2") {
	    g = get_related_genotypes_2nd_cous(n, p, sampled_markers_subpop$afr_alt[1:p])
	    pheno = get_phenotypes(g, hsq)
	    
	# Get genotypes and phenotypes for third cousins
	} else if (structure == "cousin3") {
	    g = get_related_genotypes_3rd_cous(n, p, sampled_markers_subpop$afr_alt[1:p])   
	    pheno = get_phenotypes(g, hsq)
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
		
		# Create a .bed file using plink for GCTA
		write.table(phen_file, file = get_file_name(n, p, structure, param, iteration, "test.phen"), row.names = F, col.names = F, quote = F)
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
		dicker1 = Re(dicker_2014_estimator_1_sd(g, pheno$phenotypes))
		
		# Calculate the schwartzman estimator
		schwartzman = schwartzman_2019_estimator(g, pheno$phenotypes)
		
		# Calculate the haseman elston MoM estimator without including
		# the diagonal
		MoM_no_diag_result = MoM_no_diag(g, pheno$phenotypes)
		
		# Calculate the haseman elston MoM estimator with the diagonal
		MoM_diag_result = MoM_diag(g, pheno$phenotypes)
		
		# Calculate the GCTA estimator 
		gcta_maxlike = gcta(folder)
		
		# Calculate the LDAK estimator
		ldak_res = ldak(folder)
		
		# Create a data frame with each of the estimators
		# Where applicable, keep estimates of sigma_g, sigma_e, and hsq
		names = c("gold_standard", "dicker2", "dicker1", "schwartzman", "MoM_no_diag", "MoM_diag", "gcta", "ldak")
		result_hsq = c(gold_standard[1], dicker2[1], dicker1[1], schwartzman[1], MoM_no_diag_result[1], MoM_diag_result[1], gcta_maxlike[1], ldak_res[1])
		sigmag = c(gold_standard[2], dicker2[2], dicker1[2], schwartzman[2], MoM_no_diag_result[2], MoM_diag_result[2], gcta_maxlike[2], ldak_res[2])
		sigmae = c(gold_standard[3], dicker2[3], dicker1[3], schwartzman[3], MoM_no_diag_result[3], MoM_diag_result[3], gcta_maxlike[3], ldak_res[3])
		df = data.frame(names, result_hsq, sigmae, sigmag, n, p, structure, param, iteration)
		if(clean_up) {
			system(paste("/bin/rm -r", get_file_name(n, p, structure, param, iteration, "*")))
		}
		
		# Save file
		save(df,file = get_file_name(n, p, structure, param, iteration, "alan_results.Rda")) 
	}
}


