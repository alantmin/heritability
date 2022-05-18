# This file is used to collect results and put them into one data frame


#Load libraries from other files
source("gcta.R")
source("genotypes_cor.R")
source("estimators.R")
library(mvtnorm)

# First set up the parameters that we will be testing

# Number of iterations per parameter combination
n_iter = 500

# Simulated true heritability value
hsq = .8

# Number of individuals n
n_vec = c()

# Number of parameters p
p_vec = c()

# Which LD structure to use
structure_vec = c()

# Parameter controlling level of LD
param_vec = c()

# Which iteration to test
iteration_vec = c()

# Create a data frame for all tested parameters

# For different combinations of n and p
for(np in list(c(1000, 200), c(200, 1000), c(200, 3000), c(2000, 1000))) {
	n = np[1]
	p = np[2]
	
	# For different structures of correlation 
	for(structure in c("autocorrelation", "repeat", "block")) {
		
		# Each structure gets its own parameters
		if(structure == "autocorrelation") {
			params = c(0, .2, .4, .6, .8)
		} else if (structure == "repeat") {
			params = c(0, 2, 4, 6, 8)
		} else if (structure == "block") {
			params = c(0, .2, .4, .6, .8)
		} else {
			print("Structure must be autocorrelation, repeat, or block")
			quit()
		}
		
		# For each different parameter 
		for(param in params) {
			
			# For each iteration
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

# Create a data frame with the different parameter combinations
settings = data.frame(n = n_vec,
					  p = p_vec,
					  structure = structure_vec,
					  param = param_vec,
					  iteration = iteration_vec)

# Go through each of the settings and retrieve the data
l = list()
counter = 1
for(i in 1:nrow(settings)) {
	n = settings$n[i]
	p = settings$p[i]
	structure = settings$structure[i]
	param = settings$param[i]
	iteration = settings$iteration[i]
	
	# Folder structures are made using the get_file_name function
	if(file.exists(get_file_name(n, p, structure, param, iteration, "heritability_results.Rda"))) {
		
		# Load the file and collect results 
		load(get_file_name(n, p, structure, param, iteration, "heritability_results.Rda"))
		l[[counter]] = df
		counter = counter + 1
	} else {
		print("The file did not exist")
		print(get_file_name(n, p, structure, param, iteration, "heritability_results.Rda"))
	}
	print(i)
}

# Combine results from each of the files 
df = do.call(rbind, l)

# Save results
save(df, file = "results.Rda")
