# This file is used to collect results and put them into one data frame for the cousins


source("gcta.R")
source("genotypes_cor.R")
source("estimators.R")
library(mvtnorm)

# First set up the parameters that we will be testing

# Number of iterations per parameter combination
n_iter = 1000

# Simulated true heritability value
hsq = .8
clean_up = T

# Number of individuals n
n_vec = rep(0, 80000)

# Number of parameters p
p_vec = rep(0, 80000)

# Which iteration to test
iteration_vec = rep(0, 80000)

# Which cousinship structure to use
structure_vec = rep(0, 80000)
param = 0
counter = 1

# Create a data frame for all tested parameters

# For different values of number of markers
for(p in seq(400, 8000, by = 400)) {
	n = 400
	
	# For different levels of cousinship
	for(cousinship in c(0,1,2,3)) {
		
		# For each iteration
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

l = list()
counter = 1
# Go through each of the settings and retrieve the data
for(i in 1:nrow(settings)) {
	n = settings$n[i]
	p = settings$p[i]
	structure = settings$structure[i]
	param = 0
	iteration = settings$iteration[i]
	
	# Folder structures are made using the get_file_name function
	if(file.exists(get_file_name(n, p, structure, param, iteration, "heritability_results.Rda"))) {
		load(get_file_name(n, p, structure, param, iteration, "heritability_results.Rda"))
		l[[counter]] = df
		counter = counter + 1
	} else {
		print("There was a missing file at: ")
		print(get_file_name(n, p, structure, param, iteration, "heritability_results.Rda"))
	}
	print(i)
}

# Combine results from each of the files 
df = do.call(rbind, l)

# Save results
save(df, "results_cousins.Rda")