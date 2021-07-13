source("gcta.R")
source("genotypes_cor.R")
source("estimators.R")

library(mvtnorm)
n_iter = 500
hsq = .8

n_vec = c()
p_vec = c()
structure_vec = c()
param_vec = c()
iteration_vec = c()
for(np in list(c(1000, 200), c(200, 1000), c(200, 3000), c(2000, 1000))) {
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

l = list()
counter = 1
for(i in 1:nrow(settings)) {
	n = settings$n[i]
	p = settings$p[i]
	structure = settings$structure[i]
	param = settings$param[i]
	iteration = settings$iteration[i]
	
	if(file.exists(get_file_name(n, p, structure, param, iteration, "alan_results.Rda"))) {
		load(get_file_name(n, p, structure, param, iteration, "alan_results.Rda"))
		l[[counter]] = df
		counter = counter + 1
	} else {
		print("WOW THIS DOESN'T EXIST!")
		print(get_file_name(n, p, structure, param, iteration, "alan_results.Rda"))
	}
	print(i)
}
df = do.call(rbind, l)
