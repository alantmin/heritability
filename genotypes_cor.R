library(mvtnorm)
autocorr.mat <- function(n_marker, rho) {
	mat <- diag(n_marker)
	return(rho^abs(row(mat)-col(mat)))
}

block.mat = function(n_marker, rho) {
	mat = matrix(rho, nrow = n_marker, ncol = n_marker)
	mat = mat + diag(1-rho, n_marker)
	return(mat)
}

#Mafs is a vector of MAFs that you want to have. Defaults to .5 for everything if not.
get_genotypes_autocor = function(n_indiv, n_markers, p, set_binomial = T, mafs = NULL) {
	sigma = autocorr.mat(n_markers, p)
	genotypes = rmvnorm(n = n_indiv, mean = rep(0, n_markers), sigma = sigma)
	
	if(set_binomial) {
		if(!is.null(mafs)) {
			if(length(mafs) != n_markers) {
				print("Length of mafs not equal to n_markers!")
				return(0)
			}
			lowerbounds = qnorm(mafs ^ 2)
			upperbounds = qnorm(2* mafs * (1-mafs) + mafs^2)
			lowerbounds = as.vector(matrix(rep(lowerbounds,n_indiv),
								 ncol = length(lowerbounds),
								 byrow = T))
			upperbounds = as.vector(matrix(rep(upperbounds,n_indiv),
								 ncol = length(upperbounds),
								 byrow = T))
			genvec = as.vector(genotypes)
			
			genvec = sapply(1:length(genvec), function(x) {
				if( genvec[x] < lowerbounds[x]) {
					return(0)
				} else if (genvec[x] < upperbounds[x]) {
					return(1)
				} else {
					return(2)
				}})
			genotypes = matrix(genvec, nrow = n_indiv, ncol = n_markers)
		} else {
			lowerbound = qnorm(.5 ^ 2)
			upperbound = qnorm(2* .5 * (1-.5) + .5^2)
			genvec = sapply(as.vector(genotypes), function(x) {
				if(x < lowerbound) {
					return(0)
				} else if (x < upperbound) {
					return(1)
				} else {
					return(2)
				}})
			genotypes = matrix(genvec, nrow = n_indiv, ncol = n_markers)
		}
		
	}
	return(genotypes)
}


#N_indiv is the number of individuals
#n_markers is the number of markers
#p is the parameter that defines the amount of correlation between members of a block
#set_binomial: T indicates that we convert the gaussian data into binomial data
#Mafs: Give a list of mafs that we want the markers to conform to. 
get_genotypes_block = function(n_indiv, n_markers, p, set_binomial = T, mafs = NULL, n_blocks = 10) {
	if(n_markers %% n_blocks != 0) {
		print("N markers was not divisible by n_blocks")
		quit()
		return(0)
	}
	
	#Make a big list of these fully correlated markers
	sigma = block.mat(n_markers/n_blocks, p)
	l = lapply(1:n_blocks, function(i) return( rmvnorm(n = n_indiv, mean = rep(0, n_markers/n_blocks), sigma = sigma)))
	genotypes = do.call(cbind, l)
	
	if(set_binomial) {
		if(!is.null(mafs)) {
			if(length(mafs) != n_markers) {
				print("Length of mafs not equal to n_markers!")
				return(0)
			}
			lowerbounds = qnorm(mafs ^ 2)
			upperbounds = qnorm(2* mafs * (1-mafs) + mafs^2)
			lowerbounds = as.vector(matrix(rep(lowerbounds,n_indiv),
										   ncol = length(lowerbounds),
										   byrow = T))
			upperbounds = as.vector(matrix(rep(upperbounds,n_indiv),
										   ncol = length(upperbounds),
										   byrow = T))
			genvec = as.vector(genotypes)
			
			genvec = sapply(1:length(genvec), function(x) {
				if( genvec[x] < lowerbounds[x]) {
					return(0)
				} else if (genvec[x] < upperbounds[x]) {
					return(1)
				} else {
					return(2)
				}})
			genotypes = matrix(genvec, nrow = n_indiv, ncol = n_markers)
		} else {
			stop("Error in genotypes_cor.R: This is messed up")
			genotypes = genotypes/2 + 1
			genvec = sapply(as.vector(genotypes), function(x) {
				if(x < .5) {
					return(0)
				} else if (x < 1.5) {
					return(1)
				} else {
					return(2)
				}})
			genotypes = matrix(genvec, nrow = n_indiv, ncol = n_markers)
		}
		
	}
	return(genotypes)
}

get_phenotypes_cor = function(genotypes, hsq) {
	n_indiv = nrow(genotypes)
	n_marker = ncol(genotypes)
	
	genotypes_normalized = genotypes
	for (i in 1:ncol(genotypes_normalized)) {
		genotypes_normalized[,i] = (genotypes[,i] - mean(genotypes[,i]))/sd(genotypes[,i])
	}
	
	betas = rnorm(n_marker, mean = 0, sd = sqrt(hsq/n_marker))
	true_phenotypes = apply(genotypes_normalized, 1, function(x) betas %*% x) 
	phenotypes = true_phenotypes + rnorm(n_indiv, sd = sqrt(1-hsq))
	l = list()
	l$phenotypes = phenotypes
	l$true_phenotypes = true_phenotypes
	return(l)
}



