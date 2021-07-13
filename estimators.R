library(expm)
library(reticulate)
library(qgg)

l2normsq = function(x) {
	return(sum(x ^ 2))	
}

dicker_2014_estimator_2 = function(X, y) {
	d = ncol(X) #Number of markers
	n = nrow(X) #Number of individuals
	if(n != length(y)) {
		stop("Error in dicker2014_version2: phenotype length doesn't match genotype dim")
	}
	#Center X
	gw = apply(X, 2, function(x) return((x - mean(x))))
	
	y = y - mean(y)
	
	xtx = t(gw) %*% gw
	m1 = 1/d * sum(diag(1/n * xtx))
	m2 = 1/d * sum(diag( t(1/n * xtx) %*% (1/n * xtx))) -
		1/(d*n) * (sum(diag(1/n * xtx)) ^2 )
	sigmatilde2 = (1 + d * m1^2 / ((n+1) * m2)) * 1/n * sum(y^2) -
		m1/(n*(n+1)*m2)*sum((t(gw) %*% y)^2)
	tautilde2 = -d * m1^2/(n*(n+1)*m2) * sum(y^2) +
		m1/(n*(n+1)*m2) * sum((t(gw) %*% y)^2)
	hsq_est = tautilde2 / (sigmatilde2 + tautilde2)
	return(c(hsq_est, tautilde2, sigmatilde2))
}


dicker_2014_estimator_2_sd = function(X, y) {
	d = ncol(X) #Number of markers
	n = nrow(X) #Number of individuals
	if(n != length(y)) {
		stop("Error in dicker2014_version2: phenotype length doesn't match genotype dim")
	}
	#Center X
	gw = apply(X, 2, function(x) return((x - mean(x))/sd(x)))
	
	y = y - mean(y)
	
	xtx = t(gw) %*% gw
	m1 = 1/d * sum(diag(1/n * xtx))
	m2 = 1/d * sum(diag( t(1/n * xtx) %*% (1/n * xtx))) -
		1/(d*n) * (sum(diag(1/n * xtx)) ^2 )
	sigmatilde2 = (1 + d * m1^2 / ((n+1) * m2)) * 1/n * sum(y^2) -
		m1/(n*(n+1)*m2)*sum((t(gw) %*% y)^2)
	tautilde2 = -d * m1^2/(n*(n+1)*m2) * sum(y^2) +
		m1/(n*(n+1)*m2) * sum((t(gw) %*% y)^2)
	hsq_est = tautilde2 / (sigmatilde2 + tautilde2)
	return(c(hsq_est, tautilde2, sigmatilde2))
}

dicker_2014_estimator_1 = function(genotypes, phenotypes) {
	#Make some values
	n = nrow(genotypes)
	d = ncol(genotypes)
	
	#Center X
	gw = apply(genotypes, 2, function(x) return((x - mean(x))))
	
	#Center y
	phenotypes = phenotypes - mean(phenotypes)
	
	#Assume that the LD matrix is the identity
	Shalf = diag(ncol(gw))
	
	xty = t(gw) %*% phenotypes
	
	if(rankMatrix(Shalf) == nrow(Shalf)) {
		sigma = (d + n + 1)/(n * (n + 1)) * sum(phenotypes^2) - 1/(n*(n+1)) * sum((xty)^2)
		tau = -d/(n*(n+1)) * sum(phenotypes^2) + 1/(n*(n+1)) * sum((xty)^2)
		hsq_est = tau/(sigma+tau)
	} else {
		hsq_est= -1
		sigma = -1
		tau = -1
	}
	
	return(c(hsq_est, tau, sigma))
}

dicker_2014_estimator_1_sd = function(genotypes, phenotypes) {
	#Make some values
	n = nrow(genotypes)
	d = ncol(genotypes)
	
	#Center X
	gw = apply(genotypes, 2, function(x) return((x - mean(x))/sd(x)))
	
	#Center y
	phenotypes = phenotypes - mean(phenotypes)
	
	#Assume that the LD matrix is the identity
	Shalf = diag(ncol(gw))
	
	xty = t(gw) %*% phenotypes
	
	if(rankMatrix(Shalf) == nrow(Shalf)) {
		sigma = (d + n + 1)/(n * (n + 1)) * sum(phenotypes^2) - 1/(n*(n+1)) * sum((xty)^2)
		tau = -d/(n*(n+1)) * sum(phenotypes^2) + 1/(n*(n+1)) * sum((xty)^2)
		hsq_est = tau/(sigma+tau)
	} else {
		hsq_est= -1
		sigma = -1
		tau = -1
	}
	
	return(c(hsq_est, tau, sigma))
}

schwartzman_2019_estimator = function(genotypes, phenotypes) {
	n_marker = ncol(genotypes)
	n_indiv = nrow(genotypes)
	
	#Center X
	gw = apply(genotypes, 2, function(x) return(x - mean(x)))
	
	#Renormalize columns as specified in paper
	gw = sqrt(n_indiv - 1) * apply(gw, 2, function(x) return(x / sqrt(l2normsq(x))))
	
	#Center y
	phenotypes = phenotypes - mean(phenotypes)
	
	#Renormalize y as sepcified in paper
	phenotypes = sqrt(n_indiv - 1) * phenotypes / sqrt(l2normsq(phenotypes))
	
	#Compute the correlation scores
	u = t(gw) %*% phenotypes / sqrt(n_indiv - 1)
	
	#Compute variance of correlation scores
	s2 = l2normsq(u) / n_marker
	
	#Compute the matrix
	S = t(gw) %*% gw / (n_indiv - 1)
	S2 = t(S) %*% S
	mu2 = 1/n_marker * sum(diag(S2)) - (n_marker - 1)/(n_indiv - 1)
	
	#Calculate heritability estimate
	return(hsq_est = n_marker / n_indiv / mu2 * (s2 - 1))
}

#Returns (hsq, sigma_g^2, sigma_e^2)
elizabeth_MoM_no_diag = function(genotypes, phenotypes) {
	n_marker = ncol(genotypes)
	n_indiv = nrow(genotypes)
	
	#Whole genotypes, including duplication
	#Normalized whole genotypes
	gwn = normalize_genotypes(genotypes)
	G = 1/n_marker * gwn %*% t(gwn)
	
	#Calculate MoM estimates
	numerator = 0
	denom = 0
	for(i in 2:n_indiv) {
		for(j in 1:(i-1)) {
			numerator = numerator + phenotypes[i] * phenotypes[j] *G[i,j]
			denom = denom + G[i,j]^2
		}
	}
	sigma_g = numerator/denom
	sigma_y = var(phenotypes)
	hsq_est = sigma_g/sigma_y
	return(c(hsq_est, sigma_g, sigma_y - sigma_g))
}

#Returns (hsq, sigma_g^2, sigma_e^2)
elizabeth_MoM_diag = function(genotypes, phenotypes) {
	n_marker = ncol(genotypes)
	n_indiv = nrow(genotypes)
	
	#Normalized whole genotypes
	gwn = normalize_genotypes(genotypes)
	G = 1/n_marker * gwn %*% t(gwn)
	
	#Calculate MoM estimates
	numerator = 0
	denom = 0
	for(i in 2:n_indiv) {
		for(j in 1:(i-1)) {
			numerator = numerator + phenotypes[i] * phenotypes[j] *G[i,j]
			denom = denom + G[i,j]^2
		}
	}
	
	#Calculate estimates using diagonal for sigma_e
	numerator_e = 0
	denominator_e = 0
	sigma_g = numerator/denom
	sigma_y = var(phenotypes)
	sigma_e = (sum(phenotypes^2) - (sigma_g * sum(diag(G))))/n_indiv
	hsq_est = sigma_g/(sigma_e + sigma_g)
	return(c(hsq_est, sigma_g, sigma_e))
}

#Returns (hsq, sigma_g^2, sigma_e^2)
max_like_R = function(genotypes, phenotypes) {
	source_python("~/Documents/UW Documents/gcta_proj/Scripts/20_06_24_Estimators/max_like.py")
	res = max_like(g$genotypes, p$phenotypes)
	#The python code returns the estimate for sigma_g^2 in res[1] and 
	#the estimate for sigma_e^2 in res[2]
	hsq_est = res[1] / (res[1] + res[2])
	return(c(hsq_est, res[1], res[2]))
}

max_like_greml = function(genotypes, phenotypes) {
	genotypes = normalize_genotypes(genotypes)
	grm = genotypes %*% t(genotypes) / ncol(genotypes)
	g = greml(y = phenotypes, X = rep(1, nrow(genotypes)), GRM = list(grm))
	return(c(g$theta[1]/(g$theta[2] + g$theta[1]), g$theta[1], g$theta[2]))
}

ldak = function(filepath, power = -.25) {
	cur_wd = getwd()
	setwd(filepath)
	system("~/LDAK/ldak5.1.linux --bfile test --cut-weights sections")
	system("~/LDAK/ldak5.1.linux --bfile test --calc-weights-all sections")
	system(paste("~/LDAK/ldak5.1.linux --calc-kins-direct kins --bfile test --weights sections/weights.all --power", power))
	system("~/LDAK/ldak5.1.linux --reml quant --grm kins --pheno test.phen")
	
	#Sorry this is a weird hacky way to get the heritability
	herit = as.numeric(strsplit(readLines("quant.reml")[20], split = " ")[[1]][2])
	setwd(cur_wd)
	
	return(herit)
}

gcta = function(filepath) {
	cur_wd = getwd()
	setwd(filepath)
	system("gcta64 --bfile test --make-grm --out test")
	system("gcta64 --grm test --pheno test.phen --reml --out gcta_herit")
	
	#Sorry this is a weird hacky way to get the heritability
	herit = as.numeric(unlist(strsplit(readLines("gcta_herit.hsq")[5], "\t"))[2])
	ve = as.numeric(unlist(strsplit(readLines("gcta_herit.hsq")[3], "\t"))[2])
	vg = as.numeric(unlist(strsplit(readLines("gcta_herit.hsq")[2], "\t"))[2])
	setwd(cur_wd)
	
	return(c(herit, vg, ve))
}

