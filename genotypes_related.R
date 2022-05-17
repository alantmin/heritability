library(mvtnorm)

get_related_genotypes_3rd_cous = function(n_indiv, n_marker, mafs, seglength = 3000) {
	n_founder = 242
	if(n_indiv %% 400 != 0) {
		error("n_indiv must be divisible by 400")
		return(0)
	}
	if(length(mafs) != n_marker) {
		error("Length mafs != n_marker")
		return(0)
	}
	genotype_list = list()
	
	for(i in 1:(n_indiv/400)) {
		#Create a pedigree with 80 3rd cousins, but it actually only takes the first 40
	    #I had written code to make 80 3rd cousins and I was too lazy to change it to 40
	    tmp_genotype_list = list()
	    for(j in 1:5) {
    	    pedigree = as.character(rep(1, 562))
    	    member = as.character(c(1001:1002, 2001:2160, 3001:3160, 4001:4160, 5001:5080))
    	    sex = c(1, 2, rep(1, 80), rep(2, 80), rep(1, 80), rep(2, 80), rep(1, 80), rep(2, 80), rep(1, 80))
    	    father = as.character(c(NA, NA, rep(1001, 80), rep(NA, 80), 2001:2080, rep(NA, 80), 3001:3080, rep(NA, 80), 4001:4080))
    	    mother = as.character(c(NA, NA, rep(1002, 80), rep(NA, 80), 2081:2160, rep(NA, 80), 3081:3160, rep(NA, 80), 4081:4160))
    		pedinfo = data.frame(pedigree, member, sex, father, mother, stringsAsFactors = FALSE)
    		
    		#Simulate recombination maybe?
    		inher = sim.recomb(pedinfo, seglength)
    		
			g_list = get_genotypes_haplo(n_founder * 2, n_marker, 0, 0, mafs)	
			haplotype = cbind(g_list$genotypes, g_list$genotypes_rep)
			#Make the markers evenly spaced out
			marker = seq(0, 1, length.out = n_marker)
    	
    		
    		#Simulate the snps
    		genotype = populate.snp(inher, haplotype, marker, c(483:562), output.allele = FALSE)
    		tmp_genotype_list[[j]] = genotype
	    }

	    genotype_list[[i]] = do.call(rbind, tmp_genotype_list)
	}
	return(do.call(rbind, genotype_list))
}


get_related_genotypes_2nd_cous = function(n_indiv, n_marker, mafs, seglength = 3000) {
    n_founder = 162
    if(n_indiv %% 400 != 0) {
        error("n_indiv must be divisible by 400")
        return(0)
    }
    if(length(mafs) != n_marker) {
        error("Length mafs != n_marker")
        return(0)
    }
    genotype_list = list()
    
    for(i in 1:(n_indiv/400)) {
        #Create a pedigree with 80 2nd cousins but it actually only takes the first 40
        #I had written code to make 80 2nd cousins and I was too lazy to change it to 40
        tmp_genotype_list = list()
        for(j in 1:5) {
            pedigree = as.character(rep(1, 402))
            member = as.character(c(1001:1002, 2001:2160, 3001:3160, 4001:4080))
            sex = c(1, 2, rep(1, 80), rep(2, 80), rep(1, 80), rep(2, 80), rep(1, 80))
            father = as.character(c(NA, NA, rep(1001, 80), rep(NA, 80), 2001:2080, rep(NA, 80), 3001:3080))
            mother = as.character(c(NA, NA, rep(1002, 80), rep(NA, 80), 2081:2160, rep(NA, 80), 3081:3160))
            pedinfo = data.frame(pedigree, member, sex, father, mother, stringsAsFactors = FALSE)
            
            #Simulate recombination maybe?
            inher = sim.recomb(pedinfo, seglength)
            
            g_list = get_genotypes_haplo(n_founder * 2, n_marker, 0, 0, mafs)	
            haplotype = cbind(g_list$genotypes, g_list$genotypes_rep)
            #Make the markers evenly spaced out
            marker = seq(0, 1, length.out = n_marker)
            
            
            #Simulate the snps
            genotype = populate.snp(inher, haplotype, marker, c(323:402), output.allele = FALSE)
            tmp_genotype_list[[j]] = genotype
        }
        
        genotype_list[[i]] = do.call(rbind, tmp_genotype_list)
    }
    return(do.call(rbind, genotype_list))
}

get_related_genotypes_1st_cous = function(n_indiv, n_marker, mafs, seglength = 3000) {
    n_founder = 82
    if(n_indiv %% 400 != 0) {
        error("n_indiv must be divisible by 400")
        return(0)
    }
    if(length(mafs) != n_marker) {
        error("Length mafs != n_marker")
        return(0)
    }
    genotype_list = list()
    
    for(i in 1:(n_indiv/400)) {
        #Create a pedigree with 80 2nd cousins but it actually only takes the first 40
        #I had written code to make 80 2nd cousins and I was too lazy to change it to 40
        tmp_genotype_list = list()
        for(j in 1:5) {
            pedigree = as.character(rep(1, 242))
            member = as.character(c(1001:1002, 2001:2160, 3001:3080))
            sex = c(1, 2, rep(1, 80), rep(2, 80), rep(1, 80))
            father = as.character(c(NA, NA, rep(1001, 80), rep(NA, 80), 2001:2080))
            mother = as.character(c(NA, NA, rep(1002, 80), rep(NA, 80), 2081:2160))
            pedinfo = data.frame(pedigree, member, sex, father, mother, stringsAsFactors = FALSE)
            
            #Simulate recombination maybe?
            inher = sim.recomb(pedinfo, seglength)
            
            g_list = get_genotypes_haplo(n_founder * 2, n_marker, 0, 0, mafs)	
            haplotype = cbind(g_list$genotypes, g_list$genotypes_rep)
            #Make the markers evenly spaced out
            marker = seq(0, 1, length.out = n_marker)
            
            
            #Simulate the snps
            genotype = populate.snp(inher, haplotype, marker, c(163:242), output.allele = FALSE)
            tmp_genotype_list[[j]] = genotype
        }
        
        genotype_list[[i]] = do.call(rbind, tmp_genotype_list)
    }
    return(do.call(rbind, genotype_list))
}


#Function to generate genotypes
#n_indiv is how many individuals to generate
#n_markers is how many markers to generate
#by default generates markers with frequency p=.2. 
#repeat_subset: how many markers are being repeated
#repeat_num: how many times repeat_subset is being repeated
#Returns a list with $genotypes being the original genotypes and $genotypes_rep being the repeated genotypes. 
get_genotypes_haplo = function(n_indiv, n_markers, repeat_subset, repeat_num, mafs = NULL) {
	if(is.null(mafs)) {
		error("YOU HAVE TO SPECIFY MAFS")
		return(0)
	} else {
		if(length(mafs) != n_markers) {
			print("Length of mafs not equal to n_markers!")
			return(0)
		}
		genotypes = matrix(rbinom(n_indiv * n_markers, size = 1, prob = rep(mafs, n_indiv)), nrow = n_indiv, ncol = n_markers, byrow = T) + 1 #Makes each column have the right maf
	}
	genotypes_rep = matrix(rep(genotypes[,1:repeat_subset], repeat_num), nrow = n_indiv)
	l = list()
	l$genotypes = genotypes
	l$genotypes_rep = genotypes_rep
	return(l)
}