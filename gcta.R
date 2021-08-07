# Function that creates a MAP file for PLINK
# @param chr: Either an integer or a vector of the chromosome
#	By default, every marker is from the chromosome 1
# @param rsnumber: the index of the marker
# @param genetic_distance: genetic distance. By default,
# 	this is set to be uniform from 0 to 1
# @param b_p_pos: Base pair position is also just the index
# @param ref_allele: Reference nucleotides at each position.  
#	By default, reference allele is always A.
# @param alt_allele: Alternate allele at each position.
#	By default, the alternate allele is always T
get_samp_markers = function(n_marker, chr = NULL, rs_number = NULL, genetic_distance = NULL, bp_pos = NULL, ref_allele = NULL, alt_allele = NULL) {
	#Set parameters to defaults if they are not specified 
	if(is.null(chr)) {
		chr = rep(1, n_marker)	
	}
	if(is.null(rs_number)) {
		rs_number = paste("rs", 1:n_marker, sep = "")	
	}
	if(is.null(genetic_distance)) {
		genetic_distance = seq(0, 1, length.out = n_marker)	
	}
	if(is.null(bp_pos)) {
		bp_pos = 1:n_marker
	}
	if(is.null(ref_allele)) {
		ref_allele = rep('A', n_marker)
	}
	if(is.null(alt_allele)) {
		alt_allele = rep('T', n_marker)	
	}
	
	#Check that parameters are the right length
	if(length(chr) != n_marker || length(rs_number) != n_marker || length(genetic_distance) != n_marker || length(bp_pos) != n_marker || length(ref_allele) != n_marker || length(alt_allele) != n_marker) {
		stop("get_samp_markers: The length of the parameters has to be the same as the number of markers")	
	}
	
	return(data.frame(chr, rs_number, genetic_distance, bp_pos, ref_allele, alt_allele))
}

#Takes for example (0,0,1,2), "A", "T" and returns "T T  T T  A T  A A" 
number_to_letter = function(genotype_row, ref_allele, alt_allele) {
	mat = matrix(alt_allele, nrow = length(genotype_row), ncol = 2)
	mat[genotype_row == 1, 1] = ref_allele
	mat[genotype_row == 2, c(1,2)] = ref_allele
	return(mat)
}

# Function that creates the ped file for GCTA to use.
# @param genotypes: a genotype matrix with n_indiv rows
#	and n_marker columns.
# @param phenotypes: a vector with n_indiv entries
#	corresponding the the phenotypes of each of the individuals
# @param samp_markers: the MAP file format from PLINK. 
# @param fam_id: Id numbers for the families. Each individual
#	is assumed to be in their own family unless otherwise specified
# @param letter_genotype_fn: A file name for the letter genotypes.
#	If NULL, just make a new one. Letter genotypes means it's
#	in the format A A T T, etc.
# @param save_letter_genotype: A file name to save the letter 
#	genotypes to. If NULL, do not save. This option only matters
#	if letter_genotype_fn is NULL
# @param sex: Sex of each individual. Every individual is assumed
#	to be male unless otherwise specified
get_ped_file = function(genotypes, phenotypes, samp_markers, fam_id = NULL, indiv_id = NULL, pat_id = NULL, mat_id = NULL, sex = NULL, letter_genotype_fn = NULL, save_letter_genotype = NULL) {
	#Count number of individuals and markers
	n_indiv = nrow(genotypes)
	n_marker = ncol(genotypes)
	
	#Set parameters to defaults if they are NULL
	if(is.null(fam_id)) {
		fam_id = paste("FAM", 1:n_indiv, sep = "")
	}
	if(is.null(indiv_id)) {
		indiv_id = 1:n_indiv
	}
	if(is.null(pat_id)) {
		pat_id = rep(0, n_indiv)
	}
	if(is.null(mat_id)) {
		mat_id = rep(0, n_indiv)
	}
	if(is.null(sex)) {
		sex = rep(1, n_indiv)
	}
	
	#Check that dimensions are correct
	if(length(fam_id) != n_indiv || length(indiv_id) != n_indiv || length(pat_id) != n_indiv || length(mat_id) != n_indiv || length(sex) != n_indiv) {
		stop("get_ped_file: Check dimensions of input parameters")
	}
	
	if(nrow(samp_markers) != n_marker) {
		stop("get_ped_file: Check dimensions of samp_markers")
	}
	
	if(is.null(letter_genotype_fn)) {
		#Convert the genotypes from numbers into letters
		letter_genotype = matrix("Z", nrow = n_indiv, ncol = n_marker * 2)
		ref_allele = as.character(samp_markers$ref_allele)
		alt_allele = as.character(samp_markers$alt_allele)
		for(i in 1:n_marker) {
			letter_genotype[,(i*2-1):(i*2)] = number_to_letter(genotypes[,i], ref_allele[i], alt_allele[i])
		}
		if(!is.null(save_letter_genotype)) {
			save(letter_genotype, file = save_letter_genotype)
		}
	} else {
		load(letter_genotype_fn)
	}
	return(cbind(data.frame(fam_id, indiv_id, pat_id, mat_id, sex, phenotypes), letter_genotype))
}

get_phen_file = function(phenotypes) {
	n_indiv = length(phenotypes)
	phen_df = data.frame(paste("FAM", 1:n_indiv, sep = ""), 1:n_indiv, phenotypes)
	return(phen_df)	
}

# Function to run the GCTA program, given a genotype matrix. 
# @param genotypes: Genotype matrix should have each row as
#	an individual and each column as a marker. Should consist of 
#	integer values 0, 1, or 2. 
# @param samp_markers: a matrix with information about each marker
# @param new_ped: if true, creates a new ped from scratch.
#	If false, use the letter genotypes from a previous run.
# @param plink_fn: a file name to write the samp_markers (MAP file)
#	and the ped file to (should not have a file extension)
#new_samp_markers: if true, creates a new sample marker file.
#	Otherwise doesn't do anything with the samp markers argument.
run_gcta = function(genotypes, phenotypes, samp_markers, plink_fn, new_ped = T) {
	ped_fn = paste(plink_fn, ".ped", sep = "")
	samp_markers_fn = paste(plink_fn, ".map", sep = "")
	phen_fn = paste(plink_fn, ".phen", sep = "")
	n_indiv = nrow(genotypes)
	n_marker = ncol(genotypes)
	
	#Get the ped file, and only make a new one if new_ped is true
	if(new_ped) {
		ped_file = get_ped_file(genotypes, phenotypes, samp_markers, save_letter_genotype = ped_fn)
	} else {
		ped_file = ged_ped_file(genotypes, phenotypes, samp_markers, letter_genotype_fn = ped_fn)
	}
	write.table(ped_file, file = ped_fn, quote = F, col.names = F, row.names = F)
	
	write.table(samp_markers, file = samp_markers_fn, quote = F, col.names =F, row.names = F)
	
	#Runs plink to generate some test files
	system(paste("plink -file ", plink_fn , " -out ", plink_fn, sep = ""))
	
	#Run GCTA to get the GRM files
	system(paste("gcta64 --bfile ", plink_fn," --make-grm --out test", sep = ""))
	
	#Write a .phen file for GCTA
	phen_df = data.frame(paste("FAM", 1:n_indiv, sep = ""), 1:n_indiv, phenotypes)
	write.table(phen_df, file = phen_fn, quote = F, col.names = F, row.names = F)
	
	#Run GCTA to do GREML analysis
	system(paste("gcta64 --grm ", plink_fn, " --pheno ", phen_fn, " --reml --out test"))
}

# Function that cleans up a folder 
# @param file_name: cleans all files in this folder
clean_dir = function(file_name) {
	rm_list = paste(file_name, c(".bed", ".bim", ".fam", ".grm.bin", ".grm.id", ".grim.N.bin", ".hsq", ".log", ".map", ".ped", ".phen"), sep = "")
	for (r in rm_list) {
		system(paste("rm", r))	
	}
}

# Function that generates phenotypes with a particular heritability.
# @param genotypes: non-normalized genotype matrix.
# @param hsq: true heritability that you want to simulate
# @value: returns a list with phenotypes and true phenotypes
#	phenotypes include the error term from sigma_e^2
#	true phenotypes do not have sigma_e^2
get_phenotypes = function(genotypes, hsq) {
	n_indiv = nrow(genotypes)
	n_marker = ncol(genotypes)
	
	genotypes_normalized = genotypes
	for (i in 1:ncol(genotypes_normalized)) {
		p = sum(genotypes_normalized[,i])/(2 * n_indiv)
		#if monoallelic, give an error
		if(p == 0 || p == 1) {
			stop("get_phenotypes: Careful! Monoallelic marker")
		}
		genotypes_normalized[,i] = (genotypes[,i] - 2*p)/sqrt(2*p*(1-p))
	}
	
	betas = rnorm(n_marker, mean = 0, sd = sqrt(hsq/n_marker))
	true_phenotypes = apply(genotypes_normalized, 1, function(x) betas %*% x) 
	phenotypes = true_phenotypes + rnorm(n_indiv, sd = sqrt(1-hsq))
	l = list()
	l$phenotypes = phenotypes
	l$true_phenotypes = true_phenotypes
	return(l)
}

# Function that normalizes the genotypes
# @param genotypes: a matrix with genotypes to be normalized
# @value: returns a normalized genotype matrix
normalize_genotypes = function(genotypes) {
	n_indiv = nrow(genotypes)
	n_marker = ncol(genotypes)
	genotypes_normalized = genotypes
	for (i in 1:ncol(genotypes_normalized)) {
		p = sum(genotypes_normalized[,i])/(2 * n_indiv)
		#if monoallelic, give an error
		if(p == 0 || p == 1) {
			stop("get_phenotypes: Careful! Monoallelic marker")
		}
		genotypes_normalized[,i] = (genotypes[,i] - 2*p)/sqrt(2*p*(1-p))
	}
	return(genotypes_normalized)
}

# Function to generate genotypes
# @param n_indiv: how many individuals to generate
# @param n_markers: how many markers to generate
# @param repeat_subset: how many markers are being repeated
# @param repeat_num: how many times repeat_subset is being repeated
# @param mafs: a vector or number of the minor allele frequency
#	by default generates markers with frequency p=.2. 
# @value: Returns a list with $genotypes being the original
#	genotypes and $genotypes_rep being the repeated genotypes. 
get_genotypes = function(n_indiv, n_markers, repeat_subset, repeat_num, mafs = NULL) {
	if(is.null(mafs)) {
		genotypes = matrix(rbinom(n_indiv * n_markers, size = 2, prob = .2), nrow = n_indiv, ncol = n_markers)
	} else {
		if(length(mafs) != n_markers) {
			print("Length of mafs not equal to n_markers!")
			return(0)
		}
		rafs = 1 - mafs #Reference allele frequency is what we actually have here
		genotypes = matrix(rbinom(n_indiv * n_markers, size = 2, prob = rep(rafs, n_indiv)), nrow = n_indiv, ncol = n_markers, byrow = T) #Makes each column have the right maf
	}
	genotypes_rep = matrix(rep(genotypes[,1:repeat_subset], repeat_num), nrow = n_indiv)
	l = list()
	l$genotypes = genotypes
	l$genotypes_rep = genotypes_rep
	return(l)
}


#Get the hsq value from a file
#Returns VE,VG,VG/VP
get_hsq = function(fn) {
	gcta_output = readLines(fn)
	hsq = as.numeric(strsplit(gcta_output[5], "\t")[[1]][2])
	vg = as.numeric(strsplit(gcta_output[2], "\t")[[1]][2])
	ve = as.numeric(strsplit(gcta_output[3], "\t")[[1]][2])
	return(c(hsq, vg, ve))
}

get_folder = function(n, p, structure, param, iteration) {
	return(paste("data/n", n, "p", p, "/", structure, "/param", param, "/iteration", sprintf("%03d", iteration), sep = ""))
}

get_file_name = function(n, p, structure, param, iteration, extension) {
	return(paste("data/n", n, "p", p, "/", structure, "/param", param, "/iteration", sprintf("%03d", iteration), "/", extension, sep = ""))
}











