source("gcta.R")
source("genotypes_cor.R")
source("estimators.R")
source("genotypes_related.R")
library(rres)

library(mvtnorm)
#set.seed(23)
sampled_markers_subpop <- read.csv("markers1shuffle.txt", sep="", stringsAsFactors=FALSE)

n_iter = 1000
repeat_proportion = .1
hsq = .8
clean_up = T

n_vec = rep(0, 80000)
p_vec = rep(0, 80000)
iteration_vec = rep(0, 80000)
structure_vec = rep(0, 80000)
param = 0
counter = 1
for(p in seq(400, 8000, by = 400)) {
    n = 400
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

array_index = as.integer(commandArgs(trailingOnly = T))[1]
indexes = 1:100 + (array_index - 1) * 100
for(i in indexes) {
	n = settings$n[i]
	p = settings$p[i]
	structure = settings$structure[i]
	iteration = settings$iteration[i]

	if (structure == "cousin0") {
	    g = get_genotypes(n, p, 0, 0, mafs = sampled_markers_subpop$afr_alt[1:p])$genotypes
	    pheno = get_phenotypes(g, hsq)
	} else if(structure == "cousin1") {
	    g = get_related_genotypes_1st_cous(n, p, sampled_markers_subpop$afr_alt[1:p])
	    pheno = get_phenotypes(g, hsq)
	} else if (structure == "cousin2") {
	    g = get_related_genotypes_2nd_cous(n, p, sampled_markers_subpop$afr_alt[1:p])
	    pheno = get_phenotypes(g, hsq)
	} else if (structure == "cousin3") {
	    g = get_related_genotypes_3rd_cous(n, p, sampled_markers_subpop$afr_alt[1:p])   
	    pheno = get_phenotypes(g, hsq)
	}
	
	if(!is.null(g)) {
		pheno_with_noise = pheno$phenotypes
		
		folder = get_folder(n, p, structure, param, iteration)
		system(paste("mkdir -p", folder))
		
		save(g, file = get_file_name(n, p, structure, param, iteration, "genotypes.Rda"))
		save(pheno, file = get_file_name(n, p, structure, param, iteration, "phenotypes.Rda"))
		
		map_file = get_samp_markers(ncol(g))
		write.table(map_file, file = get_file_name(n, p, structure, param, iteration, "test.map"), row.names = F, col.names = F, quote = F)
		
		ped_file = get_ped_file(g, p, map_file)
		write.table(ped_file, file = get_file_name(n, p, structure, param, iteration, "test.ped"), row.names = F, col.names = F, quote = F)
		
		phen_file = get_phen_file(pheno_with_noise)
		write.table(phen_file, file = get_file_name(n, p, structure, param, iteration, "test.phen"), row.names = F, col.names = F, quote = F)
		system(paste("plink --file", get_file_name(n, p, structure, param, iteration, "test"), "--make-bed --out", get_file_name(n, p, structure, param, iteration, "test")))
		print(paste("Just finished ", get_file_name(n, p, structure, param, iteration, "test")))
		
		
		###############################
		# DO ANALYSES #################
		###############################
		if(structure == "haplotypes") {
			next
		}
		
		folder = get_folder(n, p, structure, param, iteration)
		load(get_file_name(n, p, structure, param, iteration, "genotypes.Rda"))
		load(get_file_name(n, p, structure, param, iteration, "phenotypes.Rda"))
		
		gold_standard = var(pheno$true_phenotypes)/var(pheno$phenotypes)
		dicker2 = dicker_2014_estimator_2_sd(g, pheno$phenotypes)
		dicker1 = Re(dicker_2014_estimator_1_sd(g, pheno$phenotypes))
		schwartzman = schwartzman_2019_estimator(g, pheno$phenotypes)
		elizabeth_no_diag = elizabeth_MoM_no_diag(g, pheno$phenotypes)
		elizabeth_diag = elizabeth_MoM_diag(g, pheno$phenotypes)
		gcta_maxlike = gcta(folder)
		ldak_res = ldak(folder)
		
		names = c("gold_standard", "dicker2", "dicker1", "schwartzman", "elizabeth_no_diag", "elizabeth_diag", "gcta", "ldak")
		result_hsq = c(gold_standard[1], dicker2[1], dicker1[1], schwartzman[1], elizabeth_no_diag[1], elizabeth_diag[1], gcta_maxlike[1], ldak_res[1])
		sigmag = c(gold_standard[2], dicker2[2], dicker1[2], schwartzman[2], elizabeth_no_diag[2], elizabeth_diag[2], gcta_maxlike[2], ldak_res[2])
		sigmae = c(gold_standard[3], dicker2[3], dicker1[3], schwartzman[3], elizabeth_no_diag[3], elizabeth_diag[3], gcta_maxlike[3], ldak_res[3])
		df = data.frame(names, result_hsq, sigmae, sigmag, n, p, structure, param, iteration)
		if(clean_up) {
			system(paste("/bin/rm -r", get_file_name(n, p, structure, param, iteration, "*")))
		}
		save(df,file = get_file_name(n, p, structure, param, iteration, "alan_results.Rda")) 
	}
}


