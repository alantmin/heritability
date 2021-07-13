#This one makes 80% of the markers at a lower LD with the "noncausal" markers
#And 20% of the markers are at a higher LD with the "noncausal markers"

source("gcta.R")
source("genotypes_cor.R")
source("estimators.R")

library(mvtnorm)
#set.seed(23)
sampled_markers_subpop <- read.csv("markers1shuffle.txt", sep="", stringsAsFactors=FALSE)

n_iter = 500
repeat_proportion = .1
hsq = .8
clean_up = T

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

array_index = as.integer(commandArgs(trailingOnly = T))[1]
indexes = 1:50 + (array_index - 1) * 50 
for(i in indexes) {
	n = settings$n[i]
	p = settings$p[i]
	structure = settings$structure[i]
	param = settings$param[i]
	iteration = settings$iteration[i]
	
	if(structure == "autocorrelation") {
		g = get_genotypes_autocor(n, p, param, T, sampled_markers_subpop$afr_alt[1:p])
		pheno = get_phenotypes(g[,seq(1,p,by = 2)], hsq)
	} else if (structure == "repeat") {
		g_list = get_genotypes(n, p, round(p * repeat_proportion), param, sampled_markers_subpop$afr_alt[1:p])
		pheno = get_phenotypes(g_list$genotypes[,seq(1,p,by = 2)], hsq)
		g = cbind(g_list$genotypes, g_list$genotypes_rep)
	} else if (structure == "block") {
		g = get_genotypes_block(n, p, param, T, sampled_markers_subpop$afr_alt[1:p])
		pheno = get_phenotypes(g[,seq(1,p,by = 2)], hsq)
	} else {
		print("There is no such structure option!!!")
		quit()
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
		ldak_res_1 = ldak(folder, power=-.25)
		ldak_res_2 = ldak(folder, power=-1)
		
		names = c("gold_standard", "dicker2", "dicker1", "schwartzman", "elizabeth_no_diag", "elizabeth_diag", "gcta", "ldak1", "ldak2")
		result_hsq = c(gold_standard[1], dicker2[1], dicker1[1], schwartzman[1], elizabeth_no_diag[1], elizabeth_diag[1], gcta_maxlike[1], ldak_res_1[1], ldak_res_2[1])
		sigmag = c(gold_standard[2], dicker2[2], dicker1[2], schwartzman[2], elizabeth_no_diag[2], elizabeth_diag[2], gcta_maxlike[2], ldak_res_1[2], ldak_res_2[2])
		sigmae = c(gold_standard[3], dicker2[3], dicker1[3], schwartzman[3], elizabeth_no_diag[3], elizabeth_diag[3], gcta_maxlike[3], ldak_res_1[3], ldak_res_2[3])
		df = data.frame(names, result_hsq, sigmae, sigmag, n, p, structure, param, iteration)
		if(clean_up) {
			system(paste("rm", get_file_name(n, p, structure, param, iteration, "*")))
		}
		save(df,file = get_file_name(n, p, structure, param, iteration, "alan_results.Rda")) 
	}
}


