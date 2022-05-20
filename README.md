# Introduction

These are files used in the publication of

Min, Alan, Elizabeth A. Thompson, and Saonli Basu. "Comparing Heritability Estimators under Alternative Structures of Linkage Disequilibrium." bioRxiv (2021)

which is available at https://www.biorxiv.org/content/10.1101/2021.09.08.459523v1. 

The code assumes that `GCTA` (https://yanglab.westlake.edu.cn/software/gcta/), `LDAK` (https://dougspeed.com/ldak/), and `plink` are installed (https://zzz.bwh.harvard.edu/plink/).

The files are described below

1. `estimators.R`: The estimators are implemented in this file 
2. `gcta.R`: This file includes helper functions to run GCTA
3. `genotypes_cor.R`: This file provides functions to simulate genotypes with correlation structures. 
4. `genotypes_related.R`: This file provide functions to simulate genotypes with different levels of cousinship.
5. `simulation1.R`: This file uses the simulation functions to simulate genotypes and phenotypes using different correlation structures, then outputs a folder with results of the estimated values of heritability for various estimators.
12. `simulation2.ipynb`: This is a python notebook with code used for generating the likelihood plots.
6. `simulation3.R`: This file uses the simulation functions to simulate genotypes and phenotypes using different levels of relatedness, then outputs a folder with results of the estimated values of heritability for various estimators. 
8. `simulation1_collate_results.R`: This file looks through all of the folders that are output by the simulation1 script and puts the results into one data frame. 
9. `simulation3_collate_results.R`: This file looks through all of the folders that are output by the simulation2 script and puts the results into one data frame for the cousins data. 
10. `plotresults.Rmd`: This is a R Markdown file that makes plots given the data frame created by `collate_results.R`
11. `plotresults_cousins.Rmd`: This is a R Markdown file that makes plots given the data frame created by `collate_results_cousins.R`
13. `simulation1_results.Rda`: This is a R data file with the output of `simulation1.R`
14. `simulation3_results.Rda`: This is a R data file with the output of `simulation3.R`
7. `markers1shuffle.txt`: This file provides SNPs from 1000 genomes in a random order that are filtered so that the minor allele frequency is greater than 0.05.

# Note LDAK, GCTA, and plink
The `estimators.R` file requires that LDAK, GCTA, and plink can be run in the command line using 

```
~/LDAK/ldak5.1.linux
gcta64
plink
```
respectively.

# Usage 

The scripts `simulation1.R` and `simulation3.R` first create a data frame with a list of parameters to test. They each take as input an integer `i` which determines which subset of the parameters are simulated. Results of the simulations are then stored in a folder structure in the same directory. 

For simulation 1, the sequence of commands to reproduce the figures is as follows

```
for i in {1..600}; do
	Rscript simulation1.R $i
done 

Rscript simulation1_collate_results.R
```

Then figures can be produced the R Markdown file in simulation1_plotresults.Rmd

Likewise, for simulation 2

```
for i in {1..60}; do
	Rscript simulation3.R $i
done 

Rscript simulation3_collate_results.R
```

