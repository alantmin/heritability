# Introduction

These are files used in the publication of Comparing Heritability Estimators under Alternative Structures of Linkage Disequilibrium by Min, Basu, and Thompson (2022).

The code assumes that `GCTA` (https://yanglab.westlake.edu.cn/software/gcta/), `LDAK` (https://dougspeed.com/ldak/), and `plink` are installed (https://zzz.bwh.harvard.edu/plink/).

The files are described below

1. `estimators.R`: The estimators are implemented in this file 
2. `gcta.R`: This file includes helper functions to run GCTA
3. `genotypes_cor.R`: This file provides functions to simulate genotypes with correlation structures. 
4. `genotypes_related.R`: This file provide functions to simulate genotypes with different levels of cousinship.
5. `make_folders.R`: This file uses the simulation functions to simulate genotypes and phenotypes using different correlation structures, then outputs a folder with results of the estimated values of heritability for various estimators.
6. `make_folders_cousins.R`: This file uses uses simulation functions to simulate genotypes and phenotypes using different levels of cousinships, then outputs a folder with results of the estimated values of heritability for various estimators. 
7. `markers1shuffle.txt`: This file provides SNPs from 1000 genomes in a random order that are filtered so that the minor allele frequency is greater than 0.05. 
8. `collate_results.R`: This file looks through all of the folders that are output by the make_folders script and puts the results into one data frame. 
8. `plotresults.Rmd`: This is a RMarkdown file that makes plots given the data frame created by `collate_results.R`

# Usage 

To run a subset of the settings, use the command 

```
Rscript make_folders.R <i>
```

or 

```
Rscript make_folders_cousins.R <i>
```

Where `i` is the index of the subset of parameters to run. By default, these commands run 50 of the selected parameters.
