#### This script runs FastSpar on the DOM dataset ####
#
# Author: Felix Milke
# Date: 19.11.2026

library(tidyverse)

source("https://raw.githubusercontent.com/dermilke/ExCom/master/R/Datalist_Wrangling_Functions.R")
source("https://raw.githubusercontent.com/dermilke/ExCom/master/R/Phylogenetic_Functions.R")
source("https://raw.githubusercontent.com/dermilke/ExCom/master/R/Stats_Diversity.R")

source("R/FastSparCC_Server_Function_Studio.R")

bootstrap_num = 500
fastSpar_iterations = 25 
threads = 10

FastSparCC_Server_Function(file_input = "output/FastSpar/FastSpar_Dataprep_DOM.tsv", 
                           file_output = "output/FastSpar", 
                           fastSpar_iterations = fastSpar_iterations, 
                           bootstrap_num = bootstrap_num, 
                           threads = threads)
  
unlink("output/FastSpar/bootstrap_counts", recursive = T)
unlink("output/FastSpar/bootstrap_correlation", recursive = T)
