FastSparCC_Server_Function <- function(otu_table_input, output_folder, fastSpar_iterations, bootstrap_num, threads = 4) {
  
  system("touch tmp_fastspar_script.sh")
  
  filename <- gsub(pattern = "\\.tsv", replacement = "", 
                   x = gsub(pattern = ".*\\/.*\\/", replacement = "", x = otu_table_input))
  
  text_file <- file("tmp_fastspar_script.sh")
  
  writeLines(con = text_file, 
             c(paste0("source condaX86"),
               paste0("conda activate"),
               paste0(""),
               paste0("mkdir -p ", output_folder, "/output_files/"),
               paste0(""),
               paste0("fastspar --iterations ", fastSpar_iterations,
                      " --threads ", threads,
                      " --otu_table ", otu_table_input,
                      " --correlation ", output_folder, "/output_files/fastspar_cor_", filename, ".tsv",
                      " --covariance ", output_folder, "/output_files/fastspar_cov_", filename, ".tsv"),
               paste0(""),
               paste0("mkdir -p ", output_folder, "/bootstrap_counts/", filename),
               paste0(""),
               paste0("fastspar_bootstrap",
                      " --number ", bootstrap_num,
                      " --threads ", threads,
                      " --otu_table ", otu_table_input,
                      " --prefix ", output_folder, "/bootstrap_counts/", filename, "/boot_", filename),
               paste0(""),
               paste0("mkdir -p ", output_folder, "/bootstrap_correlation/", filename),
               paste0(""),
               paste0("parallel fastspar --threads ", threads,
                      " --otu_table {} ",
                      "--correlation ", output_folder, "/bootstrap_correlation/", filename, "/cor_{/} ",
                      "--covariance ", output_folder, "/bootstrap_correlation/", filename, "/cov_{/} ",
                      "-i 5 ::: ", output_folder, "/bootstrap_counts/", filename, "/*"),
               paste0(""),
               paste0("fastspar_pvalues", 
                      " --permutations ", bootstrap_num,
                      " --threads ", threads,
                      " --otu_table ", otu_table_input,
                      " --correlation ", output_folder, "/output_files/fastspar_cor_", filename, ".tsv",
                      " --prefix ", output_folder, "/bootstrap_correlation/", filename, "/cor_boot_", filename, "_",
                      " --outfile ", output_folder, "/output_files/fastspar_pvalues_", filename, ".tsv")))
  
  close(text_file)
  
  system("bash tmp_fastspar_script.sh")
  
  system(paste0("rm -r ", output_folder, "/bootstrap_correlation/", filename))
  system(paste0("rm -r ", output_folder, "/bootstrap_counts/", filename))

}
