### This script requires the output from make_network_cluster.R script. 
### For existing analysis, the required file (sample_clus_membership.txt) from the make_network_cluster.R script is stored in data folder
### This script will create the fasta files of CDR3 sequences for samples in each cluster

### Install Libraries
install.packages("dplyr", dependencies = TRUE)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Biostrings")

### Load libraries
library('dplyr')
library("Biostrings")

### Define the input and output folders
data_folder <- "/Users/Raghav/Pine Biotech/Epitope Project Elia/epitope_diversity/IEDB/Github/TCR-Project/data/"
fasta_output_folder <- "/Users/Raghav/Pine Biotech/Epitope Project Elia/epitope_diversity/IEDB/Github/TCR-Project/fasta_output/"

### Load the cluster number for each sample from the data folder
input_file_name = paste0(data_folder,"sample_clus_membership.txt")
ceb_df <- read.table(input_file_name, sep='\t', header=TRUE)

### Load the CDR3 Amino acid sequences for the samples
input_file_name = paste0(data_folder,"All_AAsequences3.fasta")
ALL_aaseq <- readAAStringSet(input_file_name)

### output fasta files for each cluster
for(i in 1:max(ceb$membership)){
  out_file_name = paste0(fasta_output_folder,"All_AAsequences_ceb_clus_",i,".fasta")
  writeXStringSet(ALL_aaseq[ceb$names[ceb$membership==i]], out_file_name )
}

### Use these fasta files as input for Tepi Tool prediction of ic50(binding affinity) and percentile rank (specificity)
### Also use the All_AAsequences3.fasta file as input for TCRMatch prediction of matching epitopes from source organism and its antigen

