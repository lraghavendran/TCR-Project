### This script requires the output from make_network_cluster.R script. 
### For existing analysis, the required file (sample_clus_membership.txt) from the make_network_cluster.R script is stored in data folder
### This script will create a master table compiling CDR3 sample sequence characters like, clone fraction, clone count,
### ic50, percentile rank, matching epitopes, source organism and the antigen.

library('dplyr')
library('tidytext')

### Define the input and output folders
data_folder <- "/Users/Raghav/Pine Biotech/Epitope Project Elia/epitope_diversity/IEDB/Github/TCR-Project/data/"
output_folder <- "/Users/Raghav/Pine Biotech/Epitope Project Elia/epitope_diversity/IEDB/Github/TCR-Project/output/"

### Load the Hamming's distance matrix for all the samples
input_file_name = paste0(data_folder,"ALL_Alignment_distance_matrix.csv")
ALL_distance <- read.csv(input_file_name, header = TRUE, row.names=1)

### Load the cluster number for each sample from the data folder
input_file_name = paste0(data_folder,"sample_clus_membership.txt")
ceb_df <- read.table(input_file_name, sep='\t', header=TRUE)

### Load the clone fraction and clone count
input_file_name = paste0(data_folder,"All_AAsequences_clone_fraction.txt")
ALL_clone_fraction <- read.csv(input_file_name, header = TRUE, sep="\t")

### Load the TCR match results
input_file_name = paste0(data_folder,"tcrmatch_results/tcrmatch_result_trim_cluster_1-22.tsv")
tcrmatch_ceb <- read_tsv(input_file_name)  

### Load the Tepitool results
### compile cluster-wise results from tepitool to a dataframe
i=1
tepi_tool_ceb = c()
for(i in 1:max(ceb$membership)){
  tepi_result_temp = paste0(data_folder,"tepitool_results/tepitool_complete_results_cluster_",i,".csv")
  tepi_seq_ref_temp = paste0(data_folder,"tepitool_results/tepitool_seq_ref_cluster_",i,".txt")
  tepi_result <- read_csv(tepi_result_temp)
  tepi_seq_ref <- read_tsv(tepi_seq_ref_temp)
  tepi_result_join <- left_join(tepi_result, tepi_seq_ref, by = c("seq_num" = "id")) %>% dplyr::select(Seq, Sequence, everything()) %>% rename(Seq = "Sample")
  tepi_result_join$Cluster = i
  tepi_result_join <- dplyr::select(tepi_result_join, Sample, Sequence, Cluster, everything())
  if(i==1){
    tepi_tool_ceb <- tepi_result_join
  }else {
    tepi_tool_ceb = rbind(tepi_tool_ceb, tepi_result_join)
  }
}

### merge clone fraction and clone count info with tepi_tool_ceb data frame
tepi_tool_ceb.clone_frac <- left_join(tepi_tool_ceb, ALL_clone_fraction, by = c("Sample" = "SampleID")) %>% dplyr::select(-seq_num)

#### House keeping for TCR match results

### TCRmatch results (tcrmatch_ceb) contain only the part of the CDR3 sequences that was used for prediction
### To include sampleids in the tcrmatch_ceb so that the dataframes can be merged,
### we have to create a unique sequence list from tcrmatch_ceb
### and match all sampleIDs from ALL_clone_fraction file to the input sequences in TCRMatch results
unq_tcrmatch_ceb_inp_seq <- unique(tcrmatch_ceb$input_sequence)

tcrmatch_sample_id = c()
for(i in 1:length(unq_tcrmatch_ceb_inp_seq)){
  seq = unq_tcrmatch_ceb_inp_seq[i]
  tempstr = ALL_clone_fraction[grepl(seq,ALL_clone_fraction$aaSeqCDR3),]$SampleID %>% unique()
  if(length(tempstr)==0){tempstr = "NA"}
  if(length(tempstr)>1){
    for(j in 1:length(tempstr)){
      sample_id = tempstr[j]
      tempstr1 = cbind(seq, sample_id)
      tcrmatch_sample_id = rbind(tcrmatch_sample_id, tempstr1)
    }
  }else {
    tempstr = cbind(seq,tempstr)
    tcrmatch_sample_id = rbind(tcrmatch_sample_id, tempstr)
  }
}

### Merge tcrmatch with sample ids to All clone fraction with CDR3 sequences
tcrmatch_sample_id_temp <- ALL_clone_fraction %>% dplyr::select(SampleID, aaSeqCDR3) %>% right_join(as_tibble(tcrmatch_sample_id), by = c("SampleID" = "tempstr")) %>% rename(seq = "input_sequence")

### Merge the above dataframe with TCRmatch results
tcrmatch_ceb_sampleid <- right_join(as_tibble(tcrmatch_sample_id_temp), tcrmatch_ceb, by = c("input_sequence" = "input_sequence"))

### Merge TCR match results with Tepi tool results to create a master table with all the predictions combined
tcrmatch_tepi_tool_ceb_df_clone_frac_sampleid <- full_join(tcrmatch_ceb_sampleid, tepi_tool_ceb.clone_frac, by = c("SampleID" = "Sample", "aaSeqCDR3" = "aaSeqCDR3"))

### For unknown reasons, the cluster information, clone counts and fraction information is not populated as whole
### this repopulates these information.
tcrmatch_tepi_tool_ceb_df_clone_frac_sampleid <- tcrmatch_tepi_tool_ceb_df_clone_frac_sampleid %>% 
                                                    dplyr::select(-c("Cluster", "cloneCount", "cloneFraction")) %>% 
                                                      left_join(ceb_df, by = c("SampleID" = "name"))  %>% 
                                                        rename(membership = "Cluster") %>% 
                                                          left_join(ALL_clone_fraction, by = c("SampleID" = "SampleID", "aaSeqCDR3" = "aaSeqCDR3"))


### Must be executed only once
### Change NA values in the source_organism to text format " NA " so it is easy to export
tcrmatch_tepi_tool_ceb_df_clone_frac_sampleid[is.na(tcrmatch_tepi_tool_ceb_df_clone_frac_sampleid$source_organism), ]$source_organism <- " NA "

### Change no entries observed in source_organism column to " NA " values
tcrmatch_tepi_tool_ceb_df_clone_frac_sampleid[(tcrmatch_tepi_tool_ceb_df_clone_frac_sampleid$source_organism == ","),]$source_organism <- " NA "

### Add Group column to identify control and treatment samples
tcrmatch_tepi_tool_ceb_df_clone_frac_sampleid$Group = ifelse(grepl("M", tcrmatch_tepi_tool_ceb_df_clone_frac_sampleid$SampleID), 'Treatment', 'Control')

### write out the master output table
out_file_name = paste0(data_folder,"tcrmatch_tepitool_master_table.txt")
write.table(tcrmatch_tepi_tool_ceb_df_clone_frac_sampleid, file = out_file_name, sep = "\t ",quote = FALSE,  row.names = FALSE, col.names = TRUE)



