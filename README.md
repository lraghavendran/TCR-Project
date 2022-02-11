# TCR-Project
 Codes for TCR project bioinformatics analysis
 
 Edit the path in the first few lines of the code
 
 Run in this order
 1. make_network_cluster.R # Makes network analysis and cluster identification from Hamming's distances file
 2. make_cluster_wise_fasta  # Not necessary for next step, creates fasta file to input for Tepitool match predictions
 3. make_master_table.R     # Makes a consolidated table of CDR3 sequence characteristics, tcrmatch predictions, tepitool predictions

The master table output is not included in this repository (yet), since it is big. I will host it another website and give the link for visualization that is to follow next.
