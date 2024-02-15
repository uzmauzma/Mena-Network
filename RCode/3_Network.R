
# Input file: This script utilizes the upper triangular matrix of the Pearson correlation file
# generated from the MENA pipeline (created by "2_Generate_corr.R").
# Output file: The output file consists of two columns representing the OTUs, and the third column 
# representing their correlation values,
# which will be used for network generation.
#Author: Uzma  
#Email:uzma.k.khan@glasgow.ac.uk ; uzma.1415@gmail.com
#Date: 15 Feb 2024 | Version:1.0

data<- read.table("Upper_triangular_mat23.txt", header = TRUE, sep = "\t")

data<-reshape2::melt(data,varnames = FALSE)

# The following lines of code remove the rows with null values and display only the significant correlations. 
# This is necessary because the MENA pipeline generates correlations for all OTUs, but we are interested in
# OTUs that exhibit significant correlations as mentioned in the similarity matrix.

filtered_data<-data[complete.cases(data$value),]

Network_data<-filtered_data[!filtered_data$value %in% c("0","1"),]


colnames(Network_data)[1]<-"Variable1"
colnames(Network_data)[2]<-"Variable2"

# co<-unique(c(as.character(Network_data$Variable1),as.character(Network_data$Variable2)))

#  # The file below represents the network that we will import into Cytoscape.
write.csv(Network_data,paste("Network23",".csv",sep=""),row.names = FALSE)

# Manually replace the dots with hyphens in the "network.csv" file and the "node_attributes.csv" file.
