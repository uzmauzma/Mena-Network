# Script to prepare the file for Cytoscape, displaying node attributes.
# This script classifies taxa into "peripheral_node," "connector_node," 
# "Module_hubs," and "Network_hubs" based on the file generated from the
# "centrality" and "fast greedy modularity" files obtained from MENA. 
# This classification is essential for generating the Node attribute file required for 
# visualization and clustering of taxa in Cytoscape.
# Author: Uzma
# Date: [24th April 2023 | Version: 1.0]

Node_attribute<- read.table("data/23/Mena_ready_file_week2323 0.430 centrality.txt", header = TRUE, sep = "\t")
module<- read.table("data/23/Mena_ready_file_week2323 0.430 fast_greedy modularity.txt", header = TRUE, sep = "\t",skip = 2)

# Extract the last three columns from module file
columns <- module[, c("No..module", "Zi", "Pi")]
# columns <- module[, c(module$No..module,module$Zi, module$Pi)]

# Now Combine the three columns with centrality file
Node_attribute <- cbind(Node_attribute, columns)

Z<-Node_attribute$Zi
P<-Node_attribute$Pi
OTUs<-c(length(Z))

# data<-data[complete.cases(data$Zi<= 2.5) & (data$Pi<= 0.62),]
for (l in 1:length(Z)){
     z<-Node_attribute$Zi[l]
     p<-Node_attribute$Pi[l]
     if (z <= 2.5 && p <= 0.62) {
       OTUs[l]<-"peripheral_node"
       }
     else if (z <= 2.5 && p > 0.62){
       OTUs[l]<-"connectors_node"
       }
     else if(z > 2.5 && p <= 0.62){
       OTUs[l]<-"Modules_hubs"
       }
     else  if(z > 2.5 && p > 0.62){
       OTUs[l]<-"Network_hubs"
       }
    
}

# Give name to the generated column of node 
Node_class<-list(OTUs=OTUs)
names(Node_class)<-"Nodes_classification"
Node_attribute <- cbind(Node_attribute, Node_class)

write.csv(Node_attribute,paste("Node_attribute23",".csv",sep=""),row.names=FALSE)




