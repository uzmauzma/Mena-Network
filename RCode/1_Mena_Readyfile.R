#Script: 1_Mena_Readyfile.R
#Description: This script takes a combined dataset as input and prepares it for analysis with Molecular Ecological Network Analysis (MENA) by generating the Meta-ready file.
#Author: Uzma  
#Email:uzma.k.khan@glasgow.ac.uk ; uzma.1415@gmail.com
#Date: 15 Feb 2024 | Version:1.0
# Usage:
# Ensure your combined dataset is in the correct format.
# Execute this script and provide the path to your dataset as input.
# The script will preprocess the data and generate the Meta-ready file required for MENA analysis. 
# Output:
# Meta-ready file for use with MENA.
#Note:
# If the generated dataset containing the top 75% most abundant taxa has fewer than 50 taxa, MENA will not be applied, and the threshold will be adjusted accordingly
# Adjustments to the minimum prevalence are based on the dataset size to ensure effective analysis with MENA
#For more information about MENA, visit: http://129.15.40.240/mena/
  
rm(list=ls())

library(phyloseq)
library(ggplot2)
library(viridis)
library(microbiome)
library(RColorBrewer)
library(cowplot)


#PARAMETERS ###########################
library(biomformat);b_<-read_biom("../../../Analysis/create_combined_biom_file_16S_18S/collated_feature_w_tax.biom");physeq<-merge_phyloseq(otu_table(as(biom_data(b_),"matrix"),taxa_are_rows=TRUE),tax_table(as(observation_metadata(b_),"matrix")))
meta_table<-read.csv("../../../Analysis/create_combined_biom_file_16S_18S/collated_meta_table.csv",header=T,row.names=1)
which_level<-"Genus" #Phylum Class Order Family Genus Otus
height_image_heatmap=16
legend_text_size=6
legend_title_size=8
what_detection="absolute" #absolute relative
minimum_prevalence=0.75    # 0.85
text_size=12
#/PARAMETERS ###########################

abund_table<-otu_table(physeq)
abund_table<-t(abund_table)
#Uncomment if you'd like to get rid of samples below a certain library size
abund_table<-abund_table[rowSums(abund_table)>=5000,]


OTU_taxonomy<-as.data.frame(tax_table(physeq))
colnames(OTU_taxonomy)<-c("Kingdom","Phylum","Class","Order","Family","Genus","Otus")

#Ensure that all columns of OTU_taxonomy are character and not factors
OTU_taxonomy[] <- lapply(OTU_taxonomy, function(x) as.character(x))
OTU_taxonomy[is.na(OTU_taxonomy)]<-""
OTU_taxonomy$Otus<-gsub("D_6__|s__","",OTU_taxonomy$Otus)
OTU_taxonomy$Genus<-gsub("D_5__|g__","",OTU_taxonomy$Genus)
OTU_taxonomy$Family<-gsub("D_4__|f__","",OTU_taxonomy$Family)
OTU_taxonomy$Order<-gsub("D_3__|o__","",OTU_taxonomy$Order)
OTU_taxonomy$Class<-gsub("D_2__|c__","",OTU_taxonomy$Class)
OTU_taxonomy$Phylum<-gsub("D_1__|p__","",OTU_taxonomy$Phylum)
OTU_taxonomy$Kingdom<-gsub("D_0__|d__","",OTU_taxonomy$Kingdom)

#Remove singletons and adjust OTU_taxonomy
abund_table<-abund_table[,colSums(abund_table)>1]
OTU_taxonomy<-OTU_taxonomy[colnames(abund_table),]

#get rid of contaminants with "Unassigned", "Chloroplast" and "Mitochondria" assignment", and "non classified" at Phylum level
abund_table<-abund_table[,!(OTU_taxonomy$Kingdom %in% c("Unassigned") | OTU_taxonomy$Phylum=="" | OTU_taxonomy$Order %in% c("Chloroplast") | OTU_taxonomy$Family %in% c("Mitochondria"))]
#extract subset of abund_table for which samples also exists in meta_table
abund_table<-abund_table[rownames(abund_table) %in% rownames(meta_table),]
#when reducing the abund_table, there is a high likelihood that an OTU was only present in a sample that is removed, so we shrink
#the abund_table to get rid of empty columns
abund_table<-abund_table[,colSums(abund_table)>0]
#make your meta_table smaller by only considering samples that appear in abund_table
meta_table<-meta_table[rownames(abund_table),]
#make OTU_taxonomy smaller by only considering OTUs that appear in abund_table
OTU_taxonomy<-OTU_taxonomy[colnames(abund_table),]
#At this point we have abund_table, meta_table, and OTU_taxonomy are ready and their dimensions should match
#/DATA IMPORT############################################################

# Your Hypoythesis 
#PARAMETERS CHANGE THE GROUPING COLUMN AS YOU DESIRE############################

meta_table<-meta_table[!(meta_table$Column_size) %in% c("30", "60"),]
#meta_table<-meta_table[(meta_table$Week) %in% c("5"),]
meta_table<-meta_table[(meta_table$Week) %in% c("5","9","12","23"),]
meta_table<-meta_table[(meta_table$Sample_type) %in% c("Influent", "GAC", "Effluent"),]
meta_table<-meta_table[!meta_table$Replicate %in% c("T0"),]



#PARAMETERS CHANGE THE GROUPING COLUMN AS YOU DESIRE############################

#Adjust abund_table to contain only those rows that got selected in the Hypothesis space
abund_table<-abund_table[rownames(meta_table),]
#After adjustment, get rid of OTUs that are all empty
abund_table<-abund_table[,colSums(abund_table)>0]
#Adjust OTU taxonomy
OTU_taxonomy<-OTU_taxonomy[colnames(abund_table),]

#COLLATE OTUS AT A PARTICULAR LEVEL#######################################
new_abund_table<-NULL
if(which_level=="Otus"){
  new_abund_table<-abund_table
} else {
  list<-unique(OTU_taxonomy[,which_level])
  new_abund_table<-NULL
  for(i in list){
    tmp<-data.frame(rowSums(abund_table[,rownames(OTU_taxonomy)[OTU_taxonomy[,which_level]==i],drop=FALSE]))
    if(i==""){colnames(tmp)<-c("__Unknowns__")} else {
      #colnames(tmp)<-paste("",i,sep="")
      colnames(tmp)<-gsub(";+$","",paste(sapply(OTU_taxonomy[OTU_taxonomy[,which_level]==i,][1,1:which(colnames(OTU_taxonomy)==which_level)],as.character),collapse=";"))
    }
    if(is.null(new_abund_table)){new_abund_table<-tmp} else {new_abund_table<-cbind(tmp,new_abund_table)}
  }
}



#######################  Top 75% most abundant taxa ############################
abund_table<-new_abund_table
#  Filter out thoes OTUS which are present in 75% of samples
New_abund_table<-abund_table
New_abund_table[New_abund_table>0]=1
sample_avg<-colSums(New_abund_table)/nrow(New_abund_table)
abund_table<-abund_table[,colnames(abund_table) %in% names(sample_avg[sample_avg>=taxa_pres_sample])] 
#Adjust OTU taxonomy
OTU_taxonomy<-OTU_taxonomy[colnames(abund_table),]


# calculate the relative abundance
sample_totals <- rowSums(abund_table)
abund_table<-(abund_table / sample_totals)*100

abund_table<-t(abund_table)

# Assign the "Samp_name" as the row name of the abund table .

Sample_names=rownames(abund_table)
abund_table <- cbind(Sample_names, abund_table)



# Now file ready to import into MENA tool  (OTUS * samples)
#write.csv(OTU_taxonomy,paste("OTU_taxonomy_AOA_AOBk123_",".csv",sep=""))
write.table(abund_table,paste("Mena_ready_file",".txt"), sep = "\t", row.names = FALSE, quote = FALSE)




