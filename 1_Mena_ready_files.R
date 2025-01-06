# Script for MENA-Ready File
# This script prepares the data file for use with the MENA software.
# MENA software Reference: http://129.15.40.240/mena/
# Author: Uzma; PDRA[Water and Envirnomental Group: Glasgow] 
# Date: [24 April 2023 | version:1.0 ]


rm(list=ls())

library(phyloseq)
library(ggplot2)
library(viridis)
library(microbiome)
library(RColorBrewer)
library(cowplot)


#PARAMETERS ###########################
library(biomformat);b_<-read_biom("Data/feature_w_tax.biom");physeq<-merge_phyloseq(otu_table(as(biom_data(b_),"matrix"),taxa_are_rows=TRUE),tax_table(as(observation_metadata(b_),"matrix")))
meta_table<-read.csv("Data/meta_table.csv",header=T,row.names=1)
which_level<-"Family" #Phylum Class Order Family Genus Otus
height_image_heatmap=16
legend_text_size=6
legend_title_size=8
what_detection="absolute" #absolute relative
taxa_pres_sample=0.70   # 0.75 for K1 & K2, 0.62 for K3
text_size=12
#/PARAMETERS ###########################

abund_table<-otu_table(physeq)
abund_table<-t(abund_table)
#Uncomment if you'd like to get rid of samples below a certain library size
# abund_table<-abund_table[rowSums(abund_table)>=5000,]


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


#/DATA IMPORT############################################################


#PARAMETERS CHANGE THE GROUPING COLUMN AS YOU DESIRE############################

#########COLLATE FRACTIONS TOGETHER#############
# 29 for K1, & K3, whilst 14 is for K2

meta_table<-meta_table[(meta_table$section) %in% c("TOP"),]
#meta_table<-meta_table[(meta_table$section) %in% c("BOT"),]


# day="29"  # 14  29
# activity="C13"
# type="DNA"
# site="K1" # K1 ,K2,K3
# meta_table=meta_table[meta_table$Fraction != "" & meta_table$Activity == activity & meta_table$Day == day & meta_table$Type == type & meta_table$Site == site,]
#select only the natural gradient
#Gradient="natural_gradient"
# meta_table=meta_table[c(1 ,3 ,5 ,7 ,10 ,13 ,17 ,19 ,21),]
# abund_table<-abund_table[rownames(meta_table),]



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

# calculate the relative abundance
# new_abund_table<-as.data.frame(as(new_abund_table,"matrix"))
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
write.table(abund_table,paste("Mena_ready_file_TOP",".txt"), sep = "\t", row.names = FALSE, quote = FALSE)
