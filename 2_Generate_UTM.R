# Script to Generate Upper Triangular Matrix from "Pearson Correlation" File Generated with MENA
# The MENA software generates a similarity matrix using the Pearson Correlation Coefficient, 
# which quantifies the linear relationship between pairs of variables.
# In order to determine significant relationships within the network, 
# Random Matrix Theory was employed for threshold detection.
# Author:uzma
# Date: [24 April 2023 | version:1.0 ]

# data<-read.csv("../dates/Mena_ready_file_t_deapth_relative_per Pearson Correlation.csv",header=F)
data <- read.table("BoT/MenareadyBOT Pearson Correlation.txt", header = FALSE,sep = "\t",fill = TRUE)

# # Read the "MV Estimated.txt" file to retrieve the Taxa names since they are not present in the Pearson correlation matrix.
taxa_name<- read.table("BoT/MenareadyBOT MV Estimated.txt", header = FALSE, sep = "\t")

N_rows <- nrow(data)
N_cols <- ncol(data)

data<-as.matrix(data)

mat = matrix("NaN",N_rows, N_cols)
# Number of rows and cols of correlation matrix
NR<-nrow(data)
j<-0

for (i in 1:NR){
     j<-j+1
    if (j>NR){
       break
    }
    # else if(j==2){
    #        j<-j+1
    #        mat[i,j:NR]<-data[i,1:(NR-abs(1-j))]
    #      }
    #  else {
          
           mat[i,j:NR]<-data[i,1:(NR-abs(1-j))]
          # }
    }

# print(mat)

library(data.table)



# Assign taxa names from "MV Estimated.txt" file generated by MENA

Taxon<-taxa_name$V1

colnames(mat)<-Taxon
col_names=Taxon

mat <- cbind(col_names, mat)



# write.csv(mat,paste("Upper_triangular_mat",".csv",sep=""))
write.table(mat, "Upper_triangular_matBoT.txt", sep = "\t", row.names = FALSE, quote = FALSE)
