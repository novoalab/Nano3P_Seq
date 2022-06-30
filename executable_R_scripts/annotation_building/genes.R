# Build a BED File for small RNA annotation ends

# Reading arguments from command line
args = commandArgs(trailingOnly=TRUE)

#Arguments
data.input <- read.delim(args[1],header=FALSE)  #1st variable
label <- as.character(args[2])  #1st label


#Libraries 
library(stringr)


#For Mouse Reference
data <- read.csv(args[1], header=FALSE, sep="\t")
data2<- subset(data, V3=="gene")
column<- str_split_fixed(data2$V9, ";",10)
gene_name <- column[,3]
gene_name2 <- str_split_fixed(gene_name, " ",3)
data2$gene_name<- gene_name2[,3]
gene_type <- column[,5]
gene_type2 <- str_split_fixed(gene_type, " ",3)
data2$gene_type<- gene_type2[,3]
data3 <- data2[,c("V1", "V4", "V5", "gene_name", "gene_type", "V7")]


#Export the final table
write.table(data3, file=paste(label, "genes.bed", sep="_"), quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)