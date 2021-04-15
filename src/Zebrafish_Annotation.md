###################################################
###### ANNOTATION BUILDING OF ZEBRAFISH   ###########
###################################################
########## OGUZHAN BEGIK APRIL 2020 ###############


#Manipulate GTF Files to export Gene Only BED files
```R
#Libraries 
library(stringr)
#For Mouse Reference
data <- read.csv("/users/enovoa/boguzhan/references/danio_rerio/Danio_rerio.GRCz11.103.2.gtf", header=FALSE, sep="\t")
data2<- subset(data, V3=="gene")
column<- str_split_fixed(data2$V9, ";",10)
gene_name <- column[,3]
gene_name2 <- str_split_fixed(gene_name, " ",3)
data2$gene_name<- gene_name2[,3]
gene_type <- column[,5]
gene_type2 <- str_split_fixed(gene_type, " ",3)
data2$gene_type<- gene_type2[,3]
data3 <- data2[,c("V1", "V4", "V5", "gene_name", "gene_type", "V7")]

#Ribosomal RNA
rRNA_bed <- read.delim("Zebrafish_rRNA_Annotation.bed",sep="\t", header=FALSE)
colnames(rRNA_bed) <-c("V1", "V4", "V5", "gene_name", "gene_type", "V7")







write.table(both, file="Mus_musculus.GRCm38.102_sequinv2.2_rRNA_GENE.bed", sep="\t", quote=FALSE,col.names=FALSE, row.names=FALSE)
```