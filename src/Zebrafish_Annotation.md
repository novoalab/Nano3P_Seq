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
gene_name <- column[,1]
gene_name2 <- str_split_fixed(gene_name, " ",2)
data2$gene_name<- gene_name2[,2]
gene_type <- column[,5]
gene_type2 <- str_split_fixed(gene_type, " ",3)
data2$gene_type<- gene_type2[,3]
data3 <- data2[,c("V1", "V4", "V5", "gene_name", "gene_type", "V7")]

#Ribosomal RNA
sequin_bed <- read.delim("Sequin_4_mix_Annotation.bed",sep="\t", header=FALSE)
colnames(sequin_bed) <-c("V1", "V4", "V5", "gene_name", "gene_type", "V7")


both <- rbind(data3,sequin_bed)

write.table(both, file="Danio_rerio_11_Sequin4_GENE.bed", sep="\t", quote=FALSE,col.names=FALSE, row.names=FALSE)
```

#Manipulate GTF Files to export Gene Only BED files
```R
#Libraries 
library(stringr)
#For Mouse Reference
data <- read.csv("/users/enovoa/boguzhan/references/danio_rerio/Danio_rerio.GRCz11.103.2.gtf", header=FALSE, sep="\t")
data2<- subset(data, V3=="exon")
column<- str_split_fixed(data2$V9, ";",10)
gene_name <- column[,1]
gene_name2 <- str_split_fixed(gene_name, " ",2)
data2$gene_name<- gene_name2[,2]
gene_type <- column[,8]
gene_type2 <- str_split_fixed(gene_type, " ",3)
data2$gene_type<- gene_type2[,3]
data3 <- data2[,c("V1", "V4", "V5", "gene_name", "gene_type", "V7")]

#Ribosomal RNA
sequin_bed <- read.delim("Sequin_4_mix_Annotation.bed",sep="\t", header=FALSE)
colnames(sequin_bed) <-c("V1", "V4", "V5", "gene_name", "gene_type", "V7")


both <- rbind(data3,sequin_bed)

write.table(both, file="Danio_rerio_11_Sequin4_EXON.bed", sep="\t", quote=FALSE,col.names=FALSE, row.names=FALSE)
```


## Build a BED File for small RNA annotation ends GENE
```R
#Libraries 
library(stringr)
#For Mouse Reference
data <- read.csv("/users/enovoa/boguzhan/references/danio_rerio/Danio_rerio.GRCz11.103.2.gtf", header=FALSE, sep="\t")
data2<- subset(data, V3=="gene")
column<- str_split_fixed(data2$V9, ";",10)
gene_name <- column[,1]
gene_name2 <- str_split_fixed(gene_name, " ",2)
data2$gene_name<- gene_name2[,2]
gene_type <- column[,5]
gene_type2 <- str_split_fixed(gene_type, " ",3)
data2$gene_type<- gene_type2[,3]
data3 <- data2[,c("V1", "V4", "V5", "gene_name", "gene_type", "V7")]
data4 <- subset(data3, gene_type=="snRNA" | gene_type=="snoRNA" | gene_type=="scaRNA" )
write.table(data4, file="Danio_rerio_11_SmallRNA_Gene.bed", sep="\t", quote=FALSE,col.names=FALSE, row.names=FALSE)

```

## Build a BED File for small RNA annotation ends EXONS
```R
#Libraries 
library(stringr)
#For Mouse Reference
data <- read.csv("/users/enovoa/boguzhan/references/danio_rerio/Danio_rerio.GRCz11.103.2.gtf", header=FALSE, sep="\t")
data2<- subset(data, V3=="exon")
column<- str_split_fixed(data2$V9, ";",10)
gene_name <- column[,1]
gene_name2 <- str_split_fixed(gene_name, " ",2)
data2$gene_name<- gene_name2[,2]
gene_type <- column[,8]
gene_type2 <- str_split_fixed(gene_type, " ",3)
data2$gene_type<- gene_type2[,3]
data3 <- data2[,c("V1", "V4", "V5", "gene_name", "gene_type", "V7")]
data4 <- subset(data3, gene_type=="snRNA" | gene_type=="snoRNA" | gene_type=="scaRNA" )
write.table(data4, file="Danio_rerio_11_SmallRNA_Exon.bed", sep="\t", quote=FALSE,col.names=FALSE, row.names=FALSE)

```


```R
#Manipulate the GTF File
data <- read.csv("/users/enovoa/boguzhan/references/danio_rerio/Danio_rerio.GRCz11.103.2.gtf", header=FALSE, sep="\t")
gtf_transcript <- subset(data , V3=="transcript")
column<- str_split_fixed(gtf_transcript$V9, ";",10)
gene_type <- column[,7]
gene_type2 <- str_split_fixed(gene_type, " ",3)
gtf_transcript$gene_type<- gene_type2[,3]
gtf_transcript <- subset(gtf_transcript, gene_type=="snRNA" | gene_type=="snoRNA" | gene_type=="scaRNA" )
gtf_trans_neg <- subset(gtf_transcript, V7=="-")
gtf_trans_neg2 <- gtf_trans_neg[,c("V1", "V4", "V7")]
colnames(gtf_trans_neg2) <- c("Chr", "Position1", "Strand")
gtf_trans_pos <- subset(gtf_transcript, V7=="+")
gtf_trans_pos2 <- gtf_trans_pos[,c("V1", "V5", "V7")]
colnames(gtf_trans_pos2) <- c("Chr", "Position1", "Strand")
merged <- rbind(gtf_trans_neg2,gtf_trans_pos2 )
merged$Position2 <- merged$Position1+5
merged$Position1 <- merged$Position1-5
final <- merged[,c("Chr","Position1", "Position2", "Strand")]
write.table(final, file="Danio_rerio_11_SmallRNA_TranscriptEnds.bed", sep="\t", quote=FALSE,col.names=FALSE, row.names=FALSE)
```




```R
#Libraries 
library(stringr)
#For Mouse Reference
data <- read.csv("/users/enovoa/boguzhan/references/danio_rerio/Danio_rerio.GRCz11.103.2.gtf", header=FALSE, sep="\t")
data2<- subset(data, V3=="gene")
column<- str_split_fixed(data2$V9, ";",10)
gene_name <- column[,1]
gene_name2 <- str_split_fixed(gene_name, " ",2)
data2$gene_name<- gene_name2[,2]
gene_type <- column[,5]
gene_type2 <- str_split_fixed(gene_type, " ",3)
data2$gene_type<- gene_type2[,3]
data3 <- data2[,c("V1", "V4", "V5", "gene_name", "gene_type", "V7")]
data4 <- subset(data3, gene_type=="miRNA")
write.table(data4, file="Danio_rerio_11_miRNA_Gene.bed", sep="\t", quote=FALSE,col.names=FALSE, row.names=FALSE)
```

#Manipulate GTF Files to export Gene Only BED files
```R
#Libraries 
library(stringr)
#For Mouse Reference
data <- read.csv("/users/enovoa/boguzhan/references/danio_rerio/Danio_rerio.GRCz11.103.2.gtf", header=FALSE, sep="\t")
data2<- subset(data, V3=="exon")
column<- str_split_fixed(data2$V9, ";",10)
gene_name <- column[,1]
gene_name2 <- str_split_fixed(gene_name, " ",2)
data2$gene_name<- gene_name2[,2]
gene_type <- column[,8]
gene_type2 <- str_split_fixed(gene_type, " ",3)
data2$gene_type<- gene_type2[,3]
data3 <- data2[,c("V1", "V4", "V5", "gene_name", "gene_type", "V7")]
data4 <- subset(data3,  gene_type!="miRNA" & gene_type!="snRNA" & gene_type!="snoRNA" & gene_type!="scaRNA"   )


#Ribosomal RNA
sequin_bed <- read.delim("Sequin_4_mix_Annotation.bed",sep="\t", header=FALSE)
colnames(sequin_bed) <-c("V1", "V4", "V5", "gene_name", "gene_type", "V7")


both <- rbind(data4,sequin_bed)

write.table(both, file="Danio_rerio_11_Sequin4_Rest_EXON.bed", sep="\t", quote=FALSE,col.names=FALSE, row.names=FALSE)
```



#Manipulate GTF Files to export Gene Only BED files
```R
#Libraries 
library(stringr)
#For Mouse Reference
data <- read.csv("/users/enovoa/boguzhan/references/danio_rerio/Danio_rerio.GRCz11.103.2.gtf", header=FALSE, sep="\t")
data2<- subset(data, V3=="gene")
column<- str_split_fixed(data2$V9, ";",10)
gene_name <- column[,1]
gene_name2 <- str_split_fixed(gene_name, " ",2)
data2$gene_name<- gene_name2[,2]
gene_type <- column[,5]
gene_type2 <- str_split_fixed(gene_type, " ",3)
data2$gene_type<- gene_type2[,3]
data3 <- data2[,c("V1", "V4", "V5", "gene_name", "gene_type", "V7")]
data4 <- subset(data3,  gene_type!="miRNA" & gene_type!="snRNA" & gene_type!="snoRNA" & gene_type!="scaRNA"   )

#Ribosomal RNA
sequin_bed <- read.delim("Sequin_4_mix_Annotation.bed",sep="\t", header=FALSE)
colnames(sequin_bed) <-c("V1", "V4", "V5", "gene_name", "gene_type", "V7")


both <- rbind(data4,sequin_bed)

write.table(both, file="Danio_rerio_11_Sequin4_Rest_GENE.bed", sep="\t", quote=FALSE,col.names=FALSE, row.names=FALSE)
```



# Manipulate in order to extract transcript start positions
```R
#Manipulate the GTF File
gtf <- read.csv("/users/enovoa/boguzhan/references/danio_rerio/Danio_rerio.GRCz11.103.2.gtf", header=FALSE, sep="\t")

gtf_transcript <- subset(gtf , V3=="transcript")
gtf_trans_neg <- subset(gtf_transcript, V7=="-")
gtf_trans_neg2 <- gtf_trans_neg[,c("V1", "V4", "V7")]
colnames(gtf_trans_neg2) <- c("Chr", "Position1", "Strand")
gtf_trans_pos <- subset(gtf_transcript, V7=="+")
gtf_trans_pos2 <- gtf_trans_pos[,c("V1", "V5", "V7")]
colnames(gtf_trans_pos2) <- c("Chr", "Position1", "Strand")
merged <- rbind(gtf_trans_neg2,gtf_trans_pos2 )
merged$Position2 <- merged$Position1+3
merged$Position1 <- merged$Position1-3
final <- merged[,c("Chr","Position1", "Position2", "Strand")]

write.table(final, file="Danio_rerio_11_Transcript_Ends.bed", quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)


```

