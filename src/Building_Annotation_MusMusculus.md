###################################################
###### ANNOTATION BUILDING OF MOUSE.   ###########
###################################################
########## OGUZHAN BEGIK APRIL 2020 ###############

## CONVERTING GTF TO BED FILES


## Build a BED File for small RNA annotation ends
```R
#Libraries 
library(stringr)
#For Mouse Reference
data <- read.csv("/users/enovoa/boguzhan/references/mouse_with_sequins/Mus_musculus.GRCm38.102.gtf", header=FALSE, sep="\t")
data2<- subset(data, V3=="gene")
column<- str_split_fixed(data2$V9, ";",10)
gene_name <- column[,3]
gene_name2 <- str_split_fixed(gene_name, " ",3)
data2$gene_name<- gene_name2[,3]
gene_type <- column[,5]
gene_type2 <- str_split_fixed(gene_type, " ",3)
data2$gene_type<- gene_type2[,3]
data3 <- data2[,c("V1", "V4", "V5", "gene_name", "gene_type", "V7")]
data4 <- subset(data3, gene_type=="snRNA" | gene_type=="snoRNA" | gene_type=="scaRNA" )
write.table(data4, file="Mus_musculus.GRCm38.102_sequinv2.2_SmallRNA_Gene.bed", sep="\t", quote=FALSE,col.names=FALSE, row.names=FALSE)
```
```R
#Libraries 
library(stringr)
#For Mouse Reference
data <- read.csv("/users/enovoa/boguzhan/references/mouse_with_sequins/Mus_musculus.GRCm38.102.gtf", header=FALSE, sep="\t")
data2<- subset(data, V3=="gene")
column<- str_split_fixed(data2$V9, ";",10)
gene_name <- column[,3]
gene_name2 <- str_split_fixed(gene_name, " ",3)
data2$gene_name<- gene_name2[,3]
gene_type <- column[,5]
gene_type2 <- str_split_fixed(gene_type, " ",3)
data2$gene_type<- gene_type2[,3]
data3 <- data2[,c("V1", "V4", "V5", "gene_name", "gene_type", "V7")]
data4 <- subset(data3, gene_type=="miRNA")
write.table(data4, file="Mus_musculus.GRCm38.102_sequinv2.2_miRNA_Gene.bed", sep="\t", quote=FALSE,col.names=FALSE, row.names=FALSE)
```

```R
#Manipulate the GTF File
data <- read.csv("/users/enovoa/boguzhan/references/mouse_with_sequins/Mus_musculus.GRCm38.102.gtf", header=FALSE, sep="\t")
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
merged$Position2 <- merged$Position1+3
merged$Position1 <- merged$Position1-3
final <- merged[,c("Chr","Position1", "Position2", "Strand")]
write.table(final, file="Mus_musculus.GRCm38.102_sequinv2.2_SmallRNA_TranscriptEnds.bed", sep="\t", quote=FALSE,col.names=FALSE, row.names=FALSE)
```



#Manipulate GTF Files to export Gene Only BED files
```R
#Libraries 
library(stringr)
#For Mouse Reference
data <- read.csv("/users/enovoa/boguzhan/references/mouse_with_sequins/Mus_musculus.GRCm38.102.gtf", header=FALSE, sep="\t")
data2<- subset(data, V3=="gene")
column<- str_split_fixed(data2$V9, ";",10)
gene_name <- column[,3]
gene_name2 <- str_split_fixed(gene_name, " ",3)
data2$gene_name<- gene_name2[,3]
gene_type <- column[,5]
gene_type2 <- str_split_fixed(gene_type, " ",3)
data2$gene_type<- gene_type2[,3]
data3 <- data2[,c("V1", "V4", "V5", "gene_name", "gene_type", "V7")]

#For Sequin Reference
sequin <- read.csv("/users/enovoa/boguzhan/references/mouse_with_sequins/RNAsequins.v2.2.gtf", header=FALSE, sep="\t")
sequin2<- subset(sequin, V3=="gene")
column<- str_split_fixed(sequin2$V9, ";",3)
gene_type <- column[,2]
gene_type2 <- str_split_fixed(gene_type, " ",4)
sequin2$gene_type<- gene_type2[,3]
gene_name <- column[,1]
gene_name2 <- str_split_fixed(gene_name, " ",4)
sequin2$gene_name<- gene_name2[,2]
sequin3 <- sequin2[,c("V1", "V4", "V5", "gene_name", "gene_type", "V7")]


#Ribosomal RNA
rRNA_bed <- read.delim("mouse_rRNA.bed", header=FALSE)

colnames(rRNA_bed) <-c("V1", "V4", "V5", "gene_name", "gene_type", "V7")



both <- rbind(data3, sequin3, rRNA_bed)



write.table(both, file="Mus_musculus.GRCm38.102_sequinv2.2_rRNA_GENE.bed", sep="\t", quote=FALSE,col.names=FALSE, row.names=FALSE)
```


#Manipulate GTF Files to export Exon Only BED files
```R
#Libraries 
library(stringr)

#For Mouse Reference
data <- read.csv("/users/enovoa/boguzhan/references/mouse_with_sequins/Mus_musculus.GRCm38.102.gtf", header=FALSE, sep="\t")
data2<- subset(data, V3=="exon")
column<- str_split_fixed(data2$V9, ";",10)
gene_name <- column[,6]
gene_name2 <- str_split_fixed(gene_name, " ",3)
data2$gene_name<- gene_name2[,3]
gene_type <- column[,8]
gene_type2 <- str_split_fixed(gene_type, " ",3)
data2$gene_type<- gene_type2[,3]
data3 <- data2[,c("V1", "V4", "V5", "gene_name", "gene_type", "V7")]

#For Sequin Reference
sequin <- read.csv("/users/enovoa/boguzhan/references/mouse_with_sequins/RNAsequins.v2.2.gtf", header=FALSE, sep="\t")
sequin2<- subset(sequin, V3=="exon")
column<- str_split_fixed(sequin2$V9, ";",3)
sequin2$gene_type<- "synthetic"
gene_name <- column[,1]
gene_name2 <- str_split_fixed(gene_name, " ",4)
sequin2$gene_name<- gene_name2[,2]
sequin3 <- sequin2[,c("V1", "V4", "V5", "gene_name", "gene_type", "V7")]

#Ribosomal RNA
rRNA_bed <- read.delim("mouse_rRNA.bed", header=FALSE)

colnames(rRNA_bed) <-c("V1", "V4", "V5", "gene_name", "gene_type", "V7")


both <- rbind(data3, sequin3,rRNA_bed)



write.table(both, file="Mus_musculus.GRCm38.102_sequinv2.2_rRNA_EXON.bed", sep="\t", quote=FALSE,col.names=FALSE, row.names=FALSE)
```





#Manipulate GTF Files to export Exon Only BED files
```R
#Libraries 
library(stringr)

#For Mouse Reference
data <- read.csv("/users/enovoa/boguzhan/references/mouse_with_sequins/Mus_musculus.GRCm38.102.gtf", header=FALSE, sep="\t")
data2<- subset(data, V3=="exon")
column<- str_split_fixed(data2$V9, ";",10)
gene_name <- column[,6]
gene_name2 <- str_split_fixed(gene_name, " ",3)
data2$gene_name<- gene_name2[,3]
gene_type <- column[,8]
gene_type2 <- str_split_fixed(gene_type, " ",3)
data2$gene_type<- gene_type2[,3]
data3 <- data2[,c("V1", "V4", "V5", "gene_name", "gene_type", "V7")]
data4 <- subset(data3,  gene_type!="miRNA" & gene_type!="snRNA" & gene_type!="snoRNA" & gene_type!="scaRNA"   )

#For Sequin Reference
sequin <- read.csv("/users/enovoa/boguzhan/references/mouse_with_sequins/RNAsequins.v2.2.gtf", header=FALSE, sep="\t")
sequin2<- subset(sequin, V3=="exon")
column<- str_split_fixed(sequin2$V9, ";",3)
sequin2$gene_type<- "synthetic"
gene_name <- column[,1]
gene_name2 <- str_split_fixed(gene_name, " ",4)
sequin2$gene_name<- gene_name2[,2]
sequin3 <- sequin2[,c("V1", "V4", "V5", "gene_name", "gene_type", "V7")]

#Ribosomal RNA
rRNA_bed <- read.delim("mouse_rRNA.bed", header=FALSE)

colnames(rRNA_bed) <-c("V1", "V4", "V5", "gene_name", "gene_type", "V7")




both <- rbind(data4, sequin3, rRNA_bed)



write.table(both, file="Mus_musculus.GRCm38.102_sequinv2.2_rRNA_Rest_EXON.bed", sep="\t", quote=FALSE,col.names=FALSE, row.names=FALSE)
```



# Manipulate in order to extract transcript start positions
```R
#Manipulate the GTF File
gtf <- read.csv("/users/enovoa/boguzhan/references/mouse_with_sequins/Mus_musculus.GRCm38.102.gtf", header=FALSE, sep="\t")

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


#For Sequin Reference
sequin <- read.csv("/users/enovoa/boguzhan/references/mouse_with_sequins/RNAsequins.v2.2.gtf", header=FALSE, sep="\t")
sequin_transcript <- subset(sequin , V3=="transcript")
sequin_trans_neg <- subset(sequin_transcript, V7=="-")
sequin_trans_neg2 <- sequin_trans_neg[,c("V1", "V4", "V7")]
colnames(sequin_trans_neg2) <- c("Chr", "Position1", "Strand")
sequin_trans_pos <- subset(sequin_transcript, V7=="+")
sequin_trans_pos2 <- sequin_trans_pos[,c("V1", "V5", "V7")]
colnames(sequin_trans_pos2) <- c("Chr", "Position1", "Strand")
merged_sequin <- rbind(sequin_trans_neg2,sequin_trans_pos2 )
merged_sequin$Position2 <- merged_sequin$Position1+3
merged_sequin$Position1 <- merged_sequin$Position1-3
final_sequin <- merged_sequin[,c("Chr","Position1", "Position2", "Strand")]

#For rRNA Reference
rRNA_bed <- read.delim("mouse_rRNA.bed", header=FALSE)
rRNA_bed$Chr <- rRNA_bed$V1
rRNA_bed$Position1 <- rRNA_bed$V3 - 30
rRNA_bed$Position2 <- rRNA_bed$V3
rRNA_bed$Strand <- rRNA_bed$V6
rRNA_bed2 <- rRNA_bed[,c("Chr","Position1", "Position2", "Strand")]






both <- rbind(final,final_sequin, rRNA_bed2 )


write.table(both, file="Mus_musculus.GRCm38.102_Sequin_rRNA_Transcript_Ends.bed", quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)
```

