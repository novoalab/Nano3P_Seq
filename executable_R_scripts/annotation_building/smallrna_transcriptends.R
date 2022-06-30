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


#Export the final table
write.table(final, file=paste(label, "smallrna_transcriptends.bed", sep="_"), quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)