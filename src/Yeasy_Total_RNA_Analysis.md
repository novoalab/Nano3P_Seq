##################################################
### ANALYSIS OF THE Yeast  Total RNA RUN ###
###################################################
########## OGUZHAN BEGIK APRIL 2020 ###############



#Extract start and end positions
```bash
bedtools bamtobed -i cDNA231532_SacCer_Genome_rRNA.sorted.bam > cDNA231532_SacCer_Genome_rRNA.bed
```

#Manipulate the bed file containing only one strand 
```R
data.input <- read.delim("cDNA231532_SacCer_Genome_rRNA.bed", header=FALSE)


manipulate <- function(data) {
	subs_neg <- subset(data, V6 =="-")
	subs_neg2 <- subs_neg[,c("V1", "V3", "V6", "V4")]
	colnames(subs_neg2) <- c("Chr", "Position", "Strand", "Read_ID")
	subs_pos <- subset(data, V6 =="+")
	subs_pos2 <- subs_pos[,c("V1", "V2", "V6", "V4")]
	colnames(subs_pos2) <- c("Chr", "Position", "Strand", "Read_ID")
	merged <- rbind(subs_neg2,subs_pos2 )
	merged$Position1 <- merged$Position-10
	merged$Position2 <-  merged$Position+10
	final <- merged[,c("Chr","Position1", "Position2",  "Read_ID","Strand")]
	final2 <- subset(final, Position1 > 0 & Position2 > 0)
	return(final2)
}

data.processed <- manipulate(data.input)

write.table(data.processed, file="cDNA231532_SacCer_Genome_rRNA.readstarts.bed", quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)
```


#Intersects BED files 
```bash
bedtools intersect -a cDNA231532_SacCer_Genome_rRNA.readstarts.bed -b saccer_cds_rRNA_transcript_ends.bed -wa -wb > cDNA231532_SacCer_Genome_rRNA_overlapping.bed 


#Remove duplicate reads
awk '!seen[$4]++' cDNA231532_SacCer_Genome_rRNA_overlapping.bed > cDNA231532_SacCer_Genome_rRNA_unique.bed

#Extract the Read IDs 
cut -f4 cDNA231532_SacCer_Genome_rRNA_unique.bed > cDNA231532_SacCer_Genome_rRNA_unique.readid

```

#Manipulate the GTF of SacCer3
```R
library(stringr)

#For Mouse Reference
data <- read.csv("Saccer64.gtf", header=FALSE, sep="\t")
data2<- subset(data, V3=="exon")
column<- str_split_fixed(data2$V9, ";",10)
gene_name <- column[,1]
gene_name2 <- str_split_fixed(gene_name, " ",2)
data2$gene_name<- gene_name2[,2]
gene_type <- column[,2]
gene_type2 <- str_split_fixed(gene_type, "_",3)
data2$gene_type<- gene_type2[,3]
data2$V1 <- paste("chr", data2$V1, sep="")
data3 <- data2[,c("V1", "V4", "V5", "gene_name", "gene_type", "V7")]

write.table(data3, file="Saccer64_EXON.bed", sep="\t", quote=FALSE,col.names=FALSE, row.names=FALSE)
```




#Extract the BAM files overlapping with read ends ONLY
```bash
java -jar /users/enovoa/boguzhan/Software/picard/build/libs/picard.jar FilterSamReads \
       I= cDNA231532_SacCer_Genome_rRNA.sorted.bam \
       O= cDNA231532_SacCer_Genome_rRNA.complete.bam\
       READ_LIST_FILE=cDNA231532_SacCer_Genome_rRNA_unique.readid \
       FILTER=includeReadList

samtools index cDNA231532_SacCer_Genome_rRNA.complete.bam


bedtools intersect -abam cDNA231532_SacCer_Genome_rRNA.complete.bam -b Saccer64_Exon_rRNA.bed -wa -wb -bed -split -S | awk '!seen[$4]++'>  cDNA231532_SacCer_Genome_rRNA_FINAL.bed

```