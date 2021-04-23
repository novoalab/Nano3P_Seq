##################################################
### ANALYSIS OF PRE-RNA in ALL SPECIES ###
###################################################
########## OGUZHAN BEGIK APRIL 2020 ###############



###### MOUSE 

## Map to pre-rRNA MOUSE

```bash

ref=mouse_pre_rRNA.fasta

minimap2 -ax map-ont $ref Mouse.fastq -o Mouse_preRNA.sam --MD

samtools view  -Sb -f 0x10 Mouse_preRNA.sam | samtools sort - Mouse_preRNA.sorted && samtools index Mouse_preRNA.sorted.bam


samtools view Mouse_preRNA.sorted.bam | cut -f1,3,5 > Mouse_preRNA.read_id

#Intersect
bedtools bamtobed -i Mouse_preRNA.sorted.bam> Mouse_preRNA.bed

```

#Manipulate the BED file
```R

data <- read.delim("Mouse_preRNA.bed", header=FALSE)


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

data.processed <- manipulate(data)

write.table(data.processed, file="Mouse_preRNA.processed.bed", quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)
```




#Intersects BED files  (BAM and GTF-Transcript)
```bash
bedtools intersect -a Mouse_preRNA.processed.bed -b mouse_pre_rna_18s_processed.bed  -wa -wb > Mouse_18s_Processed.bed
bedtools intersect -a Mouse_preRNA.processed.bed -b mouse_pre_rRNA_18s_pre.bed  -wa -wb >  Mouse_18s_Precursor.bed

#Remove duplicate reads
awk '!seen[$4]++' Mouse_18s_Processed.bed| cut -f4 > Mouse_18s_Processed.reads
awk '!seen[$4]++' Mouse_18s_Precursor.bed| cut -f4 > Mouse_18s_Precursor.reads

### EXTRACT THE BAM FOR SMALL RNA READS

java -jar /users/enovoa/boguzhan/Software/picard/build/libs/picard.jar FilterSamReads \
       I=Mouse_preRNA.sorted.bam \
       O=Mouse_18s_Processed.bam\
       READ_LIST_FILE=Mouse_18s_Processed.reads \
       FILTER=includeReadList

samtools index Mouse_18s_Processed.bam


java -jar /users/enovoa/boguzhan/Software/picard/build/libs/picard.jar FilterSamReads \
       I=Mouse_preRNA.sorted.bam \
       O=Mouse_18s_Precursor.bam\
       READ_LIST_FILE=Mouse_18s_Precursor.reads \
       FILTER=includeReadList

samtools index Mouse_18s_Precursor.bam




samtools view processed_18s_rRNA.bam | cut -f1,3,5 > processed_18s_rRNA.read_id
samtools view pre_18s_rRNA.bam | cut -f1,3,5 > pre_18s_rRNA.read_id

```


#Intersects BED files  (BAM and GTF-Transcript)
```bash
bedtools intersect -a Mouse_preRNA.processed.bed -b mouse_pre_rna_28s_processed.bed  -wa -wb > Mouse_28s_Processed.bed
bedtools intersect -a Mouse_preRNA.processed.bed -b mouse_pre_rRNA_28s_pre.bed  -wa -wb >  Mouse_28s_Precursor.bed

#Remove duplicate reads
awk '!seen[$4]++' Mouse_28s_Processed.bed| cut -f4 > Mouse_28s_Processed.reads
awk '!seen[$4]++' Mouse_28s_Precursor.bed| cut -f4 > Mouse_28s_Precursor.reads

### EXTRACT THE BAM FOR SMALL RNA READS

java -jar /users/enovoa/boguzhan/Software/picard/build/libs/picard.jar FilterSamReads \
       I=Mouse_preRNA.sorted.bam \
       O=Mouse_28s_Processed.bam\
       READ_LIST_FILE=Mouse_28s_Processed.reads \
       FILTER=includeReadList

samtools index Mouse_28s_Processed.bam


java -jar /users/enovoa/boguzhan/Software/picard/build/libs/picard.jar FilterSamReads \
       I=Mouse_preRNA.sorted.bam \
       O=Mouse_28s_Precursor.bam\
       READ_LIST_FILE=Mouse_28s_Precursor.reads \
       FILTER=includeReadList

samtools index Mouse_28s_Precursor.bam
```








###### YEAST 

## Map to pre-rRNA YEAST

```bash

ref=yeast_pre_RNA.fasta

minimap2 -ax map-ont $ref Yeast.fastq -o Yeast_preRNA.sam --MD

samtools view  -Sb -f 0x10 Yeast_preRNA.sam | samtools sort - Yeast_preRNA.sorted && samtools index Yeast_preRNA.sorted.bam


samtools view Yeast_preRNA.sorted.bam | cut -f1,3,5 > Yeast_preRNA.read_id

#Intersect
bedtools bamtobed -i Yeast_preRNA.sorted.bam> Yeast_preRNA.bed

```

#Manipulate the BED file
```R

data <- read.delim("Yeast_preRNA.bed", header=FALSE)


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

data.processed <- manipulate(data)

write.table(data.processed, file="Yeast_preRNA.processed.bed", quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)
```




#Intersects BED files  (BAM and GTF-Transcript)
```bash
bedtools intersect -a Yeast_preRNA.processed.bed -b yeast_pre_rRNA_18s_processed.bed  -wa -wb > Yeast_18s_Processed.bed
bedtools intersect -a Yeast_preRNA.processed.bed -b yeast_pre_rRNA_18s_precursor.bed -wa -wb >  Yeast_18s_Precursor.bed

#Remove duplicate reads
awk '!seen[$4]++' Yeast_18s_Processed.bed| cut -f4 > Yeast_18s_Processed.reads
awk '!seen[$4]++' Yeast_18s_Precursor.bed| cut -f4 > Yeast_18s_Precursor.reads

### EXTRACT THE BAM FOR SMALL RNA READS

java -jar /users/enovoa/boguzhan/Software/picard/build/libs/picard.jar FilterSamReads \
       I=Yeast_preRNA.sorted.bam \
       O=Yeast_18s_Processed.bam\
       READ_LIST_FILE=Yeast_18s_Processed.reads \
       FILTER=includeReadList

samtools index Yeast_18s_Processed.bam


java -jar /users/enovoa/boguzhan/Software/picard/build/libs/picard.jar FilterSamReads \
       I=Yeast_preRNA.sorted.bam \
       O=Yeast_18s_Precursor.bam\
       READ_LIST_FILE=Yeast_18s_Precursor.reads \
       FILTER=includeReadList

samtools index Yeast_18s_Precursor.bam




samtools view processed_18s_rRNA.bam | cut -f1,3,5 > processed_18s_rRNA.read_id
samtools view pre_18s_rRNA.bam | cut -f1,3,5 > pre_18s_rRNA.read_id

```