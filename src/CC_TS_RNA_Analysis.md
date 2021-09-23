##################################################
### ANALYSIS OF THE CURLCAKE TS RUN ###
###################################################
########## OGUZHAN BEGIK APRIL 2020 ###############



```bash
#Extract start and end positions
bedtools bamtobed -i cDNA345234.sorted.bam > TS.bed
bedtools bamtobed -i  cDNA345234_porechop.sorted.bam > TS_trimmed.bed


```


#Manipulate the bed file containing only one strand 
```R

TS <- read.delim("TS.bed", header=FALSE)


manipulate <- function(data) {
	subs_neg <- subset(data, V6 =="-")
	subs_neg2 <- subs_neg[,c("V1", "V3", "V6", "V4")]
	colnames(subs_neg2) <- c("Chr", "Position", "Strand", "Read_ID")
	subs_neg2$Position1 <- subs_neg2$Position-5
	subs_neg2$Position2 <-  subs_neg2$Position+5
	final <- subs_neg2[,c("Chr","Position1", "Position2", "Strand", "Read_ID")]
	return(final)
}

ts.processed <- manipulate(TS)

write.table(ts.processed, file="TS_processed.bed", quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)
```

#Intersects BED files 
```bash
bedtools intersect -a TS_processed.bed -b curlcake12.bed -wa -wb > TS_overlapping.bed
bedtools intersect -a TS_trimmed_processed.bed -b curlcake12.bed -wa -wb > TS_trimmed_overlapping.bed


#Remove duplicate reads
awk '!seen[$5]++' TS_overlapping.bed > TS_processed_unique.bed
awk '!seen[$5]++' TS_trimmed_overlapping.bed > TS_trimmed_processed_unique.bed

#Extract the Read IDs 
cut -f5 TS_processed_unique.bed > TS_processed_unique.readid
cut -f5 TS_trimmed_processed_unique.bed > TS_trimmed_processed_unique.readid

```


#Extract the BAM files overlapping with read ends ONLY
```bash
java -jar /users/enovoa/boguzhan/Software/picard/build/libs/picard.jar FilterSamReads \
       I= cDNA345234.sorted.bam \
       O= TS.complete.bam\
       READ_LIST_FILE=TS_processed_unique.readid \
       FILTER=includeReadList

samtools index TS.complete.bam


samtools view TS.complete.bam | cut -f1,3,5 > TS.read_id


java -jar /users/enovoa/boguzhan/Software/picard/build/libs/picard.jar FilterSamReads \
       I= cDNA345234_porechop.sorted.bam \
       O= TS_trimmed.complete.bam\
       READ_LIST_FILE=TS_trimmed_processed_unique.readid \
       FILTER=includeReadList

samtools index TS_trimmed.complete.bam


samtools view TS_trimmed.complete.bam | cut -f1,3,5 > TS_trimmed.read_id



```


