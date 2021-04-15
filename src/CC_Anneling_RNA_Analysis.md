##################################################
### ANALYSIS OF THE CURLCAKE ANNEALING RUN ###
###################################################
########## OGUZHAN BEGIK APRIL 2020 ###############



```bash
#Extract start and end positions
bedtools bamtobed -i cDNA325423.sorted.bam > Annealing.bed
```


#Manipulate the bed file containing only one strand 
```R

Annealing <- read.delim("Annealing.bed", header=FALSE)


manipulate <- function(data) {
	subs_neg <- subset(data, V6 =="-")
	subs_neg2 <- subs_neg[,c("V1", "V3", "V6", "V4")]
	colnames(subs_neg2) <- c("Chr", "Position", "Strand", "Read_ID")
	subs_neg2$Position1 <- subs_neg2$Position-5
	subs_neg2$Position2 <-  subs_neg2$Position+5
	final <- subs_neg2[,c("Chr","Position1", "Position2", "Strand", "Read_ID")]
	return(final)
}

ann.processed <- manipulate(Annealing)

write.table(ann.processed, file="Annealing_processed.bed", quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)
```

#Intersects BED files 
```bash
bedtools intersect -a Annealing_processed.bed -b curlcake12.bed -wa -wb > Annealing_overlapping.bed


#Remove duplicate reads
awk '!seen[$5]++' Annealing_overlapping.bed > Annealing_processed_unique.bed

#Extract the Read IDs 
cut -f5 Annealing_processed_unique.bed > Annealing_processed_unique.readid

```


#Extract the BAM files overlapping with read ends ONLY
```bash
java -jar /users/enovoa/boguzhan/Software/picard/build/libs/picard.jar FilterSamReads \
       I= cDNA325423.sorted.bam \
       O= Annealing.complete.bam\
       READ_LIST_FILE=Annealing_processed_unique.readid \
       FILTER=includeReadList

samtools index Annealing.complete.bam


samtools view Annealing.complete.bam | cut -f1,3,5 > Annealing.read_id
```


