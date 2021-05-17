###################################################
###### ANALYSIS OF THE NUCLEAR RNA RUN ###########
###################################################
########## OGUZHAN BEGIK APRIL 2020 ###############


# Intersect mapping with GTF Gene of miRNAs
```bash
# Convert BAM to BED
#Extract start and end positions
bedtools bamtobed -i cDNA964321_all_genome38_sequin_rrna.sorted.bam > cDNA964321_all_genome38_sequin_rrna.sorted.bed


# Intersect 

bedtools intersect -abam cDNA964321_all_genome38_sequin_rrna.sorted.bam -b Mus_musculus.GRCm38.102_sequinv2.2_miRNA_Gene.bed -wa -wb -bed -split -S > cDNA964321_all_genome38_sequin_rrna_overlapping_miRNAs.bed

awk '!seen[$4]++' cDNA964321_all_genome38_sequin_rrna_overlapping_miRNAs.bed | cut -f4 > cDNA964321_all_genome38_sequin_rrna_overlapping_miRNAs.reads


```

#Manipulate the BED file to extract Read Starts
```R

data.input <- read.delim("cDNA964321_all_genome38_sequin_rrna.sorted.bed", header=FALSE)


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

write.table(data.processed, file="cDNA964321_all_genome38_sequin_rrna.readstarts.bed", quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)
```



#Intersects BED files  (BAM and GTF-Transcript)
```bash
bedtools intersect -a cDNA964321_all_genome38_sequin_rrna.readstarts.bed -b Mus_musculus.GRCm38.102_sequinv2.2_SmallRNA_TranscriptEnds.bed -wa -wb > cDNA964321_all.genome38_sequin_smallRNAs.bed


#Remove duplicate reads
awk '!seen[$4]++' cDNA964321_all.genome38_sequin_smallRNAs.bed | cut -f4 > cDNA964321_all.genome38_sequin_smallRNAs.reads


### EXTRACT THE BAM FOR SMALL RNA READS
java -jar /users/enovoa/boguzhan/Software/picard/build/libs/picard.jar FilterSamReads \
       I=cDNA964321_all_genome38_sequin_rrna.sorted.bam \
       O=cDNA964321_all.genome38_sequin_rrna.smallRNAs.bam\
       READ_LIST_FILE=cDNA964321_all.genome38_sequin_smallRNAs.reads \
       FILTER=includeReadList

samtools sort cDNA964321_all.genome38_sequin_rrna.smallRNAs.bam cDNA964321_all.genome38_sequin_rrna.smallRNAs.sorted
samtools index cDNA964321_all.genome38_sequin_rrna.smallRNAs.sorted.bam



### EXTRACT THE REST TO INTERSECT AGAIN
java -jar /users/enovoa/boguzhan/Software/picard/build/libs/picard.jar FilterSamReads \
       I=cDNA964321_all_genome38_sequin_rrna.sorted.bam \
       O=cDNA964321_all.genome38_sequin_rrna.restRNAs.bam\
       READ_LIST_FILE=cDNA964321_all.genome38_sequin_smallRNAs.reads \
       FILTER=excludeReadList


samtools sort cDNA964321_all.genome38_sequin_rrna.restRNAs.bam cDNA964321_all.genome38_sequin_rrna.restRNAs.sorted
samtools index cDNA964321_all.genome38_sequin_rrna.restRNAs.sorted.bam
```



# Intersect the rest of the bam file with the rest of the GTF annotation

```bash
bedtools bamtobed -i cDNA964321_all.genome38_sequin_rrna.restRNAs.sorted.bam > cDNA964321_all.genome38_sequin_rrna.restRNAs.sorted.bed
```


#Manipulate the BED file
```R

data.input <- read.delim("cDNA964321_all.genome38_sequin_rrna.restRNAs.sorted.bed", header=FALSE)


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

data <- manipulate(data.input)

write.table(data, file="cDNA964321_all.genome38_sequin_rRNA.restRNAs.readstarts.bed", quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)
```


#Intersects BED files  (BAM and GTF-Transcript)
```bash

bedtools intersect -a cDNA964321_all.genome38_sequin_rRNA.restRNAs.readstarts.bed -b Mus_musculus.GRCm38.102_Sequin_rRNA_Transcript_Ends.bed -wa -wb > cDNA964321_all.genome38_sequin_rRNA.restRNAs.overlapping.bed



#Remove duplicate reads
awk '!seen[$4]++' cDNA964321_all.genome38_sequin_rRNA.restRNAs.overlapping.bed | cut -f4 > cDNA964321_all.genome38_sequin_rRNA.restRNAs.overlapping.reads



java -jar /users/enovoa/boguzhan/Software/picard/build/libs/picard.jar FilterSamReads \
       I=cDNA964321_all_genome38_sequin_rrna.sorted.bam  \
       O=cDNA964321_all.genome38_sequin_rRNA.restoverlap.bam\
       READ_LIST_FILE=cDNA964321_all.genome38_sequin_rRNA.restRNAs.overlapping.reads \
       FILTER=includeReadList

samtools sort cDNA964321_all.genome38_sequin_rRNA.restoverlap.bam cDNA964321_all.genome38_sequin_rRNA.restoverlap.sorted
samtools index cDNA964321_all.genome38_sequin_rRNA.restoverlap.sorted.bam


bedtools bamtobed -i cDNA964321_all.genome38_sequin_rRNA.restoverlap.sorted.bam > cDNA964321_all.genome38_sequin_rRNA.restoverlap.sorted.bed

```



# So in the end we have three types of reads
```bash
#DNA964321_all.rRNA.overlapping.reads : Contains rRNA intersected reads
#cDNA964321_all.genome38_sequin_overlapping_miRNAs.reads : Contains miRNA gene intersected reads
#cDNA964321_all.genome38_sequin_smallRNAs.reads : Contains smallRNA TRANSCRIPT START intersected reads
#cDNA964321_all.genome38_sequin.restRNAs.overlapping.reads : Contains rest RNA TRANSCRIPT START intersected reads


# We need to remove the reads from miRNAs that are overlapping with smallRNAs or restRNAs
diff cDNA964321_all.genome38_sequin_rRNA.restRNAs.overlapping.reads cDNA964321_all_genome38_sequin_rrna_overlapping_miRNAs.reads |grep ">"|cut -c 3- > cDNA964321_all.genome38_sequin_rRNA_overlapping_miRNAs_RestExcluded.reads


diff cDNA964321_all.genome38_sequin_smallRNAs.reads cDNA964321_all_genome38_sequin_rrna_overlapping_miRNAs.reads |grep ">"|cut -c 3- > cDNA964321_all.genome38_sequin_rRNAoverlapping_miRNAs_RestSmallExcluded.reads


# Now we processed the miRNA reads file, exluding overlaps with genes
#cDNA964321_all.genome38_sequin_overlapping_miRNAs_RestSmallExcluded.reads

java -jar /users/enovoa/boguzhan/Software/picard/build/libs/picard.jar FilterSamReads \
       I=cDNA964321_all_genome38_sequin_rrna.sorted.bam \
       O=cDNA964321_all.genome38_sequin.miRNAFINAL.bam\
       READ_LIST_FILE=cDNA964321_all.genome38_sequin_rRNAoverlapping_miRNAs_RestSmallExcluded.reads\
       FILTER=includeReadList

samtools sort cDNA964321_all.genome38_sequin.miRNAFINAL.bam cDNA964321_all.genome38_sequin.miRNAFINAL.sorted
samtools index cDNA964321_all.genome38_sequin.miRNAFINAL.sorted.bam


# Extract BAM files from these reads
#cDNA964321_all.genome38_sequin.smallRNAs.sorted.bam
#cDNA964321_all.genome38_sequin.restoverlap.sorted.bam
#cDNA964321_all.genome38_sequin.miRNAFINAL.sorted.bam



bedtools intersect -abam cDNA964321_all.genome38_sequin_rRNA.restoverlap.sorted.bam -b Mus_musculus.GRCm38.102_sequinv2.2_rRNA_Rest_EXON.bed -wa -wb -bed -split -S | awk '!seen[$4]++'>  cDNA964321_all.genome38_sequin_rRNA.restRNAs_FINAL.bed

bedtools intersect -abam cDNA964321_all.genome38_sequin_rrna.smallRNAs.sorted.bam -b Mus_musculus.GRCm38.102_sequinv2.2_SmallRNA_Gene.bed -wa -wb -bed -split -S | awk '!seen[$4]++' > cDNA964321_all.genome38_sequin_rRNA.smallRNAs_FINAL.bed

bedtools intersect -abam cDNA964321_all.genome38_sequin.miRNAFINAL.sorted.bam -b Mus_musculus.GRCm38.102_sequinv2.2_miRNA_Gene.bed -wa -wb -bed -split -S | awk '!seen[$4]++' > cDNA964321_all.genome38_sequin_rRNA.miRNAs_FINAL.bed

cat  cDNA964321_all.genome38_sequin_rRNA.restRNAs_FINAL.bed cDNA964321_all.genome38_sequin_rRNA.smallRNAs_FINAL.bed cDNA964321_all.genome38_sequin_rRNA.miRNAs_FINAL.bed > cDNA964321_all.genome38_sequin_rRNA_ALLRNAs_Merged.bed
```


