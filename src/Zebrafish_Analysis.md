##################################################
### NALYSIS OF THE pA Selected Zebrafish Run ###
###################################################
########## OGUZHAN BEGIK APRIL 2020 ###############



####################################################
## PART 3 : MAPPING     ############################
####################################################
```bash
ref=/users/enovoa/boguzhan/references/danio_rerio/rRNA_Maternal_Zygotic.fa


minimap2 -ax map-ont --MD $ref cDNA8523612.fastq | samtools view -hSb -F 3844 - >  cDNA8523612.rRNA.sam

samtools view  -f 0x10 -bq 59 cDNA8523612.rRNA.sam | samtools sort - cDNA8523612.rRNA.sorted && samtools index cDNA8523612.rRNA.sorted.bam




#Filter rRNA reads based on transcript starts
bedtools intersect -abam cDNA8523612.rRNA.sorted.bam -b Zebrafish_rRNA_Transcript_Ends.bed -wa -wb -bed > cDNA8523612.rRNA.overlapping.bed 

#Remove duplicate reads
awk '!seen[$4]++' cDNA8523612.rRNA.overlapping.bed  | cut -f4 > cDNA8523612.rRNA.overlapping.reads


java -jar /users/enovoa/boguzhan/Software/picard/build/libs/picard.jar FilterSamReads \
       I=cDNA8523612.rRNA.sorted.bam \
       O=cDNA8523612.rRNA.overlapping.bam \
       READ_LIST_FILE=cDNA8523612.rRNA.overlapping.reads \
       FILTER=includeReadList


samtools sort cDNA8523612.rRNA.overlapping.bam cDNA8523612.rRNA.overlapping.sorted
samtools index cDNA8523612.rRNA.overlapping.sorted.bam





# Exclude these reads from fastq
samtools view cDNA8523612.rRNA.sorted.bam | cut -f1 > cDNA8523612.rRNA.reads

# Excluded fastq
seqkit grep --pattern-file cDNA8523612.rRNA.reads --invert-match cDNA8523612.fastq > cDNA8523612_nonrRNA.fastq



# Map these reads to Genome + rRNA
ref=/users/enovoa/boguzhan/references/danio_rerio/zebrafish_11_sequins.fa


minimap2 -ax splice -k14 -uf --MD $ref cDNA8523612_nonrRNA.fastq | samtools view -hSb -F 3844 - >  cDNA8523612_nonrRNA.genome11_sequin.bam


samtools sort cDNA8523612_nonrRNA.genome11_sequin.bam cDNA8523612_nonrRNA.genome11_sequin.sorted && rm cDNA8523612_nonrRNA.genome11_sequin.bam && samtools index cDNA8523612_nonrRNA.genome11_sequin.sorted.bam



for i in *.sorted.bam;do samtools view -F 4 $i | cut -f1 | sort | uniq | wc -l;echo $i; done

```




# Intersect mapping with GTF Gene of miRNAs
```bash
# Convert BAM to BED
#Extract start and end positions

bedtools bamtobed -i cDNA8523612_nonrRNA.genome11_sequin.sorted.bam > cDNA8523612_nonrRNA.genome11_sequin.sorted.bed


# Intersect 

bedtools intersect -abam cDNA8523612_nonrRNA.genome11_sequin.sorted.bam -b Danio_rerio_11_miRNA_Gene.bed -wa -wb -bed -split -S > cDNA8523612_all.genome11_sequin_overlapping_miRNAs.bed

awk '!seen[$4]++' cDNA8523612_all.genome11_sequin_overlapping_miRNAs.bed | cut -f4 > cDNA8523612_all.genome11_sequin_overlapping_miRNAs.reads

```




#Manipulate the BED file to extract Read Starts
```R

data.input <- read.delim("cDNA8523612_nonrRNA.genome11_sequin.sorted.bed", header=FALSE)


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

write.table(data.processed, file="cDNA8523612_nonrRNA.genome11_sequin.readstarts.bed", quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)
```



#Intersects BED files  (BAM and GTF-Transcript)
```bash
bedtools intersect -a cDNA8523612_nonrRNA.genome11_sequin.readstarts.bed -b Danio_rerio_11_SmallRNA_TranscriptEnds.bed -wa -wb > cDNA8523612_all.genome11_sequin_smallRNAs.bed



#Remove duplicate reads
awk '!seen[$4]++' cDNA8523612_all.genome11_sequin_smallRNAs.bed | cut -f4 > cDNA8523612_all.genome11_sequin_smallRNAs.reads


### EXTRACT THE BAM FOR SMALL RNA READS
java -jar /users/enovoa/boguzhan/Software/picard/build/libs/picard.jar FilterSamReads \
       I=cDNA8523612_nonrRNA.genome11_sequin.sorted.bam \
       O=cDNA8523612.genome11_sequin.smallRNAs.bam\
       READ_LIST_FILE=cDNA8523612_all.genome11_sequin_smallRNAs.reads \
       FILTER=includeReadList

samtools sort cDNA8523612.genome11_sequin.smallRNAs.bam cDNA8523612.genome11_sequin.smallRNAs.sorted
samtools index cDNA8523612.genome11_sequin.smallRNAs.sorted.bam



### EXTRACT THE REST TO INTERSECT AGAIN
java -jar /users/enovoa/boguzhan/Software/picard/build/libs/picard.jar FilterSamReads \
       I=cDNA8523612_nonrRNA.genome11_sequin.sorted.bam \
       O=cDNA8523612.genome11_sequin.restRNAs.bam\
       READ_LIST_FILE=cDNA8523612_all.genome11_sequin_smallRNAs.reads \
       FILTER=excludeReadList



samtools sort cDNA8523612.genome11_sequin.restRNAs.bam cDNA8523612.genome11_sequin.restRNAs.sorted
samtools index cDNA8523612.genome11_sequin.restRNAs.sorted.bam
```


# Intersect the rest of the bam file with the rest of the GTF annotation

```bash
bedtools bamtobed -i cDNA8523612.genome11_sequin.restRNAs.sorted.bam > cDNA8523612.genome11_sequin.restRNAs.sorted.bed
```

#Manipulate the BED file
```R

data.input <- read.delim("cDNA8523612.genome11_sequin.restRNAs.sorted.bed", header=FALSE)


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

write.table(data, file="cDNA8523612.genome11_sequin.restRNAs.readstarts.bed", quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)
```




#Intersects BED files  (BAM and GTF-Transcript)
```bash

bedtools intersect -a cDNA8523612.genome11_sequin.restRNAs.readstarts.bed -b Danio_rerio_11_Sequin4_Transcript_Ends.bed -wa -wb > cDNA8523612.genome11_sequin.restRNAs.overlapping.bed



#Remove duplicate reads
awk '!seen[$4]++' cDNA8523612.genome11_sequin.restRNAs.overlapping.bed | cut -f4 > cDNA8523612.genome11_sequin.restRNAs.overlapping.reads



java -jar /users/enovoa/boguzhan/Software/picard/build/libs/picard.jar FilterSamReads \
       I=cDNA8523612.genome11_sequin.restRNAs.sorted.bam \
       O=cDNA8523612.genome11_sequin.restoverlap.bam\
       READ_LIST_FILE=cDNA8523612.genome11_sequin.restRNAs.overlapping.reads \
       FILTER=includeReadList

samtools sort cDNA8523612.genome11_sequin.restoverlap.bam cDNA8523612.genome11_sequin.restoverlap.sorted
samtools index cDNA8523612.genome11_sequin.restoverlap.sorted.bam

```




# So in the end we have three types of reads
```bash
#DNA964321_all.rRNA.overlapping.reads : Contains rRNA intersected reads
#cDNA964321_all.genome38_sequin_overlapping_miRNAs.reads : Contains miRNA gene intersected reads
#cDNA964321_all.genome38_sequin_smallRNAs.reads : Contains smallRNA TRANSCRIPT START intersected reads
#cDNA964321_all.genome38_sequin.restRNAs.overlapping.reads : Contains rest RNA TRANSCRIPT START intersected reads


# We need to remove the reads from miRNAs that are overlapping with smallRNAs or restRNAs
diff cDNA8523612.genome11_sequin.restRNAs.overlapping.reads cDNA8523612_all.genome11_sequin_overlapping_miRNAs.reads |grep ">"|cut -c 3- > cDNA8523612.genome11_sequin.overlapping_miRNAs_RestExcluded.reads


diff cDNA8523612_all.genome11_sequin_smallRNAs.reads cDNA8523612.genome11_sequin.overlapping_miRNAs_RestExcluded.reads |grep ">"|cut -c 3- > cDNA8523612_all.genome11_sequin.overlapping_miRNAs_RestSmallExcluded.reads


# Now we processed the miRNA reads file, exluding overlaps with genes
#cDNA964321_all.genome38_sequin_overlapping_miRNAs_RestSmallExcluded.reads

java -jar /users/enovoa/boguzhan/Software/picard/build/libs/picard.jar FilterSamReads \
       I=cDNA8523612_nonrRNA.genome11_sequin.sorted.bam \
       O=cDNA8523612.genome38_sequin.miRNAFINAL.bam\
       READ_LIST_FILE=cDNA8523612_all.genome11_sequin.overlapping_miRNAs_RestSmallExcluded.reads\
       FILTER=includeReadList

samtools sort cDNA8523612.genome38_sequin.miRNAFINAL.bam cDNA8523612.genome38_sequin.miRNAFINAL.sorted
samtools index cDNA8523612.genome38_sequin.miRNAFINAL.sorted.bam


# Extract BAM files from these reads
#cDNA964321_all.genome38_sequin.smallRNAs.sorted.bam
#cDNA964321_all.genome38_sequin.restoverlap.sorted.bam
#cDNA964321_all.genome38_sequin.miRNAFINAL.sorted.bam



bedtools intersect -abam cDNA8523612.rRNA.overlapping.sorted.bam -b Zebrafish_rRNA_Annotation.bed -wa -wb -bed -S | awk '!seen[$4]++'>  cDNA8523612.rRNA.overlapping.FINAL.bed


bedtools intersect -abam cDNA8523612.genome11_sequin.restoverlap.sorted.bam -b Danio_rerio_11_Sequin4_Rest_EXON.bed -wa -wb -bed -split -S | awk '!seen[$4]++'>  cDNA8523612.genome11_sequin.restRNAs_FINAL.bed

bedtools intersect -abam cDNA8523612.genome11_sequin.smallRNAs.sorted.bam -b Danio_rerio_11_SmallRNA_Gene.bed -wa -wb -bed -split -S | awk '!seen[$4]++' > cDNA8523612.genome11_sequin.smallRNAs_FINAL.bed

bedtools intersect -abam cDNA8523612.genome38_sequin.miRNAFINAL.sorted.bam -b Danio_rerio_11_miRNA_Gene.bed -wa -wb -bed -split -S | awk '!seen[$4]++' > cDNA8523612.genome11_sequin.miRNAs_FINAL.bed

cat  cDNA8523612.rRNA.overlapping.FINAL.bed cDNA8523612.genome11_sequin.restRNAs_FINAL.bed cDNA8523612.genome11_sequin.smallRNAs_FINAL.bed cDNA8523612.genome11_sequin.miRNAs_FINAL.bed > cDNA8523612.genome11_sequin_ALLRNAs_Merged.bed

```


