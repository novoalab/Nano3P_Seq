##################################################
### ANALYSIS OF THE pA Selected Zebrafish Run 4hpf REP1###
###################################################
########## OGUZHAN BEGIK APRIL 2020 ###############



####################################################
## PART 3 : MAPPING     ############################
####################################################
```bash
ref=/users/enovoa/boguzhan/references/danio_rerio/rRNA_Maternal_Zygotic.fa

minimap2 -ax map-ont --MD $ref cDNA786327_4hpf.fastq | samtools view -hSb -F 3844 - >  4hpf.rRNA.sam
samtools view  -f 0x10 -bq 59 4hpf.rRNA.sam | samtools sort - 4hpf.rRNA.sorted && samtools index 4hpf.rRNA.sorted.bam




#Filter rRNA reads based on transcript starts
bedtools intersect -abam 4hpf.rRNA.sorted.bam -b Zebrafish_rRNA_Transcript_Ends.bed -wa -wb -bed > 4hpf.rRNA.overlapping.bed 

#Remove duplicate reads
awk '!seen[$4]++' 4hpf.rRNA.overlapping.bed  | cut -f4 > 4hpf.rRNA.overlapping.reads


java -jar /users/enovoa/boguzhan/Software/picard/build/libs/picard.jar FilterSamReads \
       I=4hpf.rRNA.sorted.bam \
       O=4hpf.rRNA.overlapping.bam \
       READ_LIST_FILE=4hpf.rRNA.overlapping.reads \
       FILTER=includeReadList


samtools sort 4hpf.rRNA.overlapping.bam 4hpf.rRNA.overlapping.sorted
samtools index 4hpf.rRNA.overlapping.sorted.bam





# Exclude these reads from fastq
samtools view 4hpf.rRNA.sorted.bam | cut -f1 > 4hpf.rRNA.reads

# Excluded fastq
seqkit grep --pattern-file 4hpf.rRNA.reads --invert-match cDNA786327_4hpf.fastq > 4hpf_nonrRNA.fastq



# Map these reads to Genome + rRNA
ref=/users/enovoa/boguzhan/references/danio_rerio/zebrafish_11_sequins.fa


minimap2 -ax splice -k14 -uf --MD $ref 4hpf_nonrRNA.fastq | samtools view -hSb -F 3844 - >  4hpf_nonrRNA.genome11_sequin.bam


samtools sort 4hpf_nonrRNA.genome11_sequin.bam 4hpf_nonrRNA.genome11_sequin.sorted && rm 4hpf_nonrRNA.genome11_sequin.bam && samtools index 4hpf_nonrRNA.genome11_sequin.sorted.bam



for i in *.sorted.bam;do samtools view -F 4 $i | cut -f1 | sort | uniq | wc -l;echo $i; done

```




# Intersect mapping with GTF Gene of miRNAs
```bash
# Convert BAM to BED
#Extract start and end positions

bedtools bamtobed -i 4hpf_nonrRNA.genome11_sequin.sorted.bam > 4hpf_nonrRNA.genome11_sequin.sorted.bed


# Intersect 
bedtools intersect -abam 4hpf_nonrRNA.genome11_sequin.sorted.bam -b Danio_rerio_11_miRNA_Gene.bed -wa -wb -bed -split -S > 4hpf_all.genome11_sequin_overlapping_miRNAs.bed

awk '!seen[$4]++' 4hpf_all.genome11_sequin_overlapping_miRNAs.bed | cut -f4 > 4hpf_all.genome11_sequin_overlapping_miRNAs.reads

```




#Manipulate the BED file to extract Read Starts
```R

data.input <- read.delim("4hpf_nonrRNA.genome11_sequin.sorted.bed", header=FALSE)


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

write.table(data.processed, file="4hpf_nonrRNA.genome11_sequin.readstarts.bed", quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)
```



#Intersects BED files  (BAM and GTF-Transcript)
```bash
bedtools intersect -a 4hpf_nonrRNA.genome11_sequin.readstarts.bed -b Danio_rerio_11_SmallRNA_TranscriptEnds.bed -wa -wb > 4hpf_all.genome11_sequin_smallRNAs.bed



#Remove duplicate reads
awk '!seen[$4]++' 4hpf_all.genome11_sequin_smallRNAs.bed | cut -f4 > 4hpf_all.genome11_sequin_smallRNAs.reads


### EXTRACT THE BAM FOR SMALL RNA READS
java -jar /users/enovoa/boguzhan/Software/picard/build/libs/picard.jar FilterSamReads \
       I=4hpf_nonrRNA.genome11_sequin.sorted.bam \
       O=4hpf.genome11_sequin.smallRNAs.bam\
       READ_LIST_FILE=4hpf_all.genome11_sequin_smallRNAs.reads \
       FILTER=includeReadList

samtools sort 4hpf.genome11_sequin.smallRNAs.bam 4hpf.genome11_sequin.smallRNAs.sorted
samtools index 4hpf.genome11_sequin.smallRNAs.sorted.bam



### EXTRACT THE REST TO INTERSECT AGAIN
java -jar /users/enovoa/boguzhan/Software/picard/build/libs/picard.jar FilterSamReads \
       I=4hpf_nonrRNA.genome11_sequin.sorted.bam \
       O=4hpf.genome11_sequin.restRNAs.bam\
       READ_LIST_FILE=4hpf_all.genome11_sequin_smallRNAs.reads \
       FILTER=excludeReadList



samtools sort 4hpf.genome11_sequin.restRNAs.bam 4hpf.genome11_sequin.restRNAs.sorted
samtools index 4hpf.genome11_sequin.restRNAs.sorted.bam
```


# Intersect the rest of the bam file with the rest of the GTF annotation

```bash
bedtools bamtobed -i 4hpf.genome11_sequin.restRNAs.sorted.bam > 4hpf.genome11_sequin.restRNAs.sorted.bed
```

#Manipulate the BED file
```R

data.input <- read.delim("4hpf.genome11_sequin.restRNAs.sorted.bed", header=FALSE)


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

write.table(data, file="4hpf.genome11_sequin.restRNAs.readstarts.bed", quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)
```




#Intersects BED files  (BAM and GTF-Transcript)
```bash

bedtools intersect -a 4hpf.genome11_sequin.restRNAs.readstarts.bed -b Danio_rerio_11_Sequin4_Transcript_Ends.bed -wa -wb > 4hpf.genome11_sequin.restRNAs.overlapping.bed



#Remove duplicate reads
awk '!seen[$4]++' 4hpf.genome11_sequin.restRNAs.overlapping.bed | cut -f4 > 4hpf.genome11_sequin.restRNAs.overlapping.reads



java -jar /users/enovoa/boguzhan/Software/picard/build/libs/picard.jar FilterSamReads \
       I=4hpf.genome11_sequin.restRNAs.sorted.bam \
       O=4hpf.genome11_sequin.restoverlap.bam\
       READ_LIST_FILE=4hpf.genome11_sequin.restRNAs.overlapping.reads \
       FILTER=includeReadList

samtools sort 4hpf.genome11_sequin.restoverlap.bam 4hpf.genome11_sequin.restoverlap.sorted
samtools index 4hpf.genome11_sequin.restoverlap.sorted.bam

```




# So in the end we have three types of reads
```bash
#DNA964321_all.rRNA.overlapping.reads : Contains rRNA intersected reads
#cDNA964321_all.genome38_sequin_overlapping_miRNAs.reads : Contains miRNA gene intersected reads
#cDNA964321_all.genome38_sequin_smallRNAs.reads : Contains smallRNA TRANSCRIPT START intersected reads
#cDNA964321_all.genome38_sequin.restRNAs.overlapping.reads : Contains rest RNA TRANSCRIPT START intersected reads


# We need to remove the reads from miRNAs that are overlapping with smallRNAs or restRNAs
diff 4hpf.genome11_sequin.restRNAs.overlapping.reads 4hpf_all.genome11_sequin_overlapping_miRNAs.reads |grep ">"|cut -c 3- > 4hpf.genome11_sequin.overlapping_miRNAs_RestExcluded.reads


diff 4hpf_all.genome11_sequin_smallRNAs.reads 4hpf.genome11_sequin.overlapping_miRNAs_RestExcluded.reads |grep ">"|cut -c 3- > 4hpf_all.genome11_sequin.overlapping_miRNAs_RestSmallExcluded.reads


# Now we processed the miRNA reads file, exluding overlaps with genes
#cDNA964321_all.genome38_sequin_overlapping_miRNAs_RestSmallExcluded.reads

java -jar /users/enovoa/boguzhan/Software/picard/build/libs/picard.jar FilterSamReads \
       I=4hpf_nonrRNA.genome11_sequin.sorted.bam \
       O=4hpf.genome38_sequin.miRNAFINAL.bam\
       READ_LIST_FILE=4hpf_all.genome11_sequin.overlapping_miRNAs_RestSmallExcluded.reads\
       FILTER=includeReadList

samtools sort 4hpf.genome38_sequin.miRNAFINAL.bam 4hpf.genome38_sequin.miRNAFINAL.sorted
samtools index 4hpf.genome38_sequin.miRNAFINAL.sorted.bam


# Extract BAM files from these reads
#cDNA964321_all.genome38_sequin.smallRNAs.sorted.bam
#cDNA964321_all.genome38_sequin.restoverlap.sorted.bam
#cDNA964321_all.genome38_sequin.miRNAFINAL.sorted.bam



bedtools intersect -abam 4hpf.rRNA.overlapping.sorted.bam -b Zebrafish_rRNA_Annotation.bed -wa -wb -bed -S | awk '!seen[$4]++'>  4hpf.rRNA.overlapping.FINAL.bed


bedtools intersect -abam 4hpf.genome11_sequin.restoverlap.sorted.bam -b Danio_rerio_11_Sequin4_Rest_EXON.bed -wa -wb -bed -split -S | awk '!seen[$4]++'>  4hpf.genome11_sequin.restRNAs_FINAL.bed

bedtools intersect -abam 4hpf.genome11_sequin.smallRNAs.sorted.bam -b Danio_rerio_11_SmallRNA_Gene.bed -wa -wb -bed -split -S | awk '!seen[$4]++' > 4hpf.genome11_sequin.smallRNAs_FINAL.bed

bedtools intersect -abam 4hpf.genome38_sequin.miRNAFINAL.sorted.bam -b Danio_rerio_11_miRNA_Gene.bed -wa -wb -bed -split -S | awk '!seen[$4]++' > 4hpf.genome11_sequin.miRNAs_FINAL.bed

cat  4hpf.rRNA.overlapping.FINAL.bed 4hpf.genome11_sequin.restRNAs_FINAL.bed 4hpf.genome11_sequin.smallRNAs_FINAL.bed 4hpf.genome11_sequin.miRNAs_FINAL.bed > 4hpf.genome11_sequin_ALLRNAs_Merged.bed

```

