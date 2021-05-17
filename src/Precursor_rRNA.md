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

```




#Intersects BED files  (BAM and GTF-Transcript)
```bash
bedtools intersect -a Yeast_preRNA.processed.bed -b yeast_pre_rna_28s_processed.bed  -wa -wb > Yeast_28s_Processed.bed
bedtools intersect -a Yeast_preRNA.processed.bed -b yeast_pre_rRNA_28s_pre.bed  -wa -wb >  Yeast_28s_Precursor.bed
#Remove duplicate reads
awk '!seen[$4]++' Yeast_28s_Processed.bed| cut -f4 > Yeast_28s_Processed.reads
awk '!seen[$4]++' Yeast_28s_Precursor.bed| cut -f4 > Yeast_28s_Precursor.reads

### EXTRACT THE BAM FOR SMALL RNA READS

java -jar /users/enovoa/boguzhan/Software/picard/build/libs/picard.jar FilterSamReads \
       I=Yeast_preRNA.sorted.bam \
       O=Yeast_28s_Processed.bam\
       READ_LIST_FILE=Yeast_28s_Processed.reads \
       FILTER=includeReadList

samtools index Yeast_28s_Processed.bam


java -jar /users/enovoa/boguzhan/Software/picard/build/libs/picard.jar FilterSamReads \
       I=Yeast_preRNA.sorted.bam \
       O=Yeast_28s_Precursor.bam\
       READ_LIST_FILE=Yeast_28s_Precursor.reads \
       FILTER=includeReadList

samtools index Yeast_28s_Precursor.bam
```








# Running Epinano on the bam files
```bash
#Epinano.sh
java -jar /users/enovoa/boguzhan/Software/picard/build/libs/picard.jar CreateSequenceDictionary REFERENCE=mouse_pre_rRNA.fasta OUTPUT= mouse_pre_rRNA.fasta.dict


ref=mouse_pre_rRNA.fasta
sam2tsv=/users/enovoa/boguzhan/Software/jvarkit/dist/sam2tsv.jar
tsv_to_var=/users/enovoa/boguzhan/Software/epinano_basefreq_final/Epinano_base_freq_2.py
python3.7 $tsv_to_var -b Mouse_18s_Processed.bam -s $sam2tsv -R $ref -n 10 
python3.7 $tsv_to_var -b Mouse_18s_Precursor.bam -s $sam2tsv -R $ref -n 10 
python3.7 $tsv_to_var -b Mouse_28s_Processed.sorted.bam -s $sam2tsv -R $ref -n 10 
python3.7 $tsv_to_var -b Mouse_28s_Precursor.sorted.bam -s $sam2tsv -R $ref -n 10 




java -jar /users/enovoa/boguzhan/Software/picard/build/libs/picard.jar CreateSequenceDictionary REFERENCE=yeast_pre_RNA.fasta OUTPUT= yeast_pre_RNA.fasta.dict


ref=yeast_pre_RNA.fasta
sam2tsv=/users/enovoa/boguzhan/Software/jvarkit/dist/sam2tsv.jar
tsv_to_var=/users/enovoa/boguzhan/Software/epinano_basefreq_final/Epinano_base_freq_2.py
python3.7 $tsv_to_var -b Yeast_18s_Processed.bam -s $sam2tsv -R $ref -n 10 
python3.7 $tsv_to_var -b Yeast_18s_Precursor.bam -s $sam2tsv -R $ref -n 10 
```



```bash
#Cleaning the yeast pre-rRNA from the short reads
#Extract reads:
samtools view Yeast_18s_Precursor.sorted.bam | awk '{if (length($10)>300) {print $0}}' > reads.sam
#Extract header:
samtools view -H Yeast_18s_Precursor.sorted.bam > header.txt
#Concat:
cat header.txt reads.sam > sample.sam
#Transform into bam:
samtools view -Sb sample.sam > sample.bam
#Sort and index:
samtools sort sample.bam sample.sorted
samtools index sample.sorted.bam






```


# Epinano Analysis
```R

library(ggplot2)
library(ggrepel)

yeast_18s.precursor <- read.delim("Yeast_18s_Precursor.per.site.baseFreq.csv", sep=",")
yeast_18s.processed <- read.delim("Yeast_18s_Processed.per.site.baseFreq.csv", sep=",")



mouse_18s.precursor <- read.delim("Mouse_18s_Precursor.per.site.baseFreq.csv", sep=",")
mouse_18s.processed <- read.delim("Mouse_18s_Processed.per.site.baseFreq.csv", sep=",")


mouse_28s.precursor <- read.delim("Mouse_28s_Precursor.sorted.per.site.baseFreq.csv", sep=",")
mouse_28s.processed <- read.delim("Mouse_28s_Processed.sorted.per.site.baseFreq.csv", sep=",")





Mouse
5255 : m1ACP3Y
12304 : m3U
9258 : m1A




process_18s <- function(data.precursor, data.processed, label, min.pos , max.pos) {
       data.precursor <- subset(data.precursor, cov > 50)
       data.processed <- subset(data.processed, cov > 50)
       data.precursor <- data.precursor[,c("X.Ref", "pos", "mis")]
       data.processed <- data.processed[,c("X.Ref", "pos", "mis")]
       merged <- merge(data.precursor, data.processed, by.x=c("X.Ref", "pos"), by.y=c("X.Ref", "pos"))
       merged$diff.mis <- abs(merged$mis.x -merged$mis.y)
       merged$Sample <- label
       merged2 <- subset(merged, pos>min.pos  & pos<max.pos)
       return(merged2)
}

yeast_18s <- process_18s(yeast_18s.precursor,yeast_18s.processed,"Yeast_18s", 650, 2505 )
mouse_18s <- process_18s(mouse_18s.precursor,mouse_18s.processed,"Mouse_18s", 4020, 5876 )





process_28s <- function(data.precursor, data.processed, label, min.pos , max.pos) {
       data.precursor <- subset(data.precursor, cov > 30)
       data.processed <- subset(data.processed, cov > 30)
       data.precursor <- data.precursor[,c("X.Ref", "pos", "mis")]
       data.processed <- data.processed[,c("X.Ref", "pos", "mis")]
       merged <- merge(data.precursor, data.processed, by.x=c("X.Ref", "pos"), by.y=c("X.Ref", "pos"))
       merged$diff.mis <- abs(merged$mis.x -merged$mis.y)
       merged$Sample <- label
       merged2 <- subset(merged, pos>min.pos  & pos <max.pos)
       return(merged2)
}


mouse_28s <- process_28s(mouse_28s.precursor,mouse_28s.processed,"Mouse_28s", 8135, 12834 )






scatter_plot <- function(data,label){
       pdf(file=paste(label, "Mismatch_Freq_Precursor_vs_Processed.pdf", sep="_"),height=5,width=5,onefile=FALSE)
       print(ggplot(data, aes(x=mis.x, y=mis.y)) +
                     geom_point(size=1)+
                     geom_abline(slope=1, intercept=0,linetype="dashed", size=0.2, color= "black")+
                     geom_text_repel(data=subset(data, diff.mis> 0.25), aes(label=pos), colour="black",segment.size  = 0.4,segment.color = "grey50",size=5)+
                     xlim(0,1)+
                     ylim(0,1)+
                     xlab("Precursor Mismatch Frequency")+
                     ylab("Processed Mismatch Frequency") +
                     theme_bw()+
                     theme(axis.text.x = element_text(face="bold", color="black",size=11),
                             axis.text.y = element_text(face="bold", color="black", size=11),
                     plot.title = element_text(color="black", size=24, face="bold.italic",hjust = 0.5),
                     axis.title.x = element_text(color="black", size=15, face="bold"),
                     axis.title.y = element_text(color="black", size=15, face="bold"),
                     panel.background = element_blank(),
                     axis.line = element_line(colour = "black", size=0.5),
                     legend.title = element_text(color = "black", size = 20,face="bold"),
                     legend.text = element_text(color = "black", size=20)))
dev.off()
}

scatter_plot(yeast_18s, "Yeast_18s")
scatter_plot(mouse_18s, "Mouse_18s")
scatter_plot(mouse_28s, "Mouse_28s")









