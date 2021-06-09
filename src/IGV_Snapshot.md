### Getting an IGV Snapshot for figure


### Using Rep2
```bash
porechop -i 4hpf_nonrRNA.fastq> 4hpf_nonrRNA_trimmed_porechop.fastq


# Map these reads to Genome + rRNA
ref=/users/enovoa/boguzhan/references/danio_rerio/zebrafish_11_sequins.fa
minimap2 -ax splice -k14 -uf --MD $ref 2hpf_nonrRNA_trimmed_porechop.fastq | samtools view -hSb -F 3844 - >  2hpf_nonrRNA_trimmed.genome11_sequin.bam

minimap2 -ax splice -k14 -uf --MD $ref 4hpf_nonrRNA_trimmed_porechop.fastq | samtools view -hSb -F 3844 - >  4hpf_nonrRNA_trimmed.genome11_sequin.bam


minimap2 -ax splice -k14 -uf --MD $ref 4hpf_nonrRNA_trimmed.fastq | samtools view -hSb -F 3844 - >  4hpf_nonrRNA_trimmed.genome11_sequin.bam
minimap2 -ax splice -k14 -uf --MD $ref 6hpf_nonrRNA_trimmed.fastq | samtools view -hSb -F 3844 - >  6hpf_nonrRNA_trimmed.genome11_sequin.bam


samtools sort 2hpf_nonrRNA_trimmed.genome11_sequin.bam 2hpf_nonrRNA_trimmed.genome11_sequin.sorted && rm 2hpf_nonrRNA_trimmed.genome11_sequin.bam && samtools index 2hpf_nonrRNA_trimmed.genome11_sequin.sorted.bam

samtools sort 4hpf_nonrRNA_trimmed.genome11_sequin.bam 4hpf_nonrRNA_trimmed.genome11_sequin.sorted && rm 4hpf_nonrRNA_trimmed.genome11_sequin.bam && samtools index 4hpf_nonrRNA_trimmed.genome11_sequin.sorted.bam


for i in *.sorted.bam;do samtools view -F 4 $i | cut -f1 | sort | uniq | wc -l;echo $i; done






ref=pfn2.fasta

minimap2 -ax splice -k14 -uf --MD $ref 2hpf_nonrRNA_trimmed_porechop.fastq | samtools view -hSb -F 3844 - >  2hpf_nonrRNA_trimmed.pfn2.bam


samtools sort 2hpf_nonrRNA_trimmed.pfn2.bam 2hpf_nonrRNA_trimmed.pfn2.sorted && rm 2hpf_nonrRNA_trimmed.pfn2.bam && samtools index 2hpf_nonrRNA_trimmed.pfn2.sorted.bam



### FOR MOUSE



porechop -i cDNA964321.fastq > cDNA964321_porechop.fastq




# Map these reads to Genome + rRNA
ref=/users/enovoa/boguzhan/references/danio_rerio/zebrafish_11_sequins.fa


ref=/users/enovoa/boguzhan/references/mouse_with_sequins/Mus_musculus.GRCm38.dna_sm.primary_assembly_chrIS_rRNA.fa


minimap2 -ax splice -k14 -uf --MD $ref cDNA964321_porechop.fastq > cDNA964321_porechop_genome38_sequin_rrna.bam


samtools view -hSb -F 3844 cDNA964321_porechop_genome38_sequin_rrna.bam >  cDNA964321_porechop_genome38_sequin_rrna.sam


samtools sort cDNA964321_porechop_genome38_sequin_rrna.sam cDNA964321_porechop_genome38_sequin_rrna.sorted && samtools index cDNA964321_porechop_genome38_sequin_rrna.sorted.bam





#For Curlcake
porechop -i cDNA345234.fastq > cDNA345234_porechop.fastq

ref=/users/enovoa/boguzhan/references/curlcakes/curlcake_1_2.fasta


minimap2 -ax map-ont --MD $ref cDNA345234_porechop.fastq > cDNA345234_porechop.bam

samtools view -hSb -F 3844 cDNA345234_porechop.bam >  cDNA345234_porechop.sam


samtools view -hb -f 0x10 cDNA345234_porechop.sam | samtools sort - cDNA345234_porechop.sorted && samtools index cDNA345234_porechop.sorted.bam









```


# Intersect the rest of the bam file with the rest of the GTF annotation

```bash
bedtools bamtobed -i 2hpf_nonrRNA_trimmed.genome11_sequin.sorted.bam > 2hpf_nonrRNA_trimmed.genome11_sequin.sorted.bed
```

#Manipulate the BED file
```R

data.input <- read.delim("2hpf_nonrRNA_trimmed.genome11_sequin.sorted.bed", header=FALSE)


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

write.table(data, file="2hpf_nonrRNA_trimmed.genome11_sequin.readstarts.bed", quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)
```




#Intersects BED files  (BAM and GTF-Transcript)
```bash

bedtools intersect -a 2hpf_nonrRNA_trimmed.genome11_sequin.readstarts.bed -b Danio_rerio_11_Sequin4_Transcript_Ends.bed -wa -wb > 2hpf_nonrRNA_trimmed.genome11_sequin.overlapping.bed



#Remove duplicate reads
awk '!seen[$4]++' 2hpf_nonrRNA_trimmed.genome11_sequin.overlapping.bed| cut -f4 > 2hpf_nonrRNA_trimmed.genome11_sequin.overlapping.reads



java -jar /users/enovoa/boguzhan/Software/picard/build/libs/picard.jar FilterSamReads \
       I=2hpf_nonrRNA_trimmed.genome11_sequin.sorted.bam \
       O=2hpf_nonrRNA_trimmed.genome11_sequin.overlap.bam\
       READ_LIST_FILE=2hpf_nonrRNA_trimmed.genome11_sequin.overlapping.reads \
       FILTER=includeReadList

samtools sort 2hpf_nonrRNA_trimmed.genome11_sequin.overlap.bam 2hpf_nonrRNA_trimmed.genome11_sequin.overlap.sorted
samtools index 2hpf_nonrRNA_trimmed.genome11_sequin.overlap.sorted.bam

```














```R
library(stringr)
library(ggplot2)
library(ggbeeswarm)
library(ggrepel)
library(scales)
library(dplyr)
library(ggridges)

#Import tail
ribodepleted_rep2.tails <- read.delim("cDNA123791_tails.csv", sep=",")



manipulate_tail_<- function(data) { 
  data2 <- data
  data2$tail_length <- as.numeric(as.character(data2$tail_length))
  data2$tail_start <- as.numeric(as.character(data2$tail_start))
  data2$tail_end <- as.numeric(as.character(data2$tail_end))
  data2$samples_per_nt <- as.numeric(as.character(data2$samples_per_nt))
  data2$samples_per_nt <- data2$samples_per_nt*1.2
  data2$tail_length <- (data2$tail_end-data2$tail_start)/data2$samples_per_nt
  data2[["tail_length"]][is.na(data2[["tail_length"]])] <- 0
  data_filt <- subset(data2, tail_is_valid=="TRUE")
  data_filt_T <- subset(data_filt, read_type=="polyT")

  return(data_filt_T)
}

ribodepleted_rep2.tails_processed <- manipulate_tail_(ribodepleted_rep2.tails)





# Import the data
#RIBODEPLETED DATA
ribodep_hpf2_rep2.data <- read.delim("2hpf.genome11_sequin_ALLRNAs_Merged_Rep2.bed", header=FALSE)
ribodep_hpf4_rep2.data <- read.delim("4hpf.genome11_sequin_ALLRNAs_Merged_Rep2.bed", header=FALSE)
ribodep_hpf6_rep2.data <- read.delim("6hpf.genome11_sequin_ALLRNAs_Merged_Rep2.bed", header=FALSE)






# Reshape the tables and remove low quality reads
reshape<- function(data,tails,label) {
  data2 <- data[,c("V1", "V4", "V5", "V6", "V16", "V17")]
  colnames(data2) <- c("Chr","Read_ID", "Quality", "Strand", "Gene_Name", "Gene_Type")
  data3 <- subset(data2, Quality > 30)
  #Merge the data with tails
  merged <- merge(data3, tails, by.x="Read_ID",by.y=c("read_id") )
  merged <- merged[,c("Read_ID", "Chr", "Gene_Name", "Gene_Type", "tail_length")]
  ourdata <- merged[,c("Gene_Name", "tail_length")]
  ourdata2 <- do.call(data.frame,(aggregate(. ~Gene_Name, data =ourdata, FUN = function(ourdata) c(median=median(ourdata), mean = mean(ourdata), count = length(ourdata) ) )))
  colnames(ourdata2) <- c("Gene_Name", "Median_Length", "Mean_Length", "Gene_Count")
  coverage <- sum(ourdata2$Gene_Count)
  #Merge the data with the stats tabls
  merged2 <- merge(merged, ourdata2, by.x=c("Gene_Name"), by.y=c("Gene_Name"))
  merged2$Gene_Count_Norm <- merged2$Gene_Count/coverage *10000
  merged2$Sample <- label
  merged3 <-  merged2[!duplicated(merged2[c("Gene_Name")]),]
  #Create a category
  gene_type_sum <- aggregate(.~Gene_Type, merged3[,c("Gene_Type", "Gene_Count_Norm")], sum)
  gene_type_major <- subset(gene_type_sum, Gene_Count_Norm > 5.5)
  gene_type_major$Category <- gene_type_major$Gene_Type
  gene_type_minor <- subset(gene_type_sum, Gene_Count_Norm < 5.5)
  gene_type_minor$Category <- "Other"
  gene_type <- rbind(gene_type_major, gene_type_minor)
  colnames(gene_type) <- c("Gene_Type", "Gene_Type_Count_Norm", "Category")
  #MErge
  merged4 <- merge(merged2, gene_type, by.x="Gene_Type", by.y="Gene_Type")
  return(merged4)
}
ribodep_hpf2_rep2.reshape <- reshape(ribodep_hpf2_rep2.data,ribodepleted_rep2.tails_processed,"Ribodep_2hpf_rep2")
ribodep_hpf4_rep2.reshape <- reshape(ribodep_hpf4_rep2.data,ribodepleted_rep2.tails_processed,"Ribodep_4hpf_rep2")
ribodep_hpf6_rep2.reshape <- reshape(ribodep_hpf6_rep2.data,ribodepleted_rep2.tails_processed,"Ribodep_6hpf_rep2")





#### CHECK BINOMIALITY

library(diptest)
bimodality_test <- function(data) {
	bimodal_table <- vector()
	for (gene in unique(data$Gene_Name)) {
		subs <- subset(data, Gene_Name==gene)
		test <- (dip.test(subs$tail_length, simulate.p.value = FALSE, B = 2000))
		statistic <- as.numeric(as.character(test$statistic))
		p_value <- as.numeric(as.character(test$p.value))
		gene_name <- gene
		data_table <- as.data.frame(cbind(gene_name,statistic, p_value))
		data_table$statistic <- as.numeric(as.character(data_table$statistic))
		data_table$p_value <- as.numeric(as.character(data_table$p_value))
		bimodal_table <- rbind(bimodal_table,data_table)

	}
	data2 <- merge(data, bimodal_table, by.x="Gene_Name", by.y="gene_name")
	return(data2)
}




ribodep_hpf2.bimodality <- bimodality_test(ribodep_hpf2_rep2.reshape)

ribodep_hpf2_bimodal <- subset(ribodep_hpf2.bimodality, statistic > 0.05)



plot_boxplot <- function(data, label){
pdf(file=paste(label,"Bimodal_Tails.pdf",sep="_"),height=10,width=50,onefile=FALSE)
	print(ggplot(data, aes(x=Gene_Name, y=tail_length)) + 
		geom_quasirandom(varwidth = TRUE, aes(color=statistic))+
		geom_boxplot(aes(alpha=0), outlier.shape=NA)+
		stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
    		geom = "crossbar", width = 0.7, color="#c06c84")+
		theme_bw()+
		xlab("Gene_Name")+
      	ylab("Tail length") +
		theme(axis.text=element_text(size=14),strip.text = element_text(size=13),
    		axis.title=element_text(size=17,face="bold"),
    		legend.title = element_text(size = 20),
    		legend.text = element_text(color = "black", size=15)))
dev.off()
}

plot_boxplot(ribodep_hpf2_bimodal, "2HPF")






### Per-gene count / PolYA tail analysis


### Per-gene count / PolYA tail analysis
  per_gene_analysis <- function(data) {
  data <- subset(data, Gene_Count > 20)
  data$Group <- data$tail_length
  data$Group[which(data$tail_length == 0)] <- "No_PolyA"
  data$Group[which(data$tail_length < 10 & data$tail_length  > 0)]<- "Small_PolyA"
  data$Group[which(data$tail_length > 10 )]<- "Long_PolyA"
  return(data)
}


ribodep_hpf2_rep2_per_gene_tails <- per_gene_analysis(ribodep_hpf2_rep2.reshape)





export_reads <- function(data,label) {
	for (group in unique(data$Group)) {
		subs <- subset(data, Group==group)
		subs2 <- subs[,c("Read_ID")]
		write.table(subs2, file=paste(group, label,"read_id.tsv",sep="_"), sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
	}
}

export_reads(ribodep_hpf2_rep2_per_gene_tails,"Zebrafish_2HPF")
```


```bash

### EXTRACT THE BAM FOR SMALL RNA READS
java -jar /users/enovoa/boguzhan/Software/picard/build/libs/picard.jar FilterSamReads \
       I=2hpf_nonrRNA_trimmed.genome11_sequin.overlap.sorted.bam \
       O=2hpf_nonrRNA_trimmed.genome11_sequin.NoPolyA.bam \
       READ_LIST_FILE=No_PolyA_Zebrafish_2HPF_read_id.tsv \
       FILTER=includeReadList

samtools sort 2hpf_nonrRNA_trimmed.genome11_sequin.NoPolyA.bam 2hpf_nonrRNA_trimmed.genome11_sequin.NoPolyA.sorted
samtools index 2hpf_nonrRNA_trimmed.genome11_sequin.NoPolyA.sorted.bam


java -jar /users/enovoa/boguzhan/Software/picard/build/libs/picard.jar FilterSamReads \
       I=2hpf_nonrRNA_trimmed.genome11_sequin.overlap.sorted.bam \
       O=2hpf_nonrRNA_trimmed.genome11_sequin.LongPolyA.bam \
       READ_LIST_FILE=Long_PolyA_Zebrafish_2HPF_read_id.tsv \
       FILTER=includeReadList

samtools sort 2hpf_nonrRNA_trimmed.genome11_sequin.LongPolyA.bam 2hpf_nonrRNA_trimmed.genome11_sequin.LongPolyA.sorted
samtools index 2hpf_nonrRNA_trimmed.genome11_sequin.LongPolyA.sorted.bam


