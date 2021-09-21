```bash
#Trying out flAiR

flair=/users/enovoa/boguzhan/Software/flair/flair.py

ref=/users/enovoa/boguzhan/references/danio_rerio/zebrafish_11_sequins.fa
gtf=/users/enovoa/boguzhan/references/danio_rerio/Danio_rerio.GRCz11.103.2.gtf


python $flair align -g $ref -r 2hpf_nonrRNA_both.fastq -o 2hpf_nonrRNA_both_aligned
python $flair align -g $ref -r 4hpf_nonrRNA_both.fastq -o 4hpf_nonrRNA_both_aligned
python $flair align -g $ref -r 6hpf_nonrRNA_both.fastq -o 6hpf_nonrRNA_both_aligned
python $flair align -g $ref -r 246hpf_nonrRNA_both.fastq -o 246hpf_nonrRNA_both_aligned



python $flair correct -q 2hpf_nonrRNA_both_aligned.bed -g $ref -f $gtf
python $flair correct -q 4hpf_nonrRNA_both_aligned.bed -g $ref -f $gtf
python $flair correct -q 6hpf_nonrRNA_both_aligned.bed -g $ref -f $gtf
python $flair correct -q 246hpf_nonrRNA_both_aligned.bed  -g $ref -f $gtf



python $flair collapse  -g $ref -f  $gtf -r 2hpf_nonrRNA_both.fastq -q 2hpf_nonrRNA_both_aligned.bed --generate_map -o 2hpf_nonrRNA_both_aligned
python $flair collapse  -g $ref -f  $gtf -r 4hpf_nonrRNA_both.fastq -q 4hpf_nonrRNA_both_aligned.bed --generate_map -o 4hpf_nonrRNA_both_aligned
python $flair collapse -g $ref -f  $gtf -r 6hpf_nonrRNA_both.fastq -q 6hpf_nonrRNA_both_aligned.bed --generate_map -o 6hpf_nonrRNA_both_aligned
python $flair collapse -g $ref -f  $gtf -r 246hpf_nonrRNA_both.fastq -q 246hpf_nonrRNA_both_aligned.bed --generate_map -o 246hpf_nonrRNA_both_aligned






flair=/users/enovoa/boguzhan/Software/flair/flair.py

ref=/users/enovoa/boguzhan/references/danio_rerio/zebrafish_11_sequins.fa
gtf=/users/enovoa/boguzhan/references/danio_rerio/Danio_rerio.GRCz11.103.2.gtf
python $flair align -g $ref -r 246hpf_nonrRNA_both.fastq -o 246hpf_nonrRNA_both_aligned


qsub -cwd -q long-sl7 -l virtual_free=40G,h_rt=12:00:00 a246hpf.sh
qsub -cwd -q long-sl7 -l virtual_free=40G,h_rt=12:00:00 a2hpf.sh
qsub -cwd -q long-sl7 -l virtual_free=40G,h_rt=12:00:00 a6hpf.sh
qsub -cwd -q long-sl7 -l virtual_free=40G,h_rt=12:00:00 a4hpf.sh




flair=/users/enovoa/boguzhan/Software/flair/flair.py

ref=/users/enovoa/boguzhan/references/danio_rerio/zebrafish_11_sequins.fa
gtf=/users/enovoa/boguzhan/references/danio_rerio/Danio_rerio.GRCz11.103.2.gtf

python $flair collapse -t 8 -s 20 -n best_only --trust_ends --filter nosubset -g $ref -f  $gtf -r 246hpf_nonrRNA_both.fastq -q 246hpf_nonrRNA_both_aligned.bed --generate_map -o 246hpf_nonrRNA_both_aligned_strict3


qsub -cwd -q mem_512 -l virtual_free=200G,h_rt=48:00:00 test2.sh
```

```R
data <- read.csv("246hpf_nonrRNA_both_aligned_strict.isoform.read.map.txt",sep="\t",header=FALSE)

final <- vector()
for (uniq in unique(data$V1)) {
	subs <- subset(data, V1==uniq)
	reads<- data.frame(do.call('rbind', strsplit(as.character(subs$V2),',',fixed=TRUE)))
	reads2<-as.data.frame(t(reads))
	if (nrow(reads2) > 10 ){
	colnames(reads2)<- "read_id"
	reads2$uniq <- uniq
	reads2$transcript <- gsub("_.*","",reads2$uniq)
	reads2$gene <- gsub(".*_","",reads2$uniq)
	final<- rbind(final, reads2) } else{}
}
write.table(final,file="246hpf_final_reads_strict.tsv",quote=FALSE,sep="\t")

```


```python
import pandas as pd
import numpy as np

mydata = pd.read_table("246hpf_nonrRNA_both_aligned_strict.isoform.read.map.txt", sep="\t", names=['uniq', 'reads'])

df = []
for unq in mydata.uniq.unique():
	subs=mydata[mydata["uniq"] == unq]
	subs_reads = subs['reads'].str.split(',',expand=True)
	subs_reads_t = np.transpose(subs_reads)
	subs_reads_t = subs_reads_t.assign(uniq=unq)
	uniq_div = subs_reads_t['uniq'].str.split('_',expand=True)
	subs_final = pd.concat([subs_reads_t.reset_index(drop=True), uniq_div], axis=1)
	del subs_final['uniq']
	#subs_final.columns = subs_final.iloc[0]
	subs_final = subs_final.drop(subs_final.index[[0]])
	df.append(subs_final)
	df = pd.DataFrame()
print(df)

    data.columns = data.iloc[0]
    data = data.drop(data.index[[0]])


final = final.append(pd.DataFrame(data = subs_final), ignore_index=True)

```


```R
library(stringr)
library(ggplot2)
library(ggbeeswarm)
library(ggrepel)
library(scales)
library(dplyr)
library(ggridges)
library(reshape2)
library(EnvStats)

#reads
reads <- read.delim("246hpf_final_reads.tsv")


#Import tail
ribodepleted_rep1.tails <- read.delim("cDNA786327_tails.csv", sep=",")
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

ribodepleted_rep1.tails_processed <- manipulate_tail_(ribodepleted_rep1.tails)
ribodepleted_rep2.tails_processed <- manipulate_tail_(ribodepleted_rep2.tails)


ribodepleted_merged.tails_processed <- rbind(ribodepleted_rep1.tails_processed,ribodepleted_rep2.tails_processed)

ribodep_hpf2_rep1.data <- read.delim("2hpf.genome11_sequin_ALLRNAs_Merged_Rep1.bed", header=FALSE)
ribodep_hpf2_rep2.data <- read.delim("2hpf.genome11_sequin_ALLRNAs_Merged_Rep2.bed", header=FALSE)

ribodep_hpf4_rep1.data <- read.delim("4hpf.genome11_sequin_ALLRNAs_Merged_Rep1.bed", header=FALSE)
ribodep_hpf4_rep2.data <- read.delim("4hpf.genome11_sequin_ALLRNAs_Merged_Rep2.bed", header=FALSE)

ribodep_hpf6_rep1.data <- read.delim("6hpf.genome11_sequin_ALLRNAs_Merged_Rep1.bed", header=FALSE)
ribodep_hpf6_rep2.data <- read.delim("6hpf.genome11_sequin_ALLRNAs_Merged_Rep2.bed", header=FALSE)

ribodep_hpf2_merged.data<- rbind(ribodep_hpf2_rep1.data ,ribodep_hpf2_rep2.data )
ribodep_hpf4_merged.data<- rbind(ribodep_hpf4_rep1.data ,ribodep_hpf4_rep2.data )
ribodep_hpf6_merged.data<- rbind(ribodep_hpf6_rep1.data ,ribodep_hpf6_rep2.data )


# Reshape the tables and remove low quality reads
reshape<- function(data,tails,label) {
  data2 <- data[,c("V1", "V4", "V5", "V6", "V16", "V17")]
  colnames(data2) <- c("Chr","Read_ID", "Quality", "Strand", "Gene_Name", "Gene_Type")
  data3 <- subset(data2, Quality > 30)
  #Merge the data with tails
  merged <- merge(data3, tails, by.x="Read_ID",by.y=c("read_id") )
  merged <- merged[,c("Read_ID", "Chr", "Gene_Name", "Gene_Type", "tail_length")]
  return(merged)
}

ribodep_hpf2_merged.reshape <- reshape(ribodep_hpf2_merged.data,ribodepleted_merged.tails_processed,"Ribodep_2hpf_merged")
ribodep_hpf4_merged.reshape <- reshape(ribodep_hpf4_merged.data,ribodepleted_merged.tails_processed,"Ribodep_4hpf_merged")
ribodep_hpf6_merged.reshape <- reshape(ribodep_hpf6_merged.data,ribodepleted_merged.tails_processed,"Ribodep_6hpf_merged")



merge_reads <- function(data,reads,label) {
	reads_bed_tail <- merge(data, reads, by.x="Read_ID", by.y="read_id")
	reads_bed_tail2<- subset(reads_bed_tail, Gene_Type=="protein_coding")
	reads_bed_tail3 <- reads_bed_tail2[,c("Read_ID","transcript","gene","tail_length")]
    reads_bed_tail3$Gene_Initial <- substr(reads_bed_tail3$gene, 1, 3)
    reads_bed_tail4 <- subset(reads_bed_tail3,Gene_Initial=="ENS" )
    #remove transcripts less than x coverage
    df <- vector()
    for (trans in unique(reads_bed_tail4$transcript)){
    	subs <- subset(reads_bed_tail4, transcript==trans)
    	if (nrow(subs) > 10) {
    		df <- rbind(df, subs) } else {}
    	}
    #remove genes with only one isoform
    df2 <- vector()
   for (genes in unique(df$gene)){
    	subs <- subset(df, gene==genes)
    	if (length(unique(subs$transcript)) > 1) {
    		df2 <- rbind(df2,subs) } else {}
    	}
    #Count isoform
    df3 <- vector()
    for (genes in unique(df2$gene)){
    	subs <- subset(df2, gene==genes)
    	iso_count <- length(unique(subs$transcript)) 
    	subs$isoform_count <- iso_count
    	df3 <- rbind(df3, subs)
    }
    #median length per isoform
    df4 <- vector()
    for (trans in unique(df3$transcript)){
    	subs <- subset(df3, transcript==trans)
    	median_length <- median(subs$tail_length)
    	subs$median_length_isoform <- median_length 
    	df4 <- rbind(df4, subs)
    }
    write.table(df4, file=paste(label,"processed_tail_isoform.tsv", sep="_"), quote=FALSE, sep="\t")
    return(df4)
}

ribodep_hpf2_merged_READS <- merge_reads(ribodep_hpf2_merged.reshape,reads, "2HPF" )
ribodep_hpf4_merged_READS <- merge_reads(ribodep_hpf4_merged.reshape,reads, "4HPF" )
ribodep_hpf6_merged_READS <- merge_reads(ribodep_hpf6_merged.reshape,reads, "6HPF" )


ribodep_hpf2_merged_READS <- read.delim("2HPF_processed_tail_isoform.tsv")
ribodep_hpf4_merged_READS <- read.delim("4HPF_processed_tail_isoform.tsv")
ribodep_hpf6_merged_READS <- read.delim("6HPF_processed_tail_isoform.tsv")







#Kruskal-Wallis test
  compute_kruskal <- function(data) {
    df <- vector()
    for (genes in unique(data$gene)){
      subs <- subset(data, gene==genes)
      kruskal <- kruskal.test(tail_length ~ transcript, data = subs)
      p_value <- kruskal$p.value
      isoform_count <- unique(subs$isoform_count)
      table <- as.data.frame(cbind(genes,p_value,isoform_count ))
      df <- rbind(df, table)
    }
    df$p_value <- as.numeric(as.character(df$p_value))
    df$isoform_count <- as.numeric(as.character(df$isoform_count))
    df <- df[order(df$p_value),]
    return(df)
}

ribodep_hpf2_merged_READS_kruskal <- compute_kruskal(ribodep_hpf2_merged_READS)

library(ggrepel)

dotplot_kruskal <- function(data,label) {
pdf(file=paste(label, "kruskal_isoform_dotplot.pdf",sep="_"),height=5,width=12,onefile=FALSE)
      print(ggplot(data, aes(x=as.factor(isoform_count), y=log(1/p_value))) + 
        geom_quasirandom(varwidth = TRUE, aes())+
        geom_boxplot(aes(alpha=0), outlier.shape=NA)+
        stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                geom = "crossbar", width = 0.7, color="#c06c84")+
        geom_text_repel(data=subset(data, log(1/p_value) >10 ), aes(y = log(1/p_value), label = genes), color = "black", size=6)+
        stat_n_text(size=5)+
        theme_bw()+
        ggtitle(label)+
        xlab("Isoform count")+
        ylab("log(1/Kruskal-Wallis p-value) ") +
        theme(axis.text=element_text(size=14),strip.text = element_text(size=13),
                axis.title=element_text(size=17,face="bold"),
                legend.title = element_text(size = 20),
                legend.text = element_text(color = "black", size=15)))
    dev.off()
}



dotplot_kruskal(ribodep_hpf2_merged_READS_kruskal,"2hpf")




#EXTRACT READS INDIVIDUALLY `PER TRANSCRIPT
extract_reads <- function(data, genes, label) {
  subs <- subset(data, gene==genes)
  for (trans in unique(subs$transcript)) {
    subs2 <- subset(subs, transcript == trans)
    subs3 <- subs2[,"Read_ID"]
    write.table(subs3, file=paste(label, genes, trans, "reads.tsv", sep="_"), col.names=FALSE, row.names=FALSE, quote=FALSE)
  }
}


extract_reads(ribodep_hpf2_merged_READS,"ENSDARG00000056517", "2HPF")












for (Gene in unique(df4$gene)) {
		subs <- subset( df4, gene==Gene)
pdf(file=paste(label,Gene, "tail_length_per_isoform_dotplot.pdf",sep="_"),height=5,width=12,onefile=FALSE)
			print(ggplot(subs, aes(x=transcript, y=tail_length)) + 
				geom_quasirandom(varwidth = TRUE, aes())+
				geom_boxplot(aes(alpha=0), outlier.shape=NA)+
				stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
            		geom = "crossbar", width = 0.7, color="#c06c84")+
				stat_n_text(size=5)+
				theme_bw()+
				ggtitle(label)+
				xlab("Isoform")+
              	ylab("Tail length") +
				theme(axis.text=element_text(size=14),strip.text = element_text(size=13),
            		axis.title=element_text(size=17,face="bold"),
            		legend.title = element_text(size = 20),
            		legend.text = element_text(color = "black", size=15)))
		dev.off()
}

for (Gene in unique(reads_bed_tail3$gene)) {
		subs <- subset( reads_bed_tail3, gene==Gene)
		if (length(unique(subs$transcript)) >1) {
			pdf(file=paste(Gene, "tail_length_per_isoform_boxplot.pdf",sep="_"),height=6,width=6,onefile=FALSE)
			print(ggplot(subs, aes(x=transcript, y=tail_length, fill=transcript)) + 
  			 geom_boxplot() )
			dev.off() 
			} else{}
		}
```

```BASH
java -jar /users/enovoa/boguzhan/Software/picard/build/libs/picard.jar FilterSamReads \
       I=246hpf_nonrRNA_both_aligned.bam \
       O=2HPF_ENSDARG00000056517_b9fea6f3-0c9c-4a86-819f-57b446948369.bam\
       READ_LIST_FILE=2HPF_ENSDARG00000056517_b9fea6f3-0c9c-4a86-819f-57b446948369_reads.tsv\
       FILTER=includeReadList

java -jar /users/enovoa/boguzhan/Software/picard/build/libs/picard.jar FilterSamReads \
       I=246hpf_nonrRNA_both_aligned.bam \
       O=2HPF_ENSDARG00000056517_b0e9fa84-d9e8-468d-91b7-575091794040.bam\
       READ_LIST_FILE=2HPF_ENSDARG00000056517_b0e9fa84-d9e8-468d-91b7-575091794040_reads.tsv\
       FILTER=includeReadList



samtools index 2HPF_ENSDARG00000056517_b9fea6f3-0c9c-4a86-819f-57b446948369.bam
samtools index 2HPF_ENSDARG00000056517_b0e9fa84-d9e8-468d-91b7-575091794040.bam






