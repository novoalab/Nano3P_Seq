
```R
library(stringr)
library(ggplot2)
library(ggbeeswarm)
library(ggrepel)
library(scales)
library(dplyr)
library(ggridges)

#Import tail
zebrafish.tails <- read.delim("cDNA786327_tails.csv", sep=",")
yeast.tails <- read.delim("cDNA231532_tails.csv", sep=",")
mouse.tails <- read.delim("cDNA964321_tails.csv", sep=",")


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

zebrafish.tails_processed <- manipulate_tail_(zebrafish.tails)
yeast.tails_processed <- manipulate_tail_(yeast.tails)
mouse.tails_processed <- manipulate_tail_(mouse.tails)





# Import the data
# DATA
mouse.data <- read.delim("cDNA964321_all.genome38_sequin_rRNA_ALLRNAs_Merged.bed", header=FALSE)
yeast.data <- read.delim("cDNA231532_SacCer_Genome_rRNA_FINAL.bed", header=FALSE)
zebrafish.data <- read.delim("cDNA786327_6hpf.genome11_sequin_ALLRNAs_Merged.bed", header=FALSE)



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
	#Create a category
	gene_type_sum <- aggregate(.~Gene_Type, merged2[,c("Gene_Type", "Gene_Count_Norm")], sum)
	gene_type_major <- subset(gene_type_sum, Gene_Count_Norm > 250)
	gene_type_major$Category <- gene_type_major$Gene_Type
	gene_type_minor <- subset(gene_type_sum, Gene_Count_Norm < 250)
	gene_type_minor$Category <- "Other"
	gene_type <- rbind(gene_type_major, gene_type_minor)
	colnames(gene_type) <- c("Gene_Type", "Gene_Type_Count_Norm", "Category")
	#MErge
	merged3 <- merge(merged2, gene_type, by.x="Gene_Type", by.y="Gene_Type")
	return(merged3)
}


mouse.reshape <- reshape(mouse.data,mouse.tails_processed,"Mouse")
mouse.mrna <- subset(mouse.reshape, Gene_Type=="protein_coding")
yeast.reshape <- reshape(yeast.data,yeast.tails_processed,"Yeast")
yeast.mrna <- subset(yeast.reshape, Gene_Type=="mRNA")
zebrafish.reshape <- reshape(zebrafish.data,zebrafish.tails_processed,"Zebrafish")
zebrafish.mrna <- subset(zebrafish.reshape, Gene_Type=="protein_coding")

all_species <- rbind(yeast.mrna,mouse.mrna,zebrafish.mrna)




##Overall Tail comparison (Single transcript)
tail_comparison_overall <- function(data, label) {
	median_lengths <- aggregate(.~Sample, data[,c("Sample", "tail_length")], median)
	pdf(file=paste(label, "Overall_Tail_Comparison_Single_Transcript.pdf",sep="_"),height=4,width=6,onefile=FALSE)
	print(ggplot(data, aes(x=tail_length, color=Sample)) +
		geom_density()+
		geom_vline(data=median_lengths, aes(xintercept=tail_length, color=Sample),
             linetype="dashed")+
		geom_text(data=median_lengths, aes(y=0.05, label=tail_length))+
  		theme_bw()+
		xlab("PolyA Tail Length")+
        ylab("Density"))
	dev.off()
}


tail_comparison_overall(all_species, "Cross_Species")






## Tail comparison median per gene
tail_comparison_median_per_gene_protein<- function(data, label) {
	data <-  data[!duplicated(data[c("Gene_Name", "Sample")]),]
	pdf(file=paste(label, "Median_Tail_Per_Gene_Comparison_mRNA.pdf",sep="_"),height=6,width=10,onefile=FALSE)
	print(ggplot(data, aes(x=Median_Length, color=Sample)) +
		geom_density()+
  		theme_bw()+
  		xlab("Median Tail Length")+
        ylab("Density"))
  		#oord_cartesian(xlim=c(-10, 200))
	dev.off()
}
tail_comparison_median_per_gene_protein(all_species, "Cross_Species")









