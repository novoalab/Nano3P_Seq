###################################################
###### ANALYSIS OF THE Zebrafish pA TAILFINDR #########
###################################################
########## OGUZHAN BEGIK APRIL 2020 ###############


# Analysis of Tails on R 
```R

library(stringr)
library(ggplot2)
library(ggbeeswarm)
library(ggrepel)
library(scales)
library(dplyr)

#Import tail
tails <- read.delim("cDNA852361_tails.csv", sep=",")


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

tails_processed <- manipulate_tail_(tails)





# Import the data
data <- read.delim("cDNA8523612.genome11_sequin_ALLRNAs_Merged.bed", header=FALSE)
# Keep important columns
data2 <- data[,c("V1", "V2", "V3", "V4", "V5", "V6", "V16", "V17")]
colnames(data2) <- c("Chr", "Start", "End", "Read_ID", "Quality", "Strand", "Gene_Name", "Gene_Type")


# First lets cleanup the data (Remove any gene that has less than 20 reads)
cleanup<- function(data) {
	data2 <- subset(data, Quality > 30)
	final<- vector()
		for (gene in unique(data2$Gene_Name)) {
			subs <- subset(data2, Gene_Name==gene)
			if (nrow(subs) > 19) {
				final <- rbind(final, subs)
			} else {
				subs <- subs
			}
		}
		return(final)
	}
data_cleanup <- cleanup(data2)









merge_tails <- function(data, tails) { 
	merged<- merge(data, tails, by.x=c("Read_ID"),by.y=c("read_id") )
	merged$tail_is_valid <- NULL
	merged$read_type <- NULL
	return(merged)
}



data_tails <- merge_tails(data_cleanup, tails_processed)



#Filter for the genetypes
filter_gene_types <- function(data) {
	final <- vector()
	for (type in unique(data$Gene_Type)) {
		subs <- subset(data, Gene_Type == type)
		if (nrow(subs) < 50) {
		subs$Gene_Type <- "Other"
		} else {
			subs$Gene_Type <- subs$Gene_Type 
		}
		final <- rbind(final, subs)
	}
	return(final)
}

data_tails_gene_type_processed <- filter_gene_types(data_tails)



	dotplot <- function(data, label) {
		pdf(file=paste(label, "Tails_All_Genes.pdf",sep="_"),height=10,width=20,onefile=FALSE)
			print(ggplot(data, aes(x=Gene_Type, y=tail_length)) + 
				geom_quasirandom(varwidth = TRUE, aes(color=Gene_Type))+
				geom_boxplot(aes(alpha=0), outlier.shape=NA)+
				stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
            		geom = "crossbar", width = 0.7, color="#c06c84")+
				theme_bw()+
				ggtitle(label)+
				xlab("Group")+
              	ylab("Tail length") +
				theme(axis.text=element_text(size=14),strip.text = element_text(size=13),
            		axis.title=element_text(size=17,face="bold"),
            		legend.title = element_text(size = 20),
            		legend.text = element_text(color = "black", size=15)))
		dev.off()
	}

	#dotplot(data_tails,"Rep_Merged")
	dotplot(data_tails_gene_type_processed,"Zebrafish_pA")







### Per-gene count / PolYA tail analysis
	per_gene_analysis <- function(data) {
	final <- vector()
	for (gene in unique(data$Gene_Name)) { 
		subs <- subset(data, Gene_Name==gene)
		subs$Group <- subs$tail_length
		subs$Group[which(subs$tail_length == 0)] <- "No_PolyA"
		subs$Group[which(subs$tail_length < 10 & subs$tail_length  > 0)]<- "Small_PolyA"
		subs$Group[which(subs$tail_length > 10 )]<- "Long_PolyA"
		if (nrow(subs) > 20) {
			final <- rbind(final, subs)
			} else {
				final <- final
			}
	}
			return(final)
}


data_per_gene_tails <- per_gene_analysis(data_tails)



count_polyA_population <- function(data) {
final <- vector()
pre <- vector()
for (type in unique(data$Gene_Type)) { 
	subs <- subset(data, Gene_Type==type)
	for (gene in unique(subs$Gene_Name)) {
		subs2 <- subset(subs, Gene_Name==gene)
		if (nrow(subs2) > 20) {
			sub_table <- table(subs2$Group)
			pre$gene <- gene
			pre$type <- type
			pre <- as.data.frame(pre)
			pre$Long_PolyA_Count <- as.numeric(sub_table["Long_PolyA"])
			pre$No_PolyA_Count <- as.numeric(sub_table["No_PolyA"])
			final<- rbind(final,pre)
			} else {
				subs2 <- subs2
			}
		}
	}

	final[["Long_PolyA_Count"]][is.na(final[["Long_PolyA_Count"]])] <- 0
	final[["No_PolyA_Count"]][is.na(final[["No_PolyA_Count"]])] <- 0
	final$Total <- final$Long_PolyA_Count + final$No_PolyA_Count
	final$Ratio_LongPolyA <- final$Long_PolyA_Count /final$Total
	final$Ratio_NoPolyA<- final$No_PolyA_Count /final$Total
	return(final)
}


datab_polyApopcount <- count_polyA_population(data_per_gene_tails)




library(ggrepel)

scatter_plot <- function(data,label){
	pdf(file=paste(label, "PolyA_Ratio_Scatter_Plot.pdf", sep="_"),height=5,width=10,onefile=FALSE)
	print(ggplot(data, aes(x=Ratio_LongPolyA, y=Ratio_NoPolyA, colour=type)) +
			geom_point(size=1)+
			geom_abline(slope=1, intercept=0,linetype="dashed", size=0.2, color= "black")+
			geom_text_repel(data=subset(data, Ratio_NoPolyA>0.25), aes(label=gene), colour="black",segment.size  = 0.4,segment.color = "grey50",size=5)+
			xlim(0,1)+
			ylim(0,1)+
			xlab("Long PolyA Ratio")+
			ylab("No PolyA Ratio") +
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

scatter_plot(datab_polyApopcount, "Zebrafish_pA")



scatter_plot_facet <- function(data,label){
	pdf(file=paste(label, "PolyA_Ratio_Scatter_Plot_Facet.pdf", sep="_"),height=15,width=15,onefile=FALSE)
	print(ggplot(data, aes(x=Ratio_LongPolyA, y=Ratio_NoPolyA, colour=type)) +
			geom_point(size=1)+
			geom_abline(slope=1, intercept=0,linetype="dashed", size=0.2, color= "black")+
			geom_text_repel(data=subset(data, type =! "protein_coding"), aes(label=gene), colour="black",segment.size  = 0.4,segment.color = "grey50",size=5)+
			xlim(0,1)+
			ylim(0,1)+
			xlab("Long PolyA Ratio")+
			ylab("No PolyA Ratio") +
			theme_bw()+
			facet_wrap(~type)+
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
scatter_plot_facet(datab_polyApopcount, "Zebrafish_pA")



protein_coding <- subset(datab_polyApopcount, type=="protein_coding")


scatter_plot_protein_coding <- function(data,label){
	pdf(file=paste(label, "PolyA_Ratio_Scatter_Plot_Protein_Coding.pdf", sep="_"),height=5,width=8,onefile=FALSE)
	print(ggplot(data, aes(x=Ratio_LongPolyA, y=Ratio_NoPolyA, colour=type)) +
			geom_point(size=1)+
			geom_abline(slope=1, intercept=0,linetype="dashed", size=0.2, color= "black")+
			geom_text_repel(data=subset(data, Ratio_NoPolyA>0.1, ), aes(label=gene), colour="black",segment.size  = 0.4,segment.color = "grey50",size=5)+
			xlim(0,1)+
			ylim(0,1)+
			xlab("Long PolyA Ratio")+
			ylab("No PolyA Ratio") +
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

scatter_plot_protein_coding(protein_coding, "Zebrafish_pA")






	dotplot <- function(data, label) {
		pdf(file=paste(label, "tails_dotplot_proteincoding.pdf",sep="_"),height=6,width=5,onefile=FALSE)
			print(ggplot(data, aes(x=type, y=Ratio_LongPolyA)) + 
				geom_quasirandom(varwidth = TRUE, aes())+
				geom_boxplot(aes(alpha=0), outlier.shape=NA)+
				stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
            		geom = "crossbar", width = 0.7, color="#c06c84")+
				geom_text_repel(data=subset(data, Ratio_LongPolyA<0.9 ), aes(label=gene), colour="black",segment.size  = 0.4,segment.color = "grey50",size=5)+
				theme_bw()+
				ggtitle(label)+
				xlab("Group")+
              	ylab("Long Tail Fraction") +
				theme(axis.text=element_text(size=14),strip.text = element_text(size=13),
            		axis.title=element_text(size=17,face="bold"),
            		legend.title = element_text(size = 20),
            		legend.text = element_text(color = "black", size=15)))
		dev.off()
	}

	dotplot(protein_coding,"Protein_coding")

	

export_reads <- function(data,label) {
	for (group in unique(data$Group)) {
		subs <- subset(data, Group==group)
		subs2 <- subs[,c("Read_ID")]
		write.table(subs2, file=paste(group, label,"read_id.tsv",sep="_"), sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
	}
}

export_reads(data_per_gene_tails,"Zebrafish_pA")
```



```bash

### EXTRACT THE BAM FOR SMALL RNA READS
java -jar /users/enovoa/boguzhan/Software/picard/build/libs/picard.jar FilterSamReads \
       I=cDNA8523612.rRNA.sorted.bam \
       O=cDNA8523612.rRNA.NoPolyA.bam\
       READ_LIST_FILE=No_PolyA_Zebrafish_pA_read_id.tsv \
       FILTER=includeReadList

samtools sort cDNA8523612.rRNA.NoPolyA.bam cDNA8523612.rRNA.NoPolyA.sorted
samtools index cDNA8523612.rRNA.NoPolyA.sorted.bam




java -jar /users/enovoa/boguzhan/Software/picard/build/libs/picard.jar FilterSamReads \
       I=cDNA8523612.rRNA.sorted.bam \
       O=cDNA8523612.rRNA.LongPolyA.bam\
       READ_LIST_FILE=Long_PolyA_Zebrafish_pA_read_id.tsv \
       FILTER=includeReadList

samtools sort cDNA8523612.rRNA.LongPolyA.bam cDNA8523612.rRNA.LongPolyA.sorted
samtools index cDNA8523612.rRNA.LongPolyA.sorted.bam



