###################################################
###### ANALYSIS OF THE NUCLEAR TAILFINDR #########
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
tails <- read.delim("cDNA964321_tails.csv", sep=",")


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
data <- read.delim("cDNA964321_all.genome38_sequin_rRNA_ALLRNAs_Merged.bed", header=FALSE)
# Keep important columns
data2 <- data[,c("V1", "V2", "V3", "V4", "V5", "V6", "V16", "V17")]
colnames(data2) <- c("Chr", "Start", "End", "Read_ID", "Quality", "Strand", "Gene_Name", "Gene_Type")

gene_type_count <- function(data) {
	# Make a Count Barplot in order to make a count
	gene_type_counts <- as.data.frame(table(data$Gene_Type))
	colnames(gene_type_counts) <- c("Gene_Type", "Count")
	gene_type_counts_withoutsequin <- subset(gene_type_counts, Gene_Type != "synthetic")
	gene_type_counts_withoutsequin <- subset(gene_type_counts_withoutsequin, Count >0)
	gene_type_counts_withoutsequin <- gene_type_counts_withoutsequin[order(-gene_type_counts_withoutsequin$Count),]
	gene_type_counts_withoutsequin$Gene_Type <- factor(gene_type_counts_withoutsequin$Gene_Type, levels = gene_type_counts_withoutsequin$Gene_Type)
	return(gene_type_counts_withoutsequin)
}

data_gene_type_count <- gene_type_count(data2)
data_gene_type_count_min50 <- subset(data_gene_type_count, Count >50)

pdf(file="Gene_Type_Counts_BothRepMerged_Min50.pdf",height=10,width=20,onefile=FALSE)
print(ggplot(data_gene_type_count_min50, aes(x=Gene_Type, y=Count, fill=Gene_Type)) +
	theme_bw()+
  geom_bar(stat="identity"))
dev.off()






## PIE CHART
gene_type_pie_chart <- function(data) {
	gene_type_counts <- as.data.frame(table(data$Gene_Type))
	colnames(gene_type_counts) <- c("Gene_Type", "Count")
	gene_type_counts_withoutsequin <- subset(gene_type_counts, Gene_Type != "synthetic")
	gene_type_major <- subset(gene_type_counts_withoutsequin, Count > 500)
	gene_type_major$Category_New <- gene_type_major$Gene_Type
	gene_type_major$Gene_Type <- NULL
	gene_type_minor <- subset(gene_type_counts_withoutsequin, Count < 500)
	gene_type_minor2 <- as.data.frame(sum(gene_type_minor$Count))
	gene_type_minor2$Category_New <- "Other"
	colnames(gene_type_minor2) <- c("Count", "Category_New")
	Gene_Type_Fixed <- rbind(gene_type_major, gene_type_minor2)
	Gene_Type_Fixed$Category_New <- factor(Gene_Type_Fixed$Category_New, levels = unique(Gene_Type_Fixed$Category_New[order(-Gene_Type_Fixed$Count )]))
	Gene_Type_Fixed <- Gene_Type_Fixed %>% 
	  arrange(desc(Category_New)) %>%
	  mutate(prop = Count / sum(Gene_Type_Fixed$Count) *100) %>%
	  mutate(ypos = cumsum(prop)- 0.5*prop )
	  return(Gene_Type_Fixed)
	}


data_gene_type_piechart <- gene_type_pie_chart(data2)



# Basic piechart
pdf(file="Gene_Type_Counts_BothRepMerged_PieChart_Min500.pdf",height=10,width=10,onefile=FALSE)
print(ggplot(data_gene_type_piechart, aes(x="", y=prop, fill=Category_New)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() + 
  #theme(legend.position="none") +
 #geom_text_repel(aes(y = ypos, label = Category_New), color = "black", size=6)+
    geom_text_repel(aes(y = ypos, label = paste(Category_New,percent(Count/sum(data_gene_type_piechart$Count)))), color = "black", size=5))
dev.off()


























# First lets cleanup the data (Remove any gene that has less than 20 reads)
cleanup<- function(data) {
	data2 <- subset(data, Quality > 0)
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
data_cleanup <- data2




# Make a Count Barplot in order to make a count
gene_type_counts <- as.data.frame(table(data_cleanup$Gene_Type))
colnames(gene_type_counts) <- c("Gene_Type", "Count")
gene_type_counts_withoutsequin <- subset(gene_type_counts, Gene_Type != "synthetic")
gene_type_counts_withoutsequin <- subset(gene_type_counts_withoutsequin, Count >0)
gene_type_counts_withoutsequin <- gene_type_counts_withoutsequin[order(-gene_type_counts_withoutsequin$Count),]
gene_type_counts_withoutsequin$Gene_Type <- factor(gene_type_counts_withoutsequin$Gene_Type, levels = gene_type_counts_withoutsequin$Gene_Type)


pdf(file="Gene_Type_Counts_BothRepMerged_WithoutFiltering.pdf",height=10,width=20,onefile=FALSE)
print(ggplot(gene_type_counts_withoutsequin, aes(x=Gene_Type, y=Count, color=Gene_Type)) +
  geom_bar(stat="identity", fill="white"))
dev.off()





## PIE CHART
gene_type_major <- subset(gene_type_counts_withoutsequin, Count > 500)
gene_type_major$Category_New <- gene_type_major$Gene_Type
gene_type_major$Gene_Type <- NULL
gene_type_minor <- subset(gene_type_counts_withoutsequin, Count < 500)
gene_type_minor2 <- as.data.frame(sum(gene_type_minor$Count))
gene_type_minor2$Category_New <- "Other"
colnames(gene_type_minor2) <- c("Count", "Category_New")


Gene_Type_Fixed <- rbind(gene_type_major, gene_type_minor2)


Gene_Type_Fixed$Category_New <- factor(Gene_Type_Fixed$Category_New, levels = unique(Gene_Type_Fixed$Category_New[order(-Gene_Type_Fixed$Count )]))

Gene_Type_Fixed <- Gene_Type_Fixed %>% 
  arrange(desc(Category_New)) %>%
  mutate(prop = Count / sum(Gene_Type_Fixed$Count) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )
# Basic piechart
pdf(file="Gene_Type_Counts_BothRepMerged_PieChart_WithoutFiltering.pdf",height=10,width=10,onefile=FALSE)
print(ggplot(Gene_Type_Fixed, aes(x="", y=prop, fill=Category_New)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() + 
  #theme(legend.position="none") +
 #geom_text_repel(aes(y = ypos, label = Category_New), color = "black", size=6)+
    geom_text_repel(aes(y = ypos, label = paste(Category_New,percent(Count/sum(Gene_Type_Fixed$Count)))), color = "black", size=5))
dev.off()













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
data_tails_gene_type_processed_hq <- subset(data_tails_gene_type_processed, Quality > 10)



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
	dotplot(data_tails_gene_type_processed_hq,"Rep_Merged_Processed")







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

scatter_plot(datab_polyApopcount, "Both_Merged")



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
scatter_plot_facet(datab_polyApopcount, "Both_Merged")



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

scatter_plot_protein_coding(protein_coding, "Both_Merged")






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

export_reads(data_per_gene_tails,"Both_Merged")
```