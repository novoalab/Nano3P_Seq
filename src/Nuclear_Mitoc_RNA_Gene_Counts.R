###################################################
###### ANALYSIS OF THE NUCLEAR RUN GENE COUNTS#####
###################################################
########## OGUZHAN BEGIK APRIL 2020 ###############



library(stringr)
library(ggplot2)
library(ggbeeswarm)
library(ggrepel)
library(scales)
library(dplyr)

# Import the data
data <- read.delim("cDNA964321_all.genome38_sequin_rRNA_ALLRNAs_Merged.bed", header=FALSE)
# Keep important columns
data2 <- data[,c("V1", "V2", "V3", "V4", "V5", "V6", "V16", "V17")]
colnames(data2) <- c("Chr", "Start", "End", "Read_ID", "Quality", "Strand", "Gene_Name", "Gene_Type")




######################################################################
######### PLOT THE COUNTS AS THEY ARE : INCLUDE ALL THE READS #######
######################################################################

gene_type_count <- function(data) {
	# Make a Count Barplot in order to make a count
	data <- subset(data, Quality > 30)
	gene_type_counts <- as.data.frame(table(data$Gene_Type))
	colnames(gene_type_counts) <- c("Gene_Type", "Count")
	gene_type_counts2 <- subset(gene_type_counts, Gene_Type != "synthetic")
	gene_type_counts2 <- subset(gene_type_counts2, Count >0)
	gene_type_counts_major <- subset(gene_type_counts2, Count >200)
	gene_type_counts_major$Category_New  <-  gene_type_counts_major$Gene_Type
	gene_type_counts_minor <- subset(gene_type_counts2, Count <200)
	gene_type_counts_minor$Category_New  <-  "Other"
	gene_type_counts3 <- rbind(gene_type_counts_major,gene_type_counts_minor)
	gene_type_counts3 <- gene_type_counts3[order(-gene_type_counts3$Count),]
	gene_type_counts3$Gene_Type <- factor(gene_type_counts3$Gene_Type, levels = gene_type_counts3$Gene_Type)
	return(gene_type_counts3)
}

data_gene_type_count <- gene_type_count(data2)



pdf(file="Gene_Type_Counts_Barplot.pdf",height=10,width=20,onefile=FALSE)
print(ggplot(data_gene_type_count, aes(x=Gene_Type, y=Count, fill=Gene_Type)) +
	theme_bw()+
  geom_bar(stat="identity"))
dev.off()






## PIE CHART
gene_type_pie_chart <- function(data) {
	gene_type_counts <- as.data.frame(table(data$Gene_Type))
	colnames(gene_type_counts) <- c("Gene_Type", "Count")
	gene_type_counts_withoutsequin <- subset(gene_type_counts, Gene_Type != "synthetic")
	gene_type_major <- subset(gene_type_counts_withoutsequin, Count > 200)
	gene_type_major$Category_New <- gene_type_major$Gene_Type
	gene_type_major$Gene_Type <- NULL
	gene_type_minor <- subset(gene_type_counts_withoutsequin, Count < 200)
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
pdf(file="Gene_Type_Counts_PieChart_Min200.pdf",height=10,width=10,onefile=FALSE)
print(ggplot(data_gene_type_piechart, aes(x="", y=prop, fill=Category_New)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() + 
  #theme(legend.position="none") +
 #geom_text_repel(aes(y = ypos, label = Category_New), color = "black", size=6)+
    geom_text_repel(aes(y = ypos, label = paste(Category_New,percent(Count/sum(data_gene_type_piechart$Count)))), color = "black", size=5))
dev.off()







#################################################################################
######### REMOVE THE READS COMING FROM GENE THAT HAS LESS THAN 20 COVERAGE ######
#################################################################################





# First lets cleanup the data (Remove any gene that has less than 20 reads)
cleanup<- function(data) {
	final<- vector()
		for (gene in unique(data$Gene_Name)) {
			subs <- subset(data, Gene_Name==gene)
			if (nrow(subs) > 19) {
				final <- rbind(final, subs)
			} else {
				subs <- subs
			}
		}
		return(final)
	}
data_cleanup <- cleanup(data2)






gene_type_count <- function(data) {
	# Make a Count Barplot in order to make a count
	data <- subset(data, Quality > 30)
	gene_type_counts <- as.data.frame(table(data$Gene_Type))
	colnames(gene_type_counts) <- c("Gene_Type", "Count")
	gene_type_counts2 <- subset(gene_type_counts, Gene_Type != "synthetic")
	gene_type_counts2 <- subset(gene_type_counts2, Count >0)
	gene_type_counts_major <- subset(gene_type_counts2, Count >200)
	gene_type_counts_major$Category_New  <-  gene_type_counts_major$Gene_Type
	gene_type_counts_minor <- subset(gene_type_counts2, Count <200)
	gene_type_counts_minor$Category_New  <-  "Other"
	gene_type_counts3 <- rbind(gene_type_counts_major,gene_type_counts_minor)
	gene_type_counts3 <- gene_type_counts3[order(-gene_type_counts3$Count),]
	gene_type_counts3$Gene_Type <- factor(gene_type_counts3$Gene_Type, levels = gene_type_counts3$Gene_Type)
	return(gene_type_counts3)
}

data_gene_type_count_filtered <- gene_type_count(data_cleanup)




pdf(file="Gene_Type_Counts_GeneMin20Cov_Barplot.pdf",height=10,width=20,onefile=FALSE)
print(ggplot(data_gene_type_count_filtered, aes(x=Gene_Type, y=Count, fill=Gene_Type)) +
	theme_bw()+
  geom_bar(stat="identity"))
dev.off()