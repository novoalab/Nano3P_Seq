##################################################
######## ANALYSIS OF THE Mouse RUN #########
###################################################
########## OGUZHAN BEGIK APRIL 2020 ###############

```R
library(stringr)
library(ggplot2)
library(ggbeeswarm)
library(ggrepel)
library(scales)
library(dplyr)
library(ggridges)
library(EnvStats)
#IMPORT AND MANIPULATE THE TAIL DATA#############################
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

mouse.tails_processed <- manipulate_tail_(mouse.tails)
####################################################################



#IMPORT AND MANIPULATE THE BED FILE##################################
#Mouse DATA
mouse.data <- read.delim("cDNA964321_all.genome38_sequin_rRNA_ALLRNAs_Merged.bed", header=FALSE)
mouse_rep1.data  <- read.delim("cDNA964321_all.genome38_sequin_rRNA_ALLRNAs_Rep1.bed", header=FALSE) 
mouse_rep2.data  <- read.delim("cDNA964321_all.genome38_sequin_rRNA_ALLRNAs_Rep2.bed", header=FALSE) 


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
	gene_type_major <- subset(gene_type_sum, Gene_Count_Norm > 2)
	gene_type_major$Category <- gene_type_major$Gene_Type
	gene_type_minor <- subset(gene_type_sum, Gene_Count_Norm < 2)
	gene_type_minor$Category <- "Other"
	gene_type <- rbind(gene_type_major, gene_type_minor)
	colnames(gene_type) <- c("Gene_Type", "Gene_Type_Count_Norm", "Category")
	#MErge
	merged4 <- merge(merged2, gene_type, by.x="Gene_Type", by.y="Gene_Type")
	return(merged4)
}
mouse.reshape <- reshape(mouse.data,mouse.tails_processed,"mouse")
mouse_rep1.reshape <- reshape(mouse_rep1.data,mouse.tails_processed,"mouse")
mouse_rep2.reshape <- reshape(mouse_rep2.data,mouse.tails_processed,"mouse")


mouse.unique <-  mouse.reshape[!duplicated(mouse.reshape[c("Gene_Name", "Sample")]),]
write.table(mouse.unique, file="Mouse_RepsMerged_Reshaped.tsv", sep="\t", quote=FALSE,row.names=FALSE)
####################################################################


## REP1 vs REP2


mouse_rep1.unique <-  mouse_rep1.reshape[!duplicated(mouse_rep1.reshape[c("Gene_Name", "Sample")]),]
mouse_rep2.unique <-  mouse_rep2.reshape[!duplicated(mouse_rep2.reshape[c("Gene_Name", "Sample")]),]




columns1 <- c("Gene_Type", "Gene_Name", "Median_Length", "Gene_Count", "Gene_Count_Norm")
columns2 <- c("Gene_Type",  "Gene_Name")


mouse_rep1_vs_rep2 <- merge(mouse_rep1.unique[,columns1] ,mouse_rep2.unique[,columns1], by.x=columns2, by.y=columns2  )


#SCATTER PLOT 
Scatter_Comparison_Gene_Count <- function(data, label, gene_type){
   data$Gene_Count_Norm.x <- log(data$Gene_Count_Norm.x+1)
   data$Gene_Count_Norm.y <- log(data$Gene_Count_Norm.y+1)
   data$Gene_Count_Norm.x <- rescale(data$Gene_Count_Norm.x, to = c(0, 1))
   data$Gene_Count_Norm.y <- rescale(data$Gene_Count_Norm.y, to = c(0, 1))
   subs <- subset(data, Gene_Type==gene_type)
   corr <- cor.test(subs$Gene_Count_Norm.x, subs$Gene_Count_Norm.y, method = "pearson", conf.level = 0.95)
   value <- as.numeric(corr$estimate)
    pdf(file=paste(label,gene_type,"Normalized_Gene_Count_ScatterPlot.pdf", sep="_"),height=5,width=5,onefile=FALSE)
    print(ggplot(subs, aes(x=Gene_Count_Norm.x, y=Gene_Count_Norm.y)) + 
    theme_bw()+ 
    ggtitle(paste("Pearson correlation : ", value))+
    xlab("Rep1 Normalized Gene Count")+
    ylab("Rep2  Normalized Gene Count")+
      geom_point())
      dev.off()
}
Scatter_Comparison_Gene_Count(mouse_rep1_vs_rep2, "Mouse_rep1_vs_rep2", "synthetic")






mouse_rep1_vs_rep2_sequin <- subset(mouse_rep1_vs_rep2, Gene_Type=="synthetic")
mouse_rep1_vs_rep2_sequin2 <- subset(mouse_rep1_vs_rep2_sequin, Gene_Count.x>20 & Gene_Count.y>20 )



mouse_rep1_vs_rep2_sequin2$Gene_Count_Norm.x <- log(mouse_rep1_vs_rep2_sequin2$Gene_Count_Norm.x+1)
   mouse_rep1_vs_rep2_sequin2$Gene_Count_Norm.y <- log(mouse_rep1_vs_rep2_sequin2$Gene_Count_Norm.y+1)
   mouse_rep1_vs_rep2_sequin2$Gene_Count_Norm.x <- rescale(mouse_rep1_vs_rep2_sequin2$Gene_Count_Norm.x, to = c(0, 1))
   mouse_rep1_vs_rep2_sequin2$Gene_Count_Norm.y <- rescale(mouse_rep1_vs_rep2_sequin2$Gene_Count_Norm.y, to = c(0, 1))



plot_denscols_with_corr_pearson<-function(pdfname,my_x,my_y,xlab,ylab) {
  pdf(file=paste(pdfname, "_pearson.pdf",sep=""), height=6, width=6)
  dcols<-densCols(my_x,my_y, colramp=colorRampPalette(blues9[-(1:3)]))
  plot(my_x,my_y,col=dcols,cex=1, cex.lab=1,cex.main=3,lwd=5,pch=20,xlab=xlab,ylab=ylab)
  title(main=pdfname, col.main="black", font.main=4)
  abline(0,1,lwd=3)
  # Correlation
  test<-cor.test(my_x,my_y, method="pearson")
  print(test)
  cor222<-paste("Pearson's r =",round(as.numeric(test$estimate),3))
  #pval<-paste("Pval =",test$p.value)
  mtext(paste(cor222))
  #mtext(paste(cor222,pval,sep=" ; ")) #Print the subtitle with the dataset correlation
  dev.off()
}

plot_denscols_with_corr_pearson("Mouse_Rep1_Rep2_MedianTail_Sequins", mouse_rep1_vs_rep2_sequin2$Median_Length.x, mouse_rep1_vs_rep2_sequin2$Median_Length.y, "Rep1_Median_pA_Length", "Rep2_Median_pA_Length" )


plot_denscols_with_corr_pearson("Mouse_Rep1_Rep2_Gene_Count_Sequins", mouse_rep1_vs_rep2_sequin2$Gene_Count_Norm.x, mouse_rep1_vs_rep2_sequin2$Gene_Count_Norm.y, "Rep1_Median_Gene_Count", "Rep2_Median_Gene_Count" )


plot_denscols_with_corr_pearson("Mouse_Rep1_Rep2_logGene_Count_Sequins", log(mouse_rep1_vs_rep2_sequin2$Gene_Count.x+1), log(mouse_rep1_vs_rep2_sequin2$Gene_Count.y+1), "Rep1_log(Count)", "Rep2_log(Count)" )








#### Gene Type Count Comparison ###
simple_barplot_grouped <- function(data, label){
	#Remove duplicates
	data2 <-  data[!duplicated(data[c("Gene_Type")]),]
	data2 <- subset(data2, Gene_Type != "synthetic")
	category_sum <- aggregate(.~Category, data2[,c("Category", "Gene_Type_Count_Norm")], sum)
	category_sum <- category_sum[order(-category_sum$Gene_Type_Count_Norm),]
	category_sum$Category <- factor(category_sum$Category, levels = unique(category_sum$Category))
	pdf(file=paste(label, "Gene_Type_Normalized_Count_Min250.pdf",sep="_"),height=6,width=10,onefile=FALSE)
	print(ggplot(category_sum, aes( y=Gene_Type_Count_Norm, x=Category)) + 
    geom_bar(position="dodge", stat="identity")+
    scale_y_continuous(trans = 'log2')+
    theme_bw())
    dev.off()
}



simple_barplot_grouped(mouse.reshape, "mouse")




## PIE CHART
gene_type_pie_chart <- function(data) {
	gene_type_counts <- as.data.frame(table(data$Gene_Type))
	colnames(gene_type_counts) <- c("Gene_Type", "Count")
	gene_type_counts_withoutsequin <- subset(gene_type_counts, Gene_Type != "synthetic")
	gene_type_major <- subset(gene_type_counts_withoutsequin, Count > 250)
	gene_type_major$Category_New <- gene_type_major$Gene_Type
	gene_type_major$Gene_Type <- NULL
	gene_type_minor <- subset(gene_type_counts_withoutsequin, Count < 250)
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


mouse_gene_type_piechart <- gene_type_pie_chart(mouse.reshape)




# Basic piechart
pdf(file="Mouse_Gene_Type_Counts_BothRepMerged_PieChart_Min250.pdf",height=10,width=10,onefile=FALSE)
print(ggplot(mouse_gene_type_piechart, aes(x="", y=prop, fill=Category_New)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() + 
  #theme(legend.position="none") +
 #geom_text_repel(aes(y = ypos, label = Category_New), color = "black", size=6)+
    geom_text_repel(aes(y = ypos, label = paste(Category_New,percent(Count/sum(mouse_gene_type_piechart$Count)))), color = "black", size=5))
dev.off()






dotplot <- function(data, label) {
		data2 <- subset(data , Gene_Type_Count_Norm > 2)
		data2 <- subset(data2, Gene_Type!="synthetic")
		data2 <- data2[order(-data2$Gene_Type_Count_Norm),]
	   data2$Gene_Type <- factor(data2$Gene_Type, levels = unique(data2$Gene_Type))
		pdf(file=paste(label, "Tail_Length_Distribution_All_Genes_Dotplot.pdf",sep="_"),height=5,width=12,onefile=FALSE)
			print(ggplot(data2, aes(x=Gene_Type, y=tail_length)) + 
				geom_quasirandom(varwidth = TRUE, aes())+
				geom_boxplot(aes(alpha=0), outlier.shape=NA)+
				stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
            		geom = "crossbar", width = 0.7, color="#c06c84")+
				stat_n_text(size=5)+
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
	dotplot(mouse.reshape,"Mouse_Nuclear")


dotplot_median <- function(data, label) {
	    data2 <- subset(data , Gene_Type_Count_Norm > 2)
		data2 <- subset(data2 , Gene_Count > 30)
		data2 <- subset(data2, Gene_Type!="synthetic")
		data3 <-  data2[!duplicated(data2[c("Gene_Name")]),]
		data3 <- data3[order(-data3$Gene_Type_Count_Norm),]
	    data3$Gene_Type <- factor(data3$Gene_Type, levels = unique(data3$Gene_Type))
		pdf(file=paste(label, "Tail_Length_Distribution_All_Genes_MEDIAN_Dotplot_Min30.pdf",sep="_"),height=5,width=12,onefile=FALSE)
			print(ggplot(data3, aes(x=Gene_Type, y=Median_Length, color=Gene_Type)) + 
				geom_quasirandom(varwidth = TRUE, aes())+
				geom_boxplot(aes(alpha=0), outlier.shape=NA)+
				stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
            		geom = "crossbar", width = 0.7, color="#c06c84")+
				stat_n_text(size=5)+
				theme_bw()+
				ggtitle(label)+
				xlab("Group")+
              	ylab("Median Tail length") +
				theme(axis.text=element_text(size=14),strip.text = element_text(size=13),
            		axis.title=element_text(size=17,face="bold"),
            		legend.title = element_text(size = 20),
            		legend.text = element_text(color = "black", size=15)))
		dev.off()
	}

	#dotplot(data_tails,"Rep_Merged")
	dotplot_median(mouse.reshape,"Mouse_Nuclear")














##Overall Tail comparison (Single transcript)
tail_comparison_overall <- function(data, label) {
  data2 <- subset(data , Gene_Count > 30)
  data2 <- subset(data2, Gene_Type!="synthetic")
  data_means <- aggregate(.~Sample+Gene_Type, data2[,c("Sample", "Median_Length", "Gene_Type")], median)
  pdf(file=paste(label, "Overall_Tail_Comparison_Single_Transcript.pdf",sep="_"),height=6,width=10,onefile=FALSE)
  print(ggplot(data2, aes(x=tail_length, color=Sample)) +
    geom_density()+
    theme_bw()+
    #xlim(-10, 350)+
    facet_wrap(~Gene_Type, scales="free"))
  dev.off()
}


tail_comparison_overall(mouse.reshape,"Mouse_Nuclear")

























### Per-gene count / PolYA tail analysis
	per_gene_analysis <- function(data) {
	data <- subset(data, Gene_Count > 20)
	data$Group <- data$tail_length
	data$Group[which(data$tail_length == 0)] <- "No_PolyA"
	data$Group[which(data$tail_length < 10 & data$tail_length  > 0)]<- "Small_PolyA"
	data$Group[which(data$tail_length > 10 )]<- "Long_PolyA"
	return(data)
}


mouse_per_gene_tails <- per_gene_analysis(mouse.reshape)







count_polyA_population <- function(data) {
data<- subset(data, Gene_Count >10)
data$Tail_Class <- data$tail_length
data$Tail_Class[which(data$tail_length == 0)] <- "No_PolyA"
data$Tail_Class[which(data$tail_length < 10 & data$tail_length  > 0)]<- "Small_PolyA"
data$Tail_Class[which(data$tail_length >= 10 )]<- "PolyA"
final <- vector()
pre <- vector()
for (gene in unique(data$Gene_Name)) {
  subs <- subset(data, Gene_Name==gene)
  sub_table <- table(subs$Tail_Class)
  pre$gene <- gene
  pre <- as.data.frame(pre)
  pre$PolyA_Count <- as.numeric(sub_table["PolyA"])
  pre$No_PolyA_Count <- as.numeric(sub_table["No_PolyA"])
  final<- rbind(final,pre)
}

final[["PolyA_Count"]][is.na(final[["PolyA_Count"]])] <- 0
final[["No_PolyA_Count"]][is.na(final[["No_PolyA_Count"]])] <- 0
final$Total <- final$PolyA_Count + final$No_PolyA_Count
final$Ratio_PolyA <- final$PolyA_Count /final$Total
final$Ratio_NoPolyA<- final$No_PolyA_Count /final$Total
merged <- merge(data, final, by.x="Gene_Name", by.y="gene")
  return(merged)
}

mouse.tail_class_count <- count_polyA_population(mouse.reshape)




### Count pA on mtrRNAs

mouse_mtrRNA<- subset(mouse.tail_class_count, Gene_Type=="Mt_rRNA")

mouse_mtrRNA_uniq <-mouse_mtrRNA[!duplicated(mouse_mtrRNA[c("Gene_Name", "Sample")]),]









export_reads <- function(data,label) {
	for (group in unique(data$Group)) {
		subs <- subset(data, Group==group)
		subs2 <- subs[,c("Read_ID")]
		write.table(subs2, file=paste(group, label,"read_id.tsv",sep="_"), sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
	}
}

export_reads(mouse_per_gene_tails,"Mouse")








#### Sequin Expected vs Observed

sequins <- subset(mouse.unique, Gene_Type=="synthetic")
sequins_obs <- sequins[,c("Gene_Name", "Gene_Count", "Gene_Count_Norm")]

sequins_exp <- read.delim("rnasequin_genes_2.4.tsv")
sequins_exp2 <- sequin_exp[,c("NAME", "LENGTH", "MIX_A" )]


sequins_exp_obs <- merge(sequins_obs,sequins_exp2, by.x="Gene_Name", by.y="NAME" )

sequins_exp_obs$Log_Observed_Counts <- log(sequins_exp_obs$Gene_Count+1)
sequins_exp_obs$Log_Expected_Counts <- log(sequins_exp_obs$MIX_A+1)

sequins_exp_obs$Observed_Counts_Scaled <- rescale(sequins_exp_obs$Log_Observed_Counts)
sequins_exp_obs$Expected_Counts_Scaled <- rescale(sequins_exp_obs$Log_Expected_Counts)




plot_denscols_with_corr_pearson<-function(pdfname,my_x,my_y,xlab,ylab) {
  pdf(file=paste(pdfname, "_pearson.pdf",sep=""), height=6, width=6)
  dcols<-densCols(my_x,my_y, colramp=colorRampPalette(blues9[-(1:3)]))
  plot(my_x,my_y,col=dcols,cex=1, cex.lab=1,cex.main=3,lwd=5,pch=20,xlab=xlab,ylab=ylab)
  title(main=pdfname, col.main="black", font.main=4)
  #abline(0,1,lwd=3)
  # Correlation
  test<-cor.test(my_x,my_y, method="pearson")
  print(test)
  cor222<-paste("Pearson's r =",round(as.numeric(test$estimate),3))
  #pval<-paste("Pval =",test$p.value)
  mtext(paste(cor222))
  #mtext(paste(cor222,pval,sep=" ; ")) #Print the subtitle with the dataset correlation
  dev.off()
}

plot_denscols_with_corr_pearson("Sequins_Exp_vs_Obs_LogCounts", sequins_exp_obs$Log_Expected_Counts, sequins_exp_obs$Log_Observed_Counts, "log(Expected Counts)", "log(Observed Counts)" )

plot_denscols_with_corr_pearson("Sequins_Exp_vs_Obs_ScaledLogCounts", sequins_exp_obs$Expected_Counts_Scaled, sequins_exp_obs$Observed_Counts_Scaled, "scaled_log(Expected Counts)", "scaled_log(Observed Counts)" )





```




```bash


### EXTRACT THE BAM FOR SMALL RNA READS
java -jar /users/enovoa/boguzhan/Software/picard/build/libs/picard.jar FilterSamReads \
       I=cDNA964321_porechop_genome38_sequin_rrna.sorted.bam \
       O=cDNA964321_porechop_genome38_sequin_rrna.NoPolyA.bam\
       READ_LIST_FILE=No_PolyA_Mouse_read_id.tsv\
       FILTER=includeReadList

samtools sort cDNA964321_porechop_genome38_sequin_rrna.NoPolyA.bam cDNA964321_porechop_genome38_sequin_rrna.NoPolyA.sorted
samtools index cDNA964321_porechop_genome38_sequin_rrna.NoPolyA.sorted.bam



java -jar /users/enovoa/boguzhan/Software/picard/build/libs/picard.jar FilterSamReads \
       I=cDNA964321_porechop_genome38_sequin_rrna.sorted.bam \
       O=cDNA964321_porechop_genome38_sequin_rrna.LongPolyA.bam\
       READ_LIST_FILE=Long_PolyA_Mouse_read_id.tsv\
       FILTER=includeReadList

samtools sort cDNA964321_porechop_genome38_sequin_rrna.LongPolyA.bam cDNA964321_porechop_genome38_sequin_rrna.LongPolyA.sorted
samtools index cDNA964321_porechop_genome38_sequin_rrna.LongPolyA.sorted.bam














