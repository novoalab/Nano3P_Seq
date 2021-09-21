########################################################
######## ANALYSIS OF THE ZEBRAFISH pA vs Ribodep ####
########################################################
########## OGUZHAN BEGIK APRIL 2020 ####################

```R
library(stringr)
library(ggplot2)
library(ggbeeswarm)
library(ggrepel)
library(scales)
library(dplyr)
library(ggridges)

#Import tail
ribodepleted_rep1.tails <- read.delim("cDNA786327_tails.csv", sep=",")
ribodepleted_rep2.tails <- read.delim("cDNA123791_tails.csv", sep=",")

polyA.tails <- read.delim("cDNA8523612_tails.csv", sep=",")


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

polyA.tails_processed <- manipulate_tail_(polyA.tails)


ribodepleted_merged.tails_processed <- rbind(ribodepleted_rep1.tails_processed,ribodepleted_rep2.tails_processed)



# Import the data
#RIBODEPLETED DATA
ribodep_hpf4_rep1.data <- read.delim("4hpf.genome11_sequin_ALLRNAs_Merged_Rep1.bed", header=FALSE)

ribodep_hpf4_rep2.data <- read.delim("4hpf.genome11_sequin_ALLRNAs_Merged_Rep2.bed", header=FALSE)

#Merge Reps
ribodep_hpf4_merged.data<- rbind(ribodep_hpf4_rep1.data ,ribodep_hpf4_rep2.data )


#POLYA SELECTED
polyA_hpf4.data <- read.delim("4hpf_pAselected.genome11_sequin_ALLRNAs_Merged.bed", header=FALSE)





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
  merged4$Sample <- label
  return(merged4)
}



ribodep_hpf4_rep1.reshape <- reshape(ribodep_hpf4_rep1.data,ribodepleted_rep1.tails_processed,"Ribodep_4hpf_merged")
ribodep_hpf4_rep2.reshape <- reshape(ribodep_hpf4_rep2.data,ribodepleted_rep2.tails_processed,"Ribodep_4hpf_merged")


ribodep_hpf4_merged.reshape <- reshape(ribodep_hpf4_merged.data,ribodepleted_merged.tails_processed,"Ribodep_4hpf_merged")


polyA_hpf4.reshape <- reshape(polyA_hpf4.data,polyA.tails_processed,"PolyA_4hpf")




polyA_vs_ribodep_4hpf <- rbind(ribodep_hpf4_merged.reshape,polyA_hpf4.reshape)






#### Gene Type Count Comparison ###
simple_barplot_grouped <- function(data, label){
  #Remove duplicates
  data2 <-  data[!duplicated(data[c("Gene_Type", "Sample")]),]
  data2 <- subset(data2, Gene_Type != "synthetic" & Gene_Type != "ensembl_havana"  & Gene_Type != "ensembl" & Gene_Type != "havana" & Gene_Type != "polymorphic_pseudogene" & Gene_Type != "TEC" & Gene_Type != "unprocessed_pseudogene" & Gene_Type != "transcribed_unprocessed_pseudogene" & Gene_Type != "snRNA"& Gene_Type != "TR_J_gene")
  category_sum <- aggregate(.~Gene_Type+Sample, data2[,c("Gene_Type","Sample",  "Gene_Type_Count_Norm")], sum)
  category_sum <- category_sum[order(-category_sum$Gene_Type_Count_Norm),]
  category_sum$Gene_Type <- factor(category_sum$Gene_Type, levels = unique(category_sum$Gene_Type))
  pdf(file=paste(label, "Gene_Type_Normalized_Count.pdf",sep="_"),height=6,width=15,onefile=FALSE)
  print(ggplot(category_sum, aes(fill=Sample,  y=log(Gene_Type_Count_Norm+1), x=Gene_Type)) + 
    geom_bar(position="dodge", stat="identity")+
    #scale_y_continuous(trans = 'log2')+
    theme_bw())
    dev.off()
}



simple_barplot_grouped(polyA_vs_ribodep_4hpf, "PolyA_vs_Ribodep")








##Overall Tail comparison (Single transcript)
tail_comparison_overall <- function(data, label) {
  data2 <- subset(data, Gene_Type=="protein_coding")
  data_means <- aggregate(.~Sample+Gene_Type, data2[,c("Sample", "Mean_Length", "Gene_Type")], mean)
  print(data_means)
  pdf(file=paste(label, "Overall_Tail_Comparison_Single_Transcript.pdf",sep="_"),height=3,width=6,onefile=FALSE)
  print(ggplot(data2, aes(x=tail_length, color=Sample)) +
    geom_density()+
     stat_compare_means(method = "wilcox")+
    theme_bw()+
    #xlim(-10, 350)+
    facet_wrap(~Gene_Type, scales="free"))
  dev.off()
}


tail_comparison_overall(polyA_vs_ribodep_4hpf, "PolyA_vs_Ribodep_4HPF")









  ribodep_uniq <-  ribodep_hpf4_merged.reshape[!duplicated(ribodep_hpf4_merged.reshape[c("Gene_Name")]),]
   polya_uniq <-  polyA_hpf4.reshape[!duplicated(polyA_hpf4.reshape[c("Gene_Name")]),]


wilcox_test_data <- merge(ribodep_uniq,polya_uniq, by.x="Gene_Name", by.y="Gene_Name")
wilcox_test_data_mrna <- subset(wilcox_test_data, Gene_Type.x=="protein_coding")
wilcox_test_data_mrna2 <-  wilcox_test_data_mrna[,c("Gene_Name","Median_Length.x","Median_Length.y")]


res <- wilcox.test(wilcox_test_data_mrna2$Median_Length.x, wilcox_test_data_mrna2$Median_Length.y)








## Tail comparison median per gene
library(tidyverse)

tail_comparison_mean_per_gene_protein<- function(data, label) {
  data2 <- subset(data, Gene_Type=="protein_coding")
  data2 <-  data2[!duplicated(data2[c("Gene_Name", "Sample")]),]
  data_means <- aggregate(.~Sample, data2[,c("Sample", "Mean_Length")], mean)
  print(data_means)
  pdf(file=paste(label, "Mean_Tail_Per_Gene_Comparison_mRNA.pdf",sep="_"),height=6,width=10,onefile=FALSE)
  print(ggplot(data2, aes(x=Mean_Length, color=Sample)) +
    geom_density()+
      theme_bw()+
      coord_cartesian(xlim=c(-10, 200)))
  dev.off()
}



tail_comparison_mean_per_gene_protein(polyA_vs_ribodep_4hpf, "PolyA_vs_Ribodep_4HPF")





ribodep_hpf4_rep1.unique <-  ribodep_hpf4_rep1.reshape[!duplicated(ribodep_hpf4_rep1.reshape[c("Gene_Name", "Sample")]),]
ribodep_hpf4_rep2.unique <-  ribodep_hpf4_rep2.reshape[!duplicated(ribodep_hpf4_rep2.reshape[c("Gene_Name", "Sample")]),]



ribodep_hpf4_merged.unique <-  ribodep_hpf4_merged.reshape[!duplicated(ribodep_hpf4_merged.reshape[c("Gene_Name", "Sample")]),]
polyA_hpf4.unique <-  polyA_hpf4.reshape[!duplicated(polyA_hpf4.reshape[c("Gene_Name", "Sample")]),]



columns1 <- c("Gene_Type", "Gene_Name", "Mean_Length", "Gene_Count", "Gene_Count_Norm")
columns2 <- c("Gene_Type",  "Gene_Name")

polyA_vs_ribodep_4hpf_merged <- merge(ribodep_hpf4_merged.unique[,columns1] ,polyA_hpf4.unique[,columns1], by.x=columns2, by.y=columns2  )

ripodep_rep1_vs_rep2 <- merge(ribodep_hpf4_rep1.unique[,columns1] ,ribodep_hpf4_rep2.unique[,columns1], by.x=columns2, by.y=columns2  )


#SCATTER PLOT 
Scatter_Comparison_Gene_Count <- function(data, label){
   corr <- cor.test(data$Gene_Count_Norm.x, data$Gene_Count_Norm.y, method = "pearson", conf.level = 0.95)
   value <- as.numeric(corr$estimate)
    pdf(file=paste(label,"Normalized_Gene_Count_ScatterPlot_LogTransformated.pdf", sep="_"),height=4,width=4,onefile=FALSE)
    print(ggplot(data, aes(x=Gene_Count_Norm.x, y=Gene_Count_Norm.y)) + 
    theme_bw()+ 
    scale_x_continuous(trans = 'log2') +
    scale_y_continuous(trans = 'log2')+
   geom_text_repel(data=subset(data, abs(Gene_Count_Norm.x-Gene_Count_Norm.y) > 15), aes(label=Gene_Type))+
    ggtitle(label)+
    xlab("Ribodepletion Normalized Gene Count")+
    ylab("PolyA Selection Normalized Gene Count")+
      geom_point())
      dev.off()
}
Scatter_Comparison_Gene_Count(polyA_vs_ribodep_4hpf_merged, "Ribodepletion_PolyA")
Scatter_Comparison_Gene_Count(ripodep_rep1_vs_rep2, "ripodep_rep1_vs_rep2")












polyA_vs_ribodep_4hpf_merged_protein <- subset(polyA_vs_ribodep_4hpf_merged, Gene_Type=="protein_coding")
polyA_vs_ribodep_4hpf_merged_protein2 <- subset(polyA_vs_ribodep_4hpf_merged_protein, Gene_Count.x>20 & Gene_Count.y>20 )

polyA_vs_ribodep_4hpf_merged_protein_40 <- subset(polyA_vs_ribodep_4hpf_merged_protein2, Mean_Length.x <= 40 & Mean_Length.y <= 40)


ripodep_rep1_vs_rep2_protein <- subset(ripodep_rep1_vs_rep2, Gene_Type=="protein_coding")
ripodep_rep1_vs_rep2_protein2 <- subset(ripodep_rep1_vs_rep2_protein, Gene_Count.x>20 & Gene_Count.y>20 )

ripodep_rep1_vs_rep2_protein_40 <- subset(ripodep_rep1_vs_rep2_protein2, Mean_Length.x <= 40 & Mean_Length.y <= 40)


plot_denscols_with_corr_pearson<-function(pdfname,my_x,my_y,xlab,ylab) {
  pdf(file=paste(pdfname, "_pearson.pdf",sep=""), height=6, width=6)
  dcols<-densCols(my_x,my_y, colramp=colorRampPalette(blues9[-(1:3)]))
  plot(my_x,my_y,col=dcols,cex=1, cex.lab=1,cex.main=3,lwd=5,pch=20,xlab=xlab,ylab=ylab)
  title(main=pdfname, col.main="black", font.main=4)
  abline(0,1,lwd=3)
  # Correlation
  test<-cor.test(my_x,my_y, method="pearson")
  print(test)
  cor222<-paste("Pearson's rho =",round(as.numeric(test$estimate),3))
  #pval<-paste("Pval =",test$p.value)
  mtext(paste(cor222))
  #mtext(paste(cor222,pval,sep=" ; ")) #Print the subtitle with the dataset correlation
  dev.off()
}

plot_denscols_with_corr_pearson("Ribodep_vs_polyA_Mean_Tail", polyA_vs_ribodep_4hpf_merged_protein2$Mean_Length.x, polyA_vs_ribodep_4hpf_merged_protein2$Mean_Length.y, "Ribodep_Mean_pA_Length", "PolyA_Mean_pA_Length" )

plot_denscols_with_corr_pearson("Rep1_vs_Rep2_Mean_Tail", ripodep_rep1_vs_rep2_protein2$Mean_Length.x, ripodep_rep1_vs_rep2_protein2$Mean_Length.y, "Rep1_Mean_pA_Length", "Rep2_Mean_pA_Length" )


plot_denscols_with_corr_pearson("Ribodep_vs_polyA_Gene_Count", polyA_vs_ribodep_4hpf_merged_protein2$Gene_Count_Norm.x, polyA_vs_ribodep_4hpf_merged_protein2$Gene_Count_Norm.y, "Ribodep_Gene_Count", "PolyA_Gene_Count" )

plot_denscols_with_corr_pearson("Rep1_vs_Rep2_Gene_Count", ripodep_rep1_vs_rep2_protein2$Gene_Count_Norm.x, ripodep_rep1_vs_rep2_protein2$Gene_Count_Norm.y, "Rep1_Gene_Count", "Rep2_Gene_Count" )







