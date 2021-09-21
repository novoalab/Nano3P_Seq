########################################################
######## ZEBRAFISH dRNA Tail  Comparison ##############
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
library(reshape2)
library(EnvStats)
library(ggpubr)

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

ribodepleted_merged.tails_processed <- rbind(ribodepleted_rep1.tails_processed,ribodepleted_rep2.tails_processed)

polyA.tails_processed <- manipulate_tail_(polyA.tails)


# Import the data
#RIBODEPLETED DATA
ribodep_hpf2_rep1.data <- read.delim("2hpf.genome11_sequin_ALLRNAs_Merged_Rep1.bed", header=FALSE)
ribodep_hpf4_rep1.data <- read.delim("4hpf.genome11_sequin_ALLRNAs_Merged_Rep1.bed", header=FALSE)
ribodep_hpf6_rep1.data <- read.delim("6hpf.genome11_sequin_ALLRNAs_Merged_Rep1.bed", header=FALSE)

ribodep_hpf2_rep2.data <- read.delim("2hpf.genome11_sequin_ALLRNAs_Merged_Rep2.bed", header=FALSE)
ribodep_hpf4_rep2.data <- read.delim("4hpf.genome11_sequin_ALLRNAs_Merged_Rep2.bed", header=FALSE)
ribodep_hpf6_rep2.data <- read.delim("6hpf.genome11_sequin_ALLRNAs_Merged_Rep2.bed", header=FALSE)

#Merge Reps
ribodep_hpf2_merged.data <- rbind(ribodep_hpf2_rep1.data ,ribodep_hpf2_rep2.data )
ribodep_hpf4_merged.data<- rbind(ribodep_hpf4_rep1.data ,ribodep_hpf4_rep2.data )
ribodep_hpf6_merged.data <- rbind(ribodep_hpf6_rep1.data ,ribodep_hpf6_rep2.data )

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


ribodep_hpf2_merged.reshape <- reshape(ribodep_hpf2_merged.data,ribodepleted_merged.tails_processed,"Ribodep_2hpf_merged")
ribodep_hpf4_merged.reshape <- reshape(ribodep_hpf4_merged.data,ribodepleted_merged.tails_processed,"Ribodep_4hpf_merged")
ribodep_hpf6_merged.reshape <- reshape(ribodep_hpf6_merged.data,ribodepleted_merged.tails_processed,"Ribodep_6hpf_merged")

ribodep_all_merged <- rbind(ribodep_hpf2_merged.reshape, ribodep_hpf4_merged.reshape, ribodep_hpf6_merged.reshape)

polyA_hpf4.reshape <- reshape(polyA_hpf4.data,polyA.tails_processed,"PolyA_4hpf")





### COMPARE WITH LITERATURE 
dRNA.data <- read.delim("dRNA_PerGeneMedian.tsv")

dRNA.data$dRNA_2hpf <-(dRNA.data$WT.2hpf.Rep1+dRNA.data$WT.2hpf.Rep2)/2
dRNA.data$dRNA_4hpf <-(dRNA.data$WT.4hpf.Rep1+dRNA.data$WT.4hpf.Rep2)/2
dRNA.data$dRNA_6hpf <-(dRNA.data$WT.6hpf.Rep1+dRNA.data$WT.6hpf.Rep2)/2







merging <- function(data_dRNA,data_us, timepoint, threshold) {
   data_dRNA2 <- data_dRNA[,c("Gene", timepoint)] 
   #Take unique gene names for our data
   data_us2 <-  data_us[!duplicated(data_us[c("Gene_Name", "Sample")]),]
   data_us3 <- subset(data_us2, Gene_Type =="protein_coding")
   data_us4 <- data_us3[,c("Gene_Name","Median_Length","Mean_Length","Gene_Count","Gene_Count_Norm","Sample")]
   comparison <- merge(data_us4, data_dRNA2, by.x="Gene_Name", by.y="Gene")
   comparison_filtered <- subset(comparison, Gene_Count > threshold)
   comparison_complete <-  comparison_filtered[complete.cases(comparison_filtered), ]
   return(comparison_complete)
  }

dRNA_hpf2_ribodep.processed <- merging(dRNA.data, ribodep_hpf2_merged.reshape,"dRNA_2hpf" ,30)
dRNA_hpf4_ribodep.processed <- merging(dRNA.data, ribodep_hpf4_merged.reshape,"dRNA_4hpf" ,30)
dRNA_hpf6_ribodep.processed <- merging(dRNA.data, ribodep_hpf6_merged.reshape,"dRNA_6hpf" ,30)

dRNA_hpf4_pA.processed <- merging(dRNA.data, polyA_hpf4.reshape,"dRNA_4hpf" ,30)




#SCATTER PLOT 
scatter_bartel_comparison <- function(data, label){
   corr <- cor.test(data$Median_Length, data$Median.T, method = "pearson", conf.level = 0.95)
   value <- as.numeric(corr$estimate)
      pdf(file=paste(label,"Bartel_vs_OurStudy_TailMedian_Min30.pdf", sep="_"),height=4,width=4,onefile=FALSE)
    print(ggplot(data, aes(x=Median_Length, y=Median.TL)) + 
    theme_bw()+ 
        ylim(0,150)+
        xlim(0,150)+
        ggtitle(label)+
        annotate(geom="text", x=50, y=150, label=paste("Pearson Correlation =", value),
              color="red", size=2)+
    xlab("NanoTailCapture (this work) Median Tail Length")+
        ylab("PAL-Seq (Subtelny et al 2014) Median Tail Length")+
      geom_point())
      dev.off()
}
scatter_bartel_comparison(bartel_hpf2.processed, "2HPF")
scatter_bartel_comparison(bartel_hpf4.processed, "4HPF")
scatter_bartel_comparison(bartel_hpf6.processed, "6HPF")







plot_denscols_with_corr_pearson<-function(pdfname,my_x,my_y,xlab,ylab) {
  pdf(file=paste(pdfname, "_pearson.pdf",sep=""), height=6, width=6)
  dcols<-densCols(my_x,my_y, colramp=colorRampPalette(blues9[-(1:3)]))
  plot(my_x,my_y,col=dcols,cex=1, cex.lab=1,cex.main=3,lwd=5,pch=20,xlab=xlab,ylab=ylab,xlim = c(0,160),ylim = c(0,160) )
  title(main=pdfname, col.main="black", font.main=4)
  abline(a=0,b=1)
  # Correlation
  test<-cor.test(my_x,my_y, method="pearson")
  print(test)
  cor222<-paste("Pearson's p =",round(as.numeric(test$estimate),3))
  #pval<-paste("Pval =",test$p.value)
  mtext(paste(cor222))
  #mtext(paste(cor222,pval,sep=" ; ")) #Print the subtitle with the dataset correlation
  dev.off()
}


plot_denscols_with_corr_pearson("dRNA_vs_Nano3PSeq_2HPF", dRNA_hpf2_ribodep.processed$Median_Length , dRNA_hpf2_ribodep.processed$dRNA_2hpf, "Nano3P_seq (this work) Median Tail Length 2HPF ", "dRNA Median Tail Length 2HPF" )


plot_denscols_with_corr_pearson("dRNA_vs_Nano3PSeq_4HPF", dRNA_hpf4_ribodep.processed$Median_Length , dRNA_hpf4_ribodep.processed$dRNA_4hpf, "Nano3P_seq (this work) Median Tail Length 4HPF ", "dRNA Median Tail Length 4HPF" )



plot_denscols_with_corr_pearson("dRNA_vs_Nano3PSeq_6HPF", dRNA_hpf6_ribodep.processed$Median_Length , dRNA_hpf6_ribodep.processed$dRNA_6hpf, "Nano3P_seq (this work) Median Tail Length 6HPF ", "dRNA Median Tail Length 6HPF" )


plot_denscols_with_corr_pearson("dRNA_vs_Nano3PSeq_4HPF_pA", dRNA_hpf4_pA.processed$Median_Length , dRNA_hpf4_pA.processed$dRNA_4hpf, "Nano3P_seq (this work) Median Tail Length 4HPF", "dRNA Median Tail Length 4HPF" )






### DOTPLOTS 


bartel_hpf2.processed$Time_point <- "2HPF"
bartel_hpf4.processed$Time_point <- "4HPF"
bartel_hpf6.processed$Time_point <- "6HPF"

comparison_all <- rbind(bartel_hpf2.processed,bartel_hpf4.processed,bartel_hpf6.processed )



  
    pdf(file= "Bartel_Dotplot.pdf",height=10,width=20,onefile=FALSE)
      print(ggplot(comparison_all, aes(x=Time_point, y=Median.TL)) + 
        geom_quasirandom(varwidth = TRUE, aes(color=Time_point))+
        geom_boxplot(aes(alpha=0), outlier.shape=NA)+
        stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                geom = "crossbar", width = 0.7, color="#c06c84")+
        theme_bw()+
        ylim(0,170)+
        #facet_wrap(~Time_point,nrow=1)+
        ggtitle("Zebrafish Embryo")+
        xlab("Time Points")+
                ylab("Tail length") +
        theme(axis.text=element_text(size=14),strip.text = element_text(size=13),
                axis.title=element_text(size=17,face="bold"),
                legend.title = element_text(size = 20),
                legend.text = element_text(color = "black", size=15)))
    dev.off()
  
    pdf(file= "Begik_Dotplot.pdf",height=10,width=20,onefile=FALSE)
      print(ggplot(comparison_all, aes(x=Time_point, y=Median_Length)) + 
        geom_quasirandom(varwidth = TRUE, aes(color=Time_point))+
        geom_boxplot(aes(alpha=0), outlier.shape=NA)+
        stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                geom = "crossbar", width = 0.7, color="#c06c84")+
        theme_bw()+
        ylim(0,170)+
        #facet_wrap(~Time_point,nrow=1)+
        ggtitle("Zebrafish Embryo")+
        xlab("Time Points")+
        ylab("Tail length") +
        theme(axis.text=element_text(size=14),strip.text = element_text(size=13),
                axis.title=element_text(size=17,face="bold"),
                legend.title = element_text(size = 20),
                legend.text = element_text(color = "black", size=15)))
    dev.off()















### 3 Groups dRNA



#### INVESTIGATE THE MATERNAL RNAs TAIL LENGTHS
GENE_ID <- read.delim("Zebrafish_ID_Name_Conversion.txt")


dRNA.data.complete <-  dRNA.data[complete.cases(dRNA.data), ]
dRNA.data.complete2 <- dRNA.data.complete[,c("Gene", "dRNA_2hpf","dRNA_4hpf","dRNA_6hpf")]




#Gene profile
maternal <- read.delim("allTranscripts_riboZero_rnaSeq.maternal.txt", header=FALSE)
maternal$Group <- "Maternal"
zyfir <- read.delim("allTranscripts_riboZero_rnaSeq.zyfir.txt", header=FALSE)
zyfir$Group <- "Zyfir"

mir430 <- read.delim("allTranscripts_riboZero_rnaSeq.mir430.txt", header=FALSE)
mir430$Group <- "mir430"

id_convert <- read.delim("Transcript_ID_To_Gene_ID.txt")

groups <- rbind(maternal, zyfir, mir430)
groups_id <- merge(groups,id_convert, by.x="V1", by.y="Transcript.stable.ID")
groups_id$V1 <- NULL


dRNA.data.complete3 <- merge(dRNA.data.complete2, groups_id, all.x=TRUE, all.y=TRUE, by.x="Gene", by.y="Gene.stable.ID")
dRNA.data.complete3 $Group[is.na(dRNA.data.complete3 $Group)] <- "Rest"

dRNA.data.melted <- melt(dRNA.data.complete3)



my_comparisons <- list( c("dRNA_2hpf", "dRNA_4hpf"), c("dRNA_2hpf", "dRNA_6hpf"), c("dRNA_4hpf", "dRNA_6hpf") )

pdf(file="Zebrafish_Embryos_TailLength_Boxplot_Maternal_vs_Rest_3Groups_log_dRNA.pdf",height=5,width=20,onefile=FALSE)
    print(ggplot(dRNA.data.melted, aes(x = variable, y = log(value+1) )) + 
      geom_boxplot(aes(fill = variable),position = position_dodge(0.9)) +
      ylab("log2(MedianTail)")+
      stat_compare_means(comparisons = my_comparisons, label.y = c(5, 5.5, 6))+
       facet_wrap(~Group,nrow=1)+
      #stat_n_text() + 
      #scale_fill_manual(values=colors)+
      #coord_cartesian(ylim = c(0,175))+
      theme_bw())
  dev.off()


pdf(file="Zebrafish_Embryos_TailLength_Boxplot_Maternal_vs_Rest_3Groups_dRNA.pdf",height=5,width=20,onefile=FALSE)
    print(ggplot(dRNA.data.melted, aes(x = variable, y = value )) + 
      geom_boxplot(aes(fill = variable),position = position_dodge(0.9)) +
      ylab("MedianTail")+
      stat_compare_means(comparisons = my_comparisons, label.y = c(150, 160, 170))+
       facet_wrap(~Group,nrow=1)+
      #stat_n_text() + 
      #scale_fill_manual(values=colors)+
      #coord_cartesian(ylim = c(0,175))+
      theme_bw())
  dev.off()











################## MAKING DENSITY PLOTS WITH Nano3P vs dRNA



merging_keepingNA <- function(data_dRNA,data_us, timepoint, threshold) {
   data_dRNA2 <- data_dRNA[,c("Gene", timepoint)] 
   #Take unique gene names for our data
   data_us2 <-  data_us[!duplicated(data_us[c("Gene_Name", "Sample")]),]
   data_us3 <- subset(data_us2, Gene_Type =="protein_coding")
   data_us4 <- data_us3[,c("Gene_Name","Median_Length","Mean_Length","Gene_Count","Gene_Count_Norm","Sample")]
   comparison <- merge(data_us4, data_dRNA2, all.x=TRUE, by.x="Gene_Name", by.y="Gene")
   comparison2 <- comparison[,c("Gene_Name", "Median_Length", "dRNA_2hpf")]
   comparison_filtered <- subset(comparison, Gene_Count > threshold)
   #comparison_complete <-  comparison_filtered[complete.cases(comparison_filtered), ]


   return(comparison_filtered)
  }

dRNA_hpf2_ribodep.processed <- merging_keepingNA(dRNA.data, ribodep_hpf2_merged.reshape,"dRNA_2hpf" ,30)
dRNA_hpf4_ribodep.processed <- merging_keepingNA(dRNA.data, ribodep_hpf4_merged.reshape,"dRNA_4hpf" ,30)
dRNA_hpf6_ribodep.processed <- merging_keepingNA(dRNA.data, ribodep_hpf6_merged.reshape,"dRNA_6hpf" ,30)




















  #### COMPARE ANNA'S 4HPF AND MINE 

dRNA_counts <- read.delim("RawGeneCounts_ZF-MZT.tsv")

dRNA_counts_REP2 <-  dRNA_counts[,c("Gene.Name", "WT.2hpf.Rep2", "WT.4hpf.Rep2", "WT.6hpf.Rep2")]

dRNA_counts_clean <- subset(dRNA_counts_REP2, Gene.Name!="__alignment_not_unique" & Gene.Name!="__no_feature"& Gene.Name!="__ambiguous"& Gene.Name!="__not_aligned"& Gene.Name!="__too_low_aQual")

dRNA_counts_4hpf_anna <- dRNA_counts_clean[,c("Gene.Name", "WT.4hpf.Rep2")]
dRNA_counts_4hpf_anna2 <- subset(dRNA_counts_4hpf_anna, WT.4hpf.Rep2>0)




# Our data
drna_hpf4.data <- read.delim("PDBN042841_wt_4hpf_rep2_nonrRNA.genome11_sequin_ALLRNAs_Merged.bed", header=FALSE)


# Reshape the tables and remove low quality reads
reshape<- function(data,label) {
  data2 <- data[,c("V1", "V4", "V5", "V6", "V16", "V17")]
  colnames(data2) <- c("Chr","Read_ID", "Quality", "Strand", "Gene_Name", "Gene_Type")
  data3 <- subset(data2, Quality > 30)
  #Merge the data with tails
  ourdata <- data3[,c("Gene_Name", "Quality")]
  ourdata2 <- do.call(data.frame,(aggregate(. ~Gene_Name, data =ourdata, FUN = function(ourdata) c( count = length(ourdata) ) )))
  colnames(ourdata2) <- c("Gene_Name", "Gene_Count")
  coverage <- sum(ourdata2$Gene_Count)
  #Merge the data with the stats tabls
  merged2 <- merge(data3, ourdata2, by.x=c("Gene_Name"), by.y=c("Gene_Name"))
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
  merged5 <- subset(merged4,Category=="protein_coding")
  merged6 <- merged5[,c("Gene_Name", "Gene_Count")]
  merged7 <-  merged6 [!duplicated(merged6 [c("Gene_Name")]),]
  return(merged7)
}


drna_hpf4.us <- reshape(drna_hpf4.data,"dRNA")




drna_merged <- merge(dRNA_counts_4hpf_anna2,drna_hpf4.us, by.x="Gene.Name", by.y="Gene_Name" )




plot_denscols_with_corr_pearson<-function(pdfname,my_x,my_y,xlab,ylab) {
  pdf(file=paste(pdfname, "_pearson.pdf",sep=""), height=6, width=6)
  dcols<-densCols(my_x,my_y, colramp=colorRampPalette(blues9[-(1:3)]))
  plot(my_x,my_y,col=dcols,cex=1, cex.lab=1,cex.main=3,lwd=5,pch=20,xlab=xlab,ylab=ylab)
  title(main=pdfname, col.main="black", font.main=4)
  abline(a=0,b=1)
  # Correlation
  test<-cor.test(my_x,my_y, method="pearson")
  print(test)
  cor222<-paste("Pearson's p =",round(as.numeric(test$estimate),3))
  #pval<-paste("Pval =",test$p.value)
  mtext(paste(cor222))
  #mtext(paste(cor222,pval,sep=" ; ")) #Print the subtitle with the dataset correlation
  dev.off()
}


plot_denscols_with_corr_pearson("dRNA_4hpf_GeneCounts_Anna_vs_Oz", log(drna_merged$WT.4hpf.Rep2+1), log(drna_merged$Gene_Count+1), "log Gene Counts Anna ", "log Gene Counts Oz" )

plot_denscols_with_corr_pearson("abs_dRNA_4hpf_GeneCounts_Anna_vs_Oz", drna_merged$WT.4hpf.Rep2, drna_merged$Gene_Count, "Gene Counts Anna ",  "Gene Counts Oz" )





