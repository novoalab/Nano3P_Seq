########################################################
######## ZEBRAFISH 
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

################################################################################
################################################################################
######### PROCESS NANO3PSEQ DATA ###############################################
################################################################################


#Import tail
ribodepleted_rep1.tails <- read.delim("cDNA786327_tails.csv", sep=",")
ribodepleted_rep2.tails <- read.delim("cDNA123791_tails.csv", sep=",")

polyA.tails <- read.delim("cDNA8523612_tails.csv", sep=",")


manipulate_tail_nano3pseq<- function(data) { 
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

ribodepleted_rep1.tails_processed <- manipulate_tail_nano3pseq(ribodepleted_rep1.tails)
ribodepleted_rep2.tails_processed <- manipulate_tail_nano3pseq(ribodepleted_rep2.tails)

polyA.tails_processed <- manipulate_tail_nano3pseq(polyA.tails)

ribodepleted_merged.tails_processed <- rbind(ribodepleted_rep1.tails_processed,ribodepleted_rep2.tails_processed)

#POLYA SELECTED
polyA_hpf4.data <- read.delim("4hpf_pAselected.genome11_sequin_ALLRNAs_Merged.bed", header=FALSE)



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



# Reshape the tables and remove low quality reads
reshape_nano3pseq<- function(data,tails,label) {
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


ribodep_hpf2_merged.reshape <- reshape_nano3pseq(ribodep_hpf2_merged.data,ribodepleted_merged.tails_processed,"Ribodep_2hpf_merged")
ribodep_hpf4_merged.reshape <- reshape_nano3pseq(ribodep_hpf4_merged.data,ribodepleted_merged.tails_processed,"Ribodep_4hpf_merged")
ribodep_hpf6_merged.reshape <- reshape_nano3pseq(ribodep_hpf6_merged.data,ribodepleted_merged.tails_processed,"Ribodep_6hpf_merged")

polyA_hpf4.reshape <- reshape(polyA_hpf4.data,polyA.tails_processed,"PolyA_4hpf")


ribodep_all_merged <- rbind(ribodep_hpf2_merged.reshape, ribodep_hpf4_merged.reshape, ribodep_hpf6_merged.reshape)




## EXPORT DATA ###
write.table(ribodep_hpf2_merged.reshape, file="Nano3pseq_2hpf_merged_processed.tsv", sep="\t", quote=FALSE, row.names=FALSE)
write.table(ribodep_hpf4_merged.reshape, file="Nano3pseq_4hpf_merged_processed.tsv", sep="\t", quote=FALSE, row.names=FALSE)
write.table(ribodep_hpf6_merged.reshape, file="Nano3pseq_6hpf_merged_processed.tsv", sep="\t", quote=FALSE, row.names=FALSE)

write.table(polyA_hpf4.reshape, file="Nano3pseq_4hpf_pA_merged_processed.tsv", sep="\t", quote=FALSE, row.names=FALSE)


ribodep_hpf2_merged.reshape <- read.delim("Nano3pseq_2hpf_merged_processed.tsv")
ribodep_hpf4_merged.reshape <- read.delim("Nano3pseq_4hpf_merged_processed.tsv")
ribodep_hpf6_merged.reshape <- read.delim("Nano3pseq_6hpf_merged_processed.tsv")


ribodep_hpf2_merged.mrna <- subset(ribodep_hpf2_merged.reshape, Gene_Type=="protein_coding")
ribodep_hpf4_merged.mrna <- subset(ribodep_hpf4_merged.reshape, Gene_Type=="protein_coding")
ribodep_hpf6_merged.mrna <- subset(ribodep_hpf6_merged.reshape, Gene_Type=="protein_coding")













# 4. ADD TAIL INFORMATION OF 3 GROUPS
#Order the groups
gr_groups <- read.delim("gr_list.txt")



ThreeGroups_Tails <- function(hpf2, hpf4, hpf6, cov_th) {
  #Unique Columns
  hpf2.uniq <- hpf2[!duplicated(hpf2[c("Gene_Name")]),]
  hpf4.uniq <- hpf4[!duplicated(hpf4[c("Gene_Name")]),]
  hpf6.uniq <- hpf6[!duplicated(hpf6[c("Gene_Name")]),]
  #Coverage Threshold
  hpf2.mincov <- subset(hpf2.uniq, Gene_Count > cov_th)
  hpf4.mincov <- subset(hpf4.uniq, Gene_Count > cov_th)
  hpf6.mincov <- subset(hpf6.uniq, Gene_Count > cov_th)
  #Choose columns
  hpf2.mincov2 <- hpf2.mincov[,c("Gene_Name", "Median_Length")]
  hpf4.mincov2 <- hpf4.mincov[,c("Gene_Name", "Median_Length")]
  hpf6.mincov2 <- hpf6.mincov[,c("Gene_Name", "Median_Length")]
  #Merge the tables
  hpf24 <- merge(hpf2.mincov2, hpf4.mincov2, all.x=TRUE,all.y=TRUE, by.x="Gene_Name", by.y="Gene_Name")
  hpf246 <- merge(hpf24, hpf6.mincov2,all.x=TRUE,all.y=TRUE, by.x="Gene_Name", by.y="Gene_Name")
  #Column names
  colnames(hpf246) <- c("Gene_Name", "Median_Length.2hpf", "Median_Length.4hpf", "Median_Length.6hpf")
  #Merge with Groups
  hpf246_2 <- merge(hpf246, gr_groups,all.x=TRUE,  by.x="Gene_Name", "Gene_ID")
  hpf246_2$Group[is.na(hpf246_2$Group)] <- "Rest"
  hpf246_final <- melt(hpf246_2)
  hpf246_final <- hpf246_final[complete.cases(hpf246_final), ]
return(hpf246_final)
}


Nano3P_threegroups <- ThreeGroups_Tails(ribodep_hpf2_merged.mrna, ribodep_hpf4_merged.mrna,ribodep_hpf6_merged.mrna,20)




boxplots <- function(data, label) {
my_comparisons <- list( c("Median_Length.2hpf", "Median_Length.4hpf"), c("Median_Length.2hpf", "Median_Length.6hpf"), c("Median_Length.4hpf", "Median_Length.6hpf") )
pdf(file=paste(label, "Zebrafish_Embryos_TailLength_Boxplot_Maternal_vs_Rest_3Groups_log.pdf",sep="_"),height=5,width=10,onefile=FALSE)
    print(ggplot(data, aes(x = variable, y = log(value+1) )) + 
      geom_boxplot(aes(fill = variable),position = position_dodge(0.9)) +
      ylab("log2(MedianTail)")+
      stat_compare_means(comparisons = my_comparisons, label.y = c(5, 5.5, 6))+
       facet_wrap(~Group,nrow=1)+
      stat_n_text() + 
      #scale_fill_manual(values=colors)+
      #coord_cartesian(ylim = c(0,175))+
      theme_bw())
  dev.off()

pdf(file=paste(label, "Zebrafish_Embryos_TailLength_Boxplot_Maternal_vs_Rest_3Groups.pdf",sep="_"),height=5,width=10,onefile=FALSE)
    print(ggplot(data, aes(x = variable, y = value )) + 
      geom_boxplot(aes(fill = variable),position = position_dodge(0.9)) +
      ylab("MedianTail")+
      stat_compare_means(comparisons = my_comparisons, label.y = c(150, 160, 170))+
       facet_wrap(~Group,nrow=1)+
      stat_n_text() + 
      #scale_fill_manual(values=colors)+
      #coord_cartesian(ylim = c(0,175))+
      theme_bw())
  dev.off()
}



boxplots(Nano3P_threegroups, "Nano3PSeq")





#### PER READ
ThreeGroups_Tails_Perread <- function(hpf2, hpf4, hpf6) {
  #Unique Columns
  #hpf2.uniq <- hpf2[!duplicated(hpf2[c("Gene_Name")]),]
  #hpf4.uniq <- hpf4[!duplicated(hpf4[c("Gene_Name")]),]
  #hpf6.uniq <- hpf6[!duplicated(hpf6[c("Gene_Name")]),]
  #Choose columns
  hpf2_2 <- hpf2[,c("Gene_Name", "Sample", "tail_length")]
  hpf4_2 <- hpf4[,c("Gene_Name", "Sample", "tail_length")]
  hpf6_2 <- hpf6[,c("Gene_Name", "Sample", "tail_length")]
  #Samples
  hpf246 <- rbind(hpf2_2,hpf4_2,hpf6_2 )
  #Merge with Groups
  hpf246_2 <- merge(hpf246, gr_groups,all.x=TRUE,  by.x="Gene_Name", "Gene_ID")
  hpf246_2$Group[is.na(hpf246_2$Group)] <- "Rest"
return(hpf246_2)
}


Nano3P_threegroups <- ThreeGroups_Tails_Perread(ribodep_hpf2_merged.mrna, ribodep_hpf4_merged.mrna,ribodep_hpf6_merged.mrna)




boxplots_perread <- function(data, label) {
my_comparisons <- list( c("Ribodep_2hpf_merged", "Ribodep_4hpf_merged"), c("Ribodep_2hpf_merged", "Ribodep_6hpf_merged"), c("Ribodep_4hpf_merged", "Ribodep_6hpf_merged") )
pdf(file=paste(label, "Zebrafish_PERREAD_TailLength_Boxplot_Maternal_vs_Rest_3Groups_log.pdf",sep="_"),height=5,width=10,onefile=FALSE)
    print(ggplot(data, aes(x = Sample, y = log(tail_length+1) )) + 
      geom_boxplot(aes(fill = Sample),position = position_dodge(0.9)) +
      ylab("log2(Tail)")+
      stat_compare_means(comparisons = my_comparisons, label.y = c(5, 5.5, 6))+
       facet_wrap(~Group,nrow=1)+
      stat_n_text() + 
      #scale_fill_manual(values=colors)+
      #coord_cartesian(ylim = c(0,175))+
      theme_bw())
  dev.off()

pdf(file=paste(label, "Zebrafish_PERREAD_TailLength_Boxplot_Maternal_vs_Rest_3Groups.pdf",sep="_"),height=5,width=10,onefile=FALSE)
    print(ggplot(data, aes(x = Sample, y = tail_length)) + 
      geom_boxplot(aes(fill = Sample),position = position_dodge(0.9)) +
      ylab("Tail")+
      stat_compare_means(comparisons = my_comparisons, label.y = c(300 , 330, 360))+
       facet_wrap(~Group,nrow=1)+
      stat_n_text() + 
      #scale_fill_manual(values=colors)+
      #coord_cartesian(ylim = c(0,175))+
      theme_bw())
  dev.off()
}



boxplots_perread(Nano3P_threegroups, "Nano3PSeq")













# 5. ADD GENE COUNTS INFORMATION OF 3 GROUPS

three_groups_counts <- function(hpf2, hpf4, hpf6) {
  # Take the unique rows
  hpf2.uniq <-hpf2[!duplicated(hpf2[c("Gene_Name")]),]
  hpf4.uniq <-hpf4[!duplicated(hpf4[c("Gene_Name")]),]
  hpf6.uniq <-hpf6[!duplicated(hpf6[c("Gene_Name")]),]
  #Relevant columns
  hpf2.uniq <- hpf2.uniq[,c("Gene_Name", "Gene_Count_Norm")]
  hpf4.uniq <- hpf4.uniq[,c("Gene_Name", "Gene_Count_Norm")]
  hpf6.uniq <- hpf6.uniq[,c("Gene_Name", "Gene_Count_Norm")]
  # Merge the timepoints
  hpf24 <- merge(hpf2.uniq, hpf4.uniq, all.x=TRUE, all.y=TRUE, by.x="Gene_Name", by.y="Gene_Name", suffixes = c(".2hpf", ".4hpf"))
  hpf246 <- merge(hpf24, hpf6.uniq, all.x=TRUE, all.y=TRUE, by.x="Gene_Name", by.y="Gene_Name")
  #Rename Columns
  colnames(hpf246) <- c(colnames(hpf24),"Gene_Count_Norm.6hpf")
  #Merge the three groups with the data
  hpf246_2 <- merge(hpf246, gr_groups, all.x=TRUE, all.y=TRUE, by.x="Gene_Name", by.y="Gene_ID")
  hpf246_2$Group[is.na(hpf246_2$Group)] <- "Rest"
  #Rearrange the columns
  hpf246_3 <- hpf246_2[,c("Gene_Name","Group", "Gene_Count_Norm.2hpf","Gene_Count_Norm.4hpf","Gene_Count_Norm.6hpf" )]
  #Melt the data 
  hpf246_melted <- melt(hpf246_3)
  return(hpf246_melted)
}

Nano3P_threegroups_counts <- three_groups_counts(ribodep_hpf2_merged.mrna, ribodep_hpf4_merged.mrna,ribodep_hpf6_merged.mrna)



#median_counts <- aggregate(.~Group2+variable, ribodep_all_unique_246hpf_taillengths_melted_mrna_notrelativeto2hpf[,c("Group2", "variable", "value")], median)

boxplots_nano3Pseq <- function(data, label) {
my_comparisons <- list( c("Gene_Count_Norm.2hpf", "Gene_Count_Norm.4hpf"), c("Gene_Count_Norm.2hpf", "Gene_Count_Norm.6hpf"), c("Gene_Count_Norm.4hpf", "Gene_Count_Norm.6hpf") )
pdf(file=paste(label, "Zebrafish_EmbryosCounts_Boxplot_Maternal_vs_Rest_3Groups.pdf",sep="_"),height=5,width=10,onefile=FALSE)
    print(ggplot(data, aes(x = variable, y = log(value+1))) + 
      geom_boxplot(aes(fill = variable),position = position_dodge(0.9), outlier.shape=NA) +
      ylab("log2(Count)")+
      stat_compare_means(comparisons = my_comparisons, label.y = c(1.4 , 1.6, 1.8))+
      facet_wrap(~Group,nrow=1)+
      stat_n_text() + 
      #scale_fill_manual(values=colors)+
      #coord_cartesian(ylim = c(0,1.5))+
      theme_bw())
  dev.off()
}

boxplots_nano3Pseq(Nano3P_threegroups_counts, "Nano3P_Seq")






boxplots_nano3Pseqcut <- function(data, label) {
my_comparisons <- list( c("Gene_Count_Norm.2hpf", "Gene_Count_Norm.4hpf"), c("Gene_Count_Norm.2hpf", "Gene_Count_Norm.6hpf"), c("Gene_Count_Norm.4hpf", "Gene_Count_Norm.6hpf") )
pdf(file=paste(label, "Zebrafish_EmbryosCounts_Boxplot_Maternal_vs_Rest_3Groups_Cut.pdf",sep="_"),height=5,width=10,onefile=FALSE)
    print(ggplot(data, aes(x = variable, y = log(value+1))) + 
      geom_boxplot(aes(fill = variable),position = position_dodge(0.9), outlier.shape=NA) +
      ylab("log2(Count)")+
      stat_compare_means(comparisons = my_comparisons, label.y = c(1.4 , 1.6, 1.8))+
      facet_wrap(~Group,nrow=1)+
      stat_n_text() + 
      #scale_fill_manual(values=colors)+
      coord_cartesian(ylim = c(0,2))+
      theme_bw())
  dev.off()
}

boxplots_nano3Pseqcut(Nano3P_threegroups_counts, "Nano3P_Seq")













