########################################################
######## ANALYSIS OF THE ZEBRAFISH pA dRNA vs Nano3Pseq ####
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
library(ggpubr)
library(EnvStats)

################################################################################
######### PROCESS dRNA DATA ####################################################
################################################################################
## Direct RNA Data
drna_hpf2.data <- read.delim("PRPN023899_wt_2hpf_rep2_nonrRNA.genome11_sequin_ALLRNAs_Merged.bed", header=FALSE)
drna_hpf4.data <- read.delim("PDBN042841_wt_4hpf_rep2_nonrRNA.genome11_sequin_ALLRNAs_Merged.bed", header=FALSE)
drna_hpf6.data <- read.delim("PRPN039928_wt_6hpf_rep2_nonrRNA.genome11_sequin_ALLRNAs_Merged.bed", header=FALSE)




drna_hpf2.tail <- read.delim("nanopolish_tail_drna_2hpf_rep2.tsv")
drna_hpf4.tail <- read.delim("nanopolish_tail_drna_4hpf_rep2.tsv")
drna_hpf6.tail <- read.delim("nanopolish_tail_drna_6hpf_rep2.tsv")




manipulate_tail_dRNA<- function(data) { 
  data2 <- subset(data, qc_tag =="PASS")
  data2$polya_length <- as.numeric(as.character(data2$polya_length))
  data2[["polya_length"]][is.na(data2[["polya_length"]])] <- 0
  data3 <- data2[,c("readname", "polya_length")]
  colnames(data3) <- c("readname", "tail_length")
  return(data3)
}

drna_hpf2.tail_processed <- manipulate_tail_dRNA(drna_hpf2.tail)
drna_hpf4.tail_processed <- manipulate_tail_dRNA(drna_hpf4.tail)
drna_hpf6.tail_processed <- manipulate_tail_dRNA(drna_hpf6.tail)




# Reshape the tables and remove low quality reads
reshape_dRNA <- function(data,tails,label) {
  data2 <- data[,c("V1", "V4", "V5", "V6", "V16", "V17")]
  colnames(data2) <- c("Chr","Read_ID", "Quality", "Strand", "Gene_Name", "Gene_Type")
  data3 <- subset(data2, Quality > 30)
  #Merge the data with tails
  merged <- merge(data3, tails, by.x="Read_ID",by.y="readname")
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



drna_hpf2.reshape <- reshape_dRNA(drna_hpf2.data,drna_hpf2.tail_processed,"dRNA_2hpf")
drna_hpf4.reshape <- reshape_dRNA(drna_hpf4.data,drna_hpf4.tail_processed,"dRNA_4hpf")
drna_hpf6.reshape <- reshape_dRNA(drna_hpf6.data,drna_hpf6.tail_processed,"dRNA_6hpf")



write.table(drna_hpf2.reshape, file="dRNA_2hpf_rep2_processed.tsv", sep="\t", quote=FALSE, row.names=FALSE)
write.table(drna_hpf4.reshape, file="dRNA_4hpf_rep2_processed.tsv", sep="\t", quote=FALSE, row.names=FALSE)
write.table(drna_hpf6.reshape, file="dRNA_6hpf_rep2_processed.tsv", sep="\t", quote=FALSE, row.names=FALSE)


drna_hpf2.reshape <- read.delim("dRNA_2hpf_rep2_processed.tsv")
drna_hpf4.reshape <- read.delim("dRNA_4hpf_rep2_processed.tsv")
drna_hpf6.reshape <- read.delim("dRNA_6hpf_rep2_processed.tsv")


drna_hpf2.mrna <- subset(drna_hpf2.reshape, Gene_Type=="protein_coding")
drna_hpf4.mrna <- subset(drna_hpf4.reshape, Gene_Type=="protein_coding")
drna_hpf6.mrna <- subset(drna_hpf6.reshape, Gene_Type=="protein_coding")


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

polyA_hpf4.reshape <- reshape_nano3pseq(polyA_hpf4.data,polyA.tails_processed,"PolyA_4hpf")


ribodep_all_merged <- rbind(ribodep_hpf2_merged.reshape, ribodep_hpf4_merged.reshape, ribodep_hpf6_merged.reshape)




## EXPORT DATA ###
write.table(ribodep_hpf2_merged.reshape, file="Nano3pseq_2hpf_merged_processed.tsv", sep="\t", quote=FALSE, row.names=FALSE)
write.table(ribodep_hpf4_merged.reshape, file="Nano3pseq_4hpf_merged_processed.tsv", sep="\t", quote=FALSE, row.names=FALSE)
write.table(ribodep_hpf6_merged.reshape, file="Nano3pseq_6hpf_merged_processed.tsv", sep="\t", quote=FALSE, row.names=FALSE)

write.table(polyA_hpf4.reshape, file="Nano3pseq_4hpf_pA_merged_processed.tsv", sep="\t", quote=FALSE, row.names=FALSE)


ribodep_hpf2_merged.reshape <- read.delim("Nano3pseq_2hpf_merged_processed.tsv")
ribodep_hpf4_merged.reshape <- read.delim("Nano3pseq_4hpf_merged_processed.tsv")
ribodep_hpf6_merged.reshape <- read.delim("Nano3pseq_6hpf_merged_processed.tsv")

polyA_hpf4.reshape <- read.delim("Nano3pseq_4hpf_pA_merged_processed.tsv")


ribodep_hpf2_merged.mrna <- subset(ribodep_hpf2_merged.reshape, Gene_Type=="protein_coding")
ribodep_hpf4_merged.mrna <- subset(ribodep_hpf4_merged.reshape, Gene_Type=="protein_coding")
ribodep_hpf6_merged.mrna <- subset(ribodep_hpf6_merged.reshape, Gene_Type=="protein_coding")

polyA_hpf4.mrna <- subset(polyA_hpf4.reshape, Gene_Type=="protein_coding")


#######################################################
#######################################################
#######################################################



# 1. OVERALL COMPARISON OF THE TAIL LENGTH PER-READ

density_single_transcript <- function(drna, nano3p, label) {
  both <- rbind(drna, nano3p)
  pdf(file=paste(label, "dRNA_vs_Nano3Pseq_Overall_Tail_Comparison_Single_Transcript.pdf",sep="_"),height=3,width=5,onefile=FALSE)
  print(ggplot(both, aes(x=tail_length, color=Sample)) +
    geom_density()+
      theme_bw()+
      xlim(-10, 350))
  dev.off()
}

density_single_transcript(drna_hpf2.mrna, ribodep_hpf2_merged.mrna, "2hpf")
density_single_transcript(drna_hpf4.mrna, ribodep_hpf4_merged.mrna, "4hpf")
density_single_transcript(drna_hpf6.mrna, ribodep_hpf6_merged.mrna, "6hpf")

density_single_transcript(drna_hpf4.mrna, polyA_hpf4.mrna, "4hpf_pA")

# 2. OVERALL COMPARISON OF THE TAIL LENGTH MEDIAN PER GENE


density_median_transcript <- function(drna, nano3p, label) {
  both <- rbind(drna, nano3p)
  both2 <-  both[!duplicated(both[c("Gene_Name", "Sample")]),]
  pdf(file=paste(label, "dRNA_vs_Nano3Pseq_Overall_Tail_Comparison_Medianpergene.pdf",sep="_"),height=3,width=5,onefile=FALSE)
  print(ggplot(both2, aes(x=tail_length, color=Sample)) +
    geom_density()+
      theme_bw()+
      xlim(-10, 350))
  dev.off()
}

density_median_transcript(drna_hpf2.mrna, ribodep_hpf2_merged.mrna, "2hpf")
density_median_transcript(drna_hpf4.mrna, ribodep_hpf4_merged.mrna, "4hpf")
density_median_transcript(drna_hpf6.mrna, ribodep_hpf6_merged.mrna, "6hpf")

density_median_transcript(drna_hpf4.mrna, polyA_hpf4.mrna, "4hpf_pA")



# 3. SCATTER PLOT COMPARISON


merge_two <- function(drna, nano3p) {
  drna_unique <-   drna[!duplicated(drna[c("Gene_Name", "Sample")]),]
  nano3p_unique <-   nano3p[!duplicated(nano3p[c("Gene_Name", "Sample")]),]
  data_merged <- merge(drna_unique,nano3p_unique, by.x="Gene_Name", by.y="Gene_Name" )
  data_merged_mRNA <-subset(data_merged, Gene_Type.x =="protein_coding")
  data_merged_mRNA_min30 <- subset(data_merged_mRNA, Gene_Count.x > 30 & Gene_Count.y >30)
  return(data_merged_mRNA_min30)
}

data_merged_min30_2hpf <- merge_two(drna_hpf2.mrna, ribodep_hpf2_merged.mrna)
data_merged_min30_4hpf <- merge_two(drna_hpf4.mrna, ribodep_hpf4_merged.mrna)
data_merged_min30_6hpf <- merge_two(drna_hpf6.mrna, ribodep_hpf6_merged.mrna)

data_merged_min30_4hpf_pA <- merge_two(drna_hpf4.mrna, polyA_hpf4.mrna)



plot_denscols_with_corr_pearson<-function(pdfname,my_x,my_y,xlab,ylab) {
  pdf(file=paste(pdfname, "_pearson.pdf",sep=""), height=6, width=6)
  dcols<-densCols(my_x,my_y, colramp=colorRampPalette(blues9[-(1:3)]))
  plot(my_x,my_y,col=dcols,cex=1, cex.lab=1,cex.main=3,lwd=5,pch=20,xlab=xlab,ylab=ylab, xlim=c(0,160), ylim=c(0,160))
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


plot_denscols_with_corr_pearson("dRNA_vs_Nano3P_Median_Tail_2HPF_MinCov30", data_merged_min30_2hpf$Median_Length.x, data_merged_min30_2hpf$Median_Length.y, "dRNA", "Nano3P-seq" )

plot_denscols_with_corr_pearson("dRNA_vs_Nano3P_Median_Tail_4HPF_MinCov30", data_merged_min30_4hpf$Median_Length.x, data_merged_min30_4hpf$Median_Length.y, "dRNA", "Nano3P-seq" )

plot_denscols_with_corr_pearson("dRNA_vs_Nano3P_Median_Tail_6HPF_MinCov30", data_merged_min30_6hpf$Median_Length.x, data_merged_min30_6hpf$Median_Length.y, "dRNA", "Nano3P-seq" )


plot_denscols_with_corr_pearson("dRNA_vs_Nano3P_pA_Median_Tail_4HPF_MinCov30", data_merged_min30_4hpf_pA$Median_Length.x, data_merged_min30_4hpf_pA$Median_Length.y, "dRNA", "Nano3P-seq" )





# 4. ADD TAIL INFORMATION OF 3 GROUPS
#Order the groups
  maternal <- read.delim("allTranscripts_riboZero_rnaSeq.maternal.txt", header=FALSE)
  maternal$Group <- "Maternal"
  zyfir <- read.delim("allTranscripts_riboZero_rnaSeq.zyfir.txt", header=FALSE)
  zyfir$Group <- "Zygotic"
  mir430 <- read.delim("allTranscripts_riboZero_rnaSeq.mir430.txt", header=FALSE)
  mir430$Group <- "mir430"
  id_convert <- read.delim("Transcript_ID_To_Gene_ID.txt")
  groups <- rbind(maternal, zyfir, mir430)
  groups_id <- merge(groups,id_convert, by.x="V1", by.y="Transcript.stable.ID")
  groups_id$V1 <- NULL



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
  hpf246_2 <- merge(hpf246, groups_id,all.x=TRUE,  by.x="Gene_Name", "Gene.stable.ID")
  hpf246_2$Group[is.na(hpf246_2$Group)] <- "Rest"
  hpf246_final <- melt(hpf246_2)
  hpf246_final <- hpf246_final[complete.cases(hpf246_final), ]
return(hpf246_final)
}



dRNA_threegroups <- ThreeGroups_Tails(drna_hpf2.mrna, drna_hpf4.mrna,drna_hpf6.mrna, 20)
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



boxplots(dRNA_threegroups, "dRNA")
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
  hpf246_2 <- merge(hpf246, groups_id,all.x=TRUE,  by.x="Gene_Name", "Gene.stable.ID")
  hpf246_2$Group[is.na(hpf246_2$Group)] <- "Rest"
return(hpf246_2)
}


dRNA_threegroups <- ThreeGroups_Tails_Perread(drna_hpf2.mrna, drna_hpf4.mrna,drna_hpf6.mrna)
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
      stat_compare_means(comparisons = my_comparisons, label.y = c(150, 160, 170))+
       facet_wrap(~Group,nrow=1)+
      stat_n_text() + 
      #scale_fill_manual(values=colors)+
      #coord_cartesian(ylim = c(0,175))+
      theme_bw())
  dev.off()
}



boxplots_perread(dRNA_threegroups, "dRNA")
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
  #Make the NAs 0
  #hpf246_characters <- hpf246[,c(1)]
  #hpf246_numbers <- hpf246[,-c(1)]
  hpf246_final <- hpf246
  maternal <- read.delim("264_top_maternal_decay.txt")
  maternal$Group1 <- "Maternal"
  hpf246_final2 <- merge(hpf246_final, maternal,all.x=TRUE, by.x="Gene_Name", by.y="Ensembl_Gene_ID")
  hpf246_final2$Group1[is.na(hpf246_final2$Group1)] <- "Rest"
  hpf246_final2$Ensembl_Transcript_ID <- NULL
  #Import the Three Groups
  maternal <- read.delim("allTranscripts_riboZero_rnaSeq.maternal.txt", header=FALSE)
  maternal$Group2 <- "Maternal"
  zyfir <- read.delim("allTranscripts_riboZero_rnaSeq.zyfir.txt", header=FALSE)
  zyfir$Group2 <- "Zyfir"
  mir430 <- read.delim("allTranscripts_riboZero_rnaSeq.mir430.txt", header=FALSE)
  mir430$Group2 <- "mir430"
  id_convert <- read.delim("Transcript_ID_To_Gene_ID.txt")
  groups <- rbind(maternal, zyfir, mir430)
  groups_id <- merge(groups,id_convert, by.x="V1", by.y="Transcript.stable.ID")
  groups_id$V1 <- NULL
  #Merge the three groups with the data
  hpf246_final3 <- merge(hpf246_final2, groups_id, all.x=TRUE, all.y=TRUE, by.x="Gene_Name", by.y="Gene.stable.ID")
  hpf246_final3$Group2[is.na(hpf246_final3$Group2)] <- "Rest"
  #Rearrange the columns
  hpf246_final4 <- hpf246_final3[,c("Gene_Name","Group1","Group2", "Gene_Count_Norm.2hpf","Gene_Count_Norm.4hpf","Gene_Count_Norm.6hpf" )]
  #Melt the data 
  hpf246_melted <- melt(hpf246_final4)
  return(hpf246_melted)
}

dRNA_threegroups_counts <- three_groups_counts(drna_hpf2.mrna, drna_hpf4.mrna,drna_hpf6.mrna)
Nano3P_threegroups_counts <- three_groups_counts(ribodep_hpf2_merged.mrna, ribodep_hpf4_merged.mrna,ribodep_hpf6_merged.mrna)



#median_counts <- aggregate(.~Group2+variable, ribodep_all_unique_246hpf_taillengths_melted_mrna_notrelativeto2hpf[,c("Group2", "variable", "value")], median)

boxplots_nano3Pseq <- function(data, label) {
my_comparisons <- list( c("Gene_Count_Norm.2hpf", "Gene_Count_Norm.4hpf"), c("Gene_Count_Norm.2hpf", "Gene_Count_Norm.6hpf"), c("Gene_Count_Norm.4hpf", "Gene_Count_Norm.6hpf") )
pdf(file=paste(label, "Zebrafish_EmbryosCounts_Boxplot_Maternal_vs_Rest_3Groups.pdf",sep="_"),height=5,width=10,onefile=FALSE)
    print(ggplot(data, aes(x = variable, y = log(value+1))) + 
      geom_boxplot(aes(fill = variable),position = position_dodge(0.9), outlier.shape=NA) +
      ylab("log2(Count)")+
      stat_compare_means(comparisons = my_comparisons, label.y = c(0.8, 0.9, 1))+
       facet_wrap(~Group2,nrow=1)+
      #stat_n_text() + 
      #scale_fill_manual(values=colors)+
      coord_cartesian(ylim = c(0,1.5))+
      theme_bw())
  dev.off()
}

boxplots_nano3Pseq(Nano3P_threegroups_counts, "Nano3P_Seq")



boxplots_dRNAseq <- function(data, label) {
my_comparisons <- list( c("Gene_Count_Norm.2hpf", "Gene_Count_Norm.4hpf"), c("Gene_Count_Norm.2hpf", "Gene_Count_Norm.6hpf"), c("Gene_Count_Norm.4hpf", "Gene_Count_Norm.6hpf") )
pdf(file=paste(label, "Zebrafish_EmbryosCounts_Boxplot_Maternal_vs_Rest_3Groups.pdf",sep="_"),height=5,width=10,onefile=FALSE)
    print(ggplot(data, aes(x = variable, y = log(value+1))) + 
      geom_boxplot(aes(fill = variable),position = position_dodge(0.9), outlier.shape=NA) +
      ylab("log2(Count)")+
      stat_compare_means(comparisons = my_comparisons, label.y = c(1.5, 1.8, 2))+
       facet_wrap(~Group2,nrow=1)+
      #stat_n_text() + 
      #scale_fill_manual(values=colors)+
      coord_cartesian(ylim = c(0,2.5))+
      theme_bw())
  dev.off()
}


boxplots_dRNAseq(dRNA_threegroups_counts, "dRNA")











# 6. USE GR GROUPS
gr_groups <- read.delim("gr_list.txt")


GR_Groups_Tails <- function(hpf2, hpf4, hpf6, cov_th) {
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
  hpf246_2$Group[is.na(hpf246_2$Group)] <- "rest"
  hpf246_final <- melt(hpf246_2)
  hpf246_final <- hpf246_final[complete.cases(hpf246_final), ]
return(hpf246_final)
}



dRNA_gr_groups <- GR_Groups_Tails(drna_hpf2.mrna, drna_hpf4.mrna,drna_hpf6.mrna, 20)
Nano3P_gr_groups <- GR_Groups_Tails(ribodep_hpf2_merged.mrna, ribodep_hpf4_merged.mrna,ribodep_hpf6_merged.mrna,20)





boxplots_gr <- function(data, label) {
my_comparisons <- list( c("Median_Length.2hpf", "Median_Length.4hpf"), c("Median_Length.2hpf", "Median_Length.6hpf"), c("Median_Length.4hpf", "Median_Length.6hpf") )
pdf(file=paste(label, "Zebrafish_Embryos_TailLength_Boxplot_GR_3Groups_log.pdf",sep="_"),height=5,width=10,onefile=FALSE)
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

pdf(file=paste(label, "Zebrafish_Embryos_TailLength_Boxplot_GR_3Groups.pdf",sep="_"),height=5,width=10,onefile=FALSE)
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



boxplots_gr(dRNA_gr_groups, "dRNA")
boxplots_gr(Nano3P_gr_groups, "Nano3PSeq")







