########################################################
######## ZEBRAFISH Overall Tail and Gene Count #########
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





dotplot <- function(data, label) {
    data2 <- subset(data , Gene_Count > 30)
    data2 <- subset(data2, Gene_Type!="synthetic" & Gene_Type !="ensembl_havana" & Gene_Type!="havana"& Gene_Type!="Mt_tRNA")
    data2 <- data2[order(-data2$Gene_Type_Count_Norm),]
     data2$Gene_Type <- factor(data2$Gene_Type, levels = unique(data2$Gene_Type))
    pdf(file=paste(label, "Tail_Length_Distribution_All_Genes_Dotplot_Min30.pdf",sep="_"),height=5,width=12,onefile=FALSE)
      print(ggplot(data2, aes(x=Gene_Type, y=tail_length, color=Gene_Type)) + 
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
  dotplot(ribodep_hpf2_merged.reshape,"2HPF")
  dotplot(ribodep_hpf4_merged.reshape,"4HPF")
  dotplot(ribodep_hpf6_merged.reshape,"6HPF")




boxplot <- function(data, label) {
    data2 <- subset(data , Gene_Count > 30)
    data2 <- subset(data2, Gene_Type!="synthetic" & Gene_Type !="ensembl_havana" & Gene_Type!="havana"& Gene_Type!="Mt_tRNA"& Gene_Type!="snoRNA"& Gene_Type!="polymorphic_pseudogene"& Gene_Type!="unprocessed_pseudogene"& Gene_Type!="transcribed_unprocessed_pseudogene"& Gene_Type!="TEC"& Gene_Type!="scaRNA"& Gene_Type!="ribozyme")
    data2 <- data2[order(-data2$Gene_Type_Count_Norm),]
    data2$Gene_Type <- factor(data2$Gene_Type, levels = unique(data2$Gene_Type))
      pdf(file=paste(label,"Tails_Boxplot.pdf", sep="_"),height=4,width=20, onefile=FALSE)
            print(ggplot(data2, aes(x = Gene_Type, y = tail_length, fill=Sample )) + 
              geom_boxplot(aes(fill = Sample),position = position_dodge(0.9), outlier.shape = NA) +
              theme_bw()+
              coord_cartesian(ylim = c(0, 250))
              #scale_fill_manual(values = c("#91c788", "#ffaaa7"))
              )
         dev.off()
}

boxplot(ribodep_all_merged, "Ribodepletion_All_Timepoints")











##Overall Tail comparison (Single transcript)
tail_comparison_overall <- function(data, label) {
  data2 <- subset(data, Gene_Type=="lincRNA"|Gene_Type=="Mt_rRNA"|Gene_Type=="protein_coding"|Gene_Type=="rRNA" )
  data_means <- aggregate(.~Sample+Gene_Type, data2[,c("Sample", "Median_Length", "Gene_Type")], median)
  print(data_means)
  pdf(file=paste(label, "Overall_Tail_Comparison_Single_Transcript.pdf",sep="_"),height=6,width=10,onefile=FALSE)
  print(ggplot(data2, aes(x=tail_length, color=Sample)) +
    geom_density()+
    theme_bw()+
    #xlim(-10, 350)+
    facet_wrap(~Gene_Type, scales="free"))
  dev.off()
}


tail_comparison_overall(ribodep_all_merged, "Ribodepletion_All_Timepoints_Merged")






#### Gene Type Count Comparison ###
simple_barplot_grouped <- function(data, label){
  #Remove duplicates
  data2 <-  data[!duplicated(data[c("Gene_Type", "Sample")]),]
  data2 <- subset(data2, Gene_Type != "synthetic" & Gene_Type != "ensembl_havana"  & Gene_Type != "ensembl" & Gene_Type != "havana" & Gene_Type != "polymorphic_pseudogene" & Gene_Type != "TEC" & Gene_Type != "unprocessed_pseudogene" & Gene_Type != "transcribed_unprocessed_pseudogene" & Gene_Type != "snRNA"& Gene_Type != "TR_J_gene"& Gene_Type != "sense_intronic"& Gene_Type != "Mt_tRNA"& Gene_Type != "processed_transcript")
  category_sum <- aggregate(.~Gene_Type+Sample, data2[,c("Gene_Type","Sample",  "Gene_Type_Count_Norm")], sum)
  category_sum <- category_sum[order(-category_sum$Gene_Type_Count_Norm),]
  category_sum$Gene_Type <- factor(category_sum$Gene_Type, levels = unique(category_sum$Gene_Type))
  pdf(file=paste(label, "Gene_Type_Normalized_Count.pdf",sep="_"),height=6,width=20,onefile=FALSE)
  print(ggplot(category_sum, aes(fill=Sample,  y=log(Gene_Type_Count_Norm+1), x=Gene_Type)) + 
    geom_bar(position="dodge", stat="identity")+
    #scale_y_continuous(trans = 'log2')+
    theme_bw())
    dev.off()
}



simple_barplot_grouped(ribodep_all_merged, "Ribodepletion")






## Tail comparison median per gene
tail_comparison_median_per_gene_protein<- function(data, label) {
  data2 <- subset(data, Gene_Type=="protein_coding")
  data2 <-  data2[!duplicated(data2[c("Gene_Name", "Sample")]),]
  data_means <- aggregate(.~Sample, data2[,c("Sample", "Median_Length")], median)
  print(data_means)
  pdf(file=paste(label, "Median_Tail_Per_Gene_Comparison_mRNA.pdf",sep="_"),height=6,width=10,onefile=FALSE)
  print(ggplot(data2, aes(x=Median_Length, color=Sample)) +
    geom_density()+
      theme_bw()+
      coord_cartesian(xlim=c(-10, 200)))
  dev.off()
}


tail_comparison_median_per_gene_protein(ribodep_all_merged, "Ribodepletion_All_Timepoints_Merged")





  ribodep_all_unique <-ribodep_all_merged[!duplicated(ribodep_all_merged[c("Gene_Name", "Sample")]),]
  ribodep_all_unique_min30 <-  subset(ribodep_all_unique, Gene_Count > 30)
  ribodep_all_unique_min30_mRNA <- subset(ribodep_all_unique_min30, Gene_Type =="protein_coding")
  

    pdf(file= "Tails_Ribodepleted_Timepoints_ProteinCoding_Tails_Min30_merged.pdf",height=10,width=10,onefile=FALSE)
      print(ggplot(ribodep_all_unique_min30_mRNA, aes(x=Sample, y=Median_Length)) + 
        geom_quasirandom(varwidth = TRUE, aes(color=Sample))+
        geom_boxplot(aes(alpha=0), outlier.shape=NA)+
        #geom_text_repel(data=subset(ribodep_all_unique_min20_mRNA,Median_Length==0), aes(label=Gene.name))+
        stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                geom = "crossbar", width = 0.7, color="#c06c84")+
        theme_bw()+
        ggtitle("Zebrafish Embryo")+
        xlab("Time Points")+
                ylab("Tail length") +
        theme(axis.text=element_text(size=14),strip.text = element_text(size=13),
                axis.title=element_text(size=17,face="bold"),
                legend.title = element_text(size = 20),
                legend.text = element_text(color = "black", size=15)))
    dev.off()
  

