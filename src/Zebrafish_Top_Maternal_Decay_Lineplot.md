########################################################
######## ZEBRAFISH Top Maternal Decay Line Plot ########
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






#### INVESTIGATE THE MATERNAL RNAs GENE COUNT
ribodep_all_unique <-ribodep_all_merged[!duplicated(ribodep_all_merged[c("Gene_Name", "Sample")]),]


ribodep_all_unique_2hpf <- subset(ribodep_all_unique, Sample=="Ribodep_2hpf_merged")
ribodep_all_unique_4hpf <- subset(ribodep_all_unique, Sample=="Ribodep_4hpf_merged")
ribodep_all_unique_6hpf <- subset(ribodep_all_unique, Sample=="Ribodep_6hpf_merged")


ribodep_all_unique_2hpf <- ribodep_all_unique_2hpf[,c("Gene_Type", "Gene_Name", "Median_Length", "Gene_Count_Norm")]
ribodep_all_unique_4hpf <- ribodep_all_unique_4hpf[,c("Gene_Type", "Gene_Name", "Median_Length", "Gene_Count_Norm")]
ribodep_all_unique_6hpf <- ribodep_all_unique_6hpf[,c("Gene_Type", "Gene_Name", "Median_Length", "Gene_Count_Norm")]

columns <-c("Gene_Type", "Gene_Name")

ribodep_all_unique_2hpf_4hpf <- merge(ribodep_all_unique_2hpf, ribodep_all_unique_4hpf, all.x=TRUE, all.y=TRUE, by.x=columns, by.y=columns, suffixes = c(".2hpf", ".4hpf"))

ribodep_all_unique_246hpf <- merge(ribodep_all_unique_2hpf_4hpf, ribodep_all_unique_6hpf, all.x=TRUE, all.y=TRUE, by.x=columns, by.y=columns)

colnames(ribodep_all_unique_246hpf) <- c(colnames(ribodep_all_unique_2hpf_4hpf),"Median_Length.6hpf", "Gene_Count_Norm.6hpf")

ribodep_all_unique_246hpf_characters <- ribodep_all_unique_246hpf[,c(1,2)]
ribodep_all_unique_246hpf_numbers <- ribodep_all_unique_246hpf[,-c(1,2)]

ribodep_all_unique_246hpf_numbers[is.na(ribodep_all_unique_246hpf_numbers)] <- 0

ribodep_all_unique_246hpf_final <- cbind(ribodep_all_unique_246hpf_characters,ribodep_all_unique_246hpf_numbers )


maternal <- read.delim("264_top_maternal_decay.txt")
maternal$Group1 <- "Maternal"

ribodep_all_unique_246hpf_final2 <- merge(ribodep_all_unique_246hpf_final, maternal,all.x=TRUE, by.x="Gene_Name", by.y="Ensembl_Gene_ID")
ribodep_all_unique_246hpf_final2$Group1[is.na(ribodep_all_unique_246hpf_final2$Group1)] <- "Rest"
ribodep_all_unique_246hpf_final2$Ensembl_Transcript_ID <- NULL


ribodep_all_unique_246hpf_genecounts <- ribodep_all_unique_246hpf_final2[,c("Gene_Name","Gene_Type","Group1", "Gene_Count_Norm.2hpf","Gene_Count_Norm.4hpf","Gene_Count_Norm.6hpf" )]


#Normalize by 2hpf
ribodep_all_unique_246hpf_genecounts$Gene_Count_Norm2.2hpf <- (ribodep_all_unique_246hpf_genecounts$Gene_Count_Norm.2hpf+1)/(ribodep_all_unique_246hpf_genecounts$Gene_Count_Norm.2hpf+1)
ribodep_all_unique_246hpf_genecounts$Gene_Count_Norm2.4hpf <- (ribodep_all_unique_246hpf_genecounts$Gene_Count_Norm.4hpf+1)/(ribodep_all_unique_246hpf_genecounts$Gene_Count_Norm.2hpf+1)
ribodep_all_unique_246hpf_genecounts$Gene_Count_Norm2.6hpf <- (ribodep_all_unique_246hpf_genecounts$Gene_Count_Norm.6hpf+1)/(ribodep_all_unique_246hpf_genecounts$Gene_Count_Norm.2hpf+1)

ribodep_all_unique_246hpf_genecounts_relativeto2hpf <- ribodep_all_unique_246hpf_genecounts[,c("Gene_Name","Gene_Type", "Group1", "Gene_Count_Norm2.2hpf" ,"Gene_Count_Norm2.4hpf", "Gene_Count_Norm2.6hpf")]


ribodep_all_unique_246hpf_genecounts_melted_relativeto2hpf <- melt(ribodep_all_unique_246hpf_genecounts_relativeto2hpf)
ribodep_all_unique_246hpf_genecounts_melted_mrna_relativeto2hpf  <- subset(ribodep_all_unique_246hpf_genecounts_melted_relativeto2hpf, Gene_Type=="protein_coding")



pdf(file="Zebrafish_Embryos_GeneCount_Line_Plot_Maternal_vs_Rest_Relative_to_2hpf.pdf",height=4,width=5,onefile=FALSE)
print(ggplot(ribodep_all_unique_246hpf_genecounts_melted_mrna_relativeto2hpf, aes(x=variable, y=log(value), group=Gene_Name)) +
  geom_line(data=subset(ribodep_all_unique_246hpf_genecounts_melted_mrna_relativeto2hpf, Group1=="Rest"), alpha=1/5, color="grey")+
  geom_line(data=subset(ribodep_all_unique_246hpf_genecounts_melted_mrna_relativeto2hpf, Group1=="Maternal"), alpha=1/2, color="red")+
  theme_bw()+
  geom_point(data=subset(ribodep_all_unique_246hpf_genecounts_melted_mrna_relativeto2hpf, Group1=="Rest"), color="grey")+
    geom_point(data=subset(ribodep_all_unique_246hpf_genecounts_melted_mrna_relativeto2hpf, Group1=="Maternal"), color="red"))

dev.off()



