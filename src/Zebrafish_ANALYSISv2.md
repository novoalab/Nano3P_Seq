########################################################
######## ANALYSIS OF THE ZEBRAFISH RUNS REPS MERGED ####
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



ribodep_all_unique <-ribodep_all_merged[!duplicated(ribodep_all_merged[c("Gene_Name", "Sample")]),]
write.table(ribodep_all.unique, file="Zebrafish_Ribodepletion_RepsMerged_Reshaped.tsv", sep="\t", quote=FALSE,row.names=FALSE)






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

ribodep_hpf2.tail_class_count <- count_polyA_population(ribodep_hpf2_merged.reshape)
ribodep_hpf4.tail_class_count <- count_polyA_population(ribodep_hpf4_merged.reshape)
ribodep_hpf6.tail_class_count <- count_polyA_population(ribodep_hpf6_merged.reshape)







#### INVESTIGATE THE MATERNAL RNAs GENE COUNTS

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




#Gene profile
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


ribodep_all_unique_246hpf_final3 <- merge(ribodep_all_unique_246hpf_final2, groups_id, all.x=TRUE, all.y=TRUE, by.x="Gene_Name", by.y="Gene.stable.ID")
ribodep_all_unique_246hpf_final3$Group2[is.na(ribodep_all_unique_246hpf_final3$Group2)] <- "Rest"



##### ANALYSE NONCODING RNA POPULATION

ribodep_all_unique_246hpf_ncRNA <- subset(ribodep_all_unique_246hpf_final, Gene_Type=="lincRNA" |Gene_Type=="snoRNA" | Gene_Type=="snRNA"| Gene_Type=="scaRNA" )
write.table(ribodep_all_unique_246hpf_ncRNA, file="Zebrafish_ncRNA_Expression_Tail.tsv", sep="\t", row.names=FALSE, quote=FALSE)


ribodep_all_unique_246hpf_heatmap <- ribodep_all_unique_246hpf_ncRNA[,c("Gene_Type", "Gene_Name", "Gene_Count_Norm.2hpf", "Gene_Count_Norm.4hpf", "Gene_Count_Norm.6hpf")]



library(ComplexHeatmap)
library(circlize)


data <- ribodep_all_unique_246hpf_heatmap

rownames(data)<- data[,2] #assign gene names as rownames  
data2<- data[,-c(1:2)] #remove the first three columns for the heatmap
data3 <- t(scale(t(data2)))#Normalize by row (by gene)

pdf("heatmap.ncRNA_genecounts.zscaled.pdf",height=12,width=8)
Heatmap(data3, name = "z-scale Counts", 
  #col = colorRamp2(c(-3,0,4), c("cadetblue3","floralwhite", "maroon4"),space = "RGB"), 
    #cluster_rows = FALSE, 
    col = colorRamp2(c(-3,-1.5,0,1.5,3), c("#2c7bb6","#abd9e9","floralwhite","#fdae61", "#d7191c"),space = "RGB"),
    cluster_columns = FALSE,
    column_title = "Time-point", 
    column_title_gp = gpar(fontsize = 10, fontface = "bold"),
    column_names_gp = gpar(fontsize = 7, fontface = "bold"),
    row_title = "ncRNAs", row_title_rot = 90,
    row_title_gp = gpar(fontsize = 8, fontface = "bold"),
    cluster_rows = TRUE,
    show_row_names = FALSE,
    row_names_gp = gpar(fontsize = 5), #row names size
    column_order = 1:dim(data3)[2],#Keep the column order, make clustering FALSE for this
    row_dend_side = "right", #Dendogram on the right side
    #row_order = 1:dim(data4)[1], #Keep the row order, make clustering FALSE for this
    #show_column_dend = TRUE, #
    #column_dend_side = "top",
    column_names_side = "bottom",
    split = data$Gene_Type, #Splitting by Class
    gap = unit(1, "mm"), #Gap
    )
dev.off()




rownames(data)<- data[,2] #assign gene names as rownames  
data2<- data[,-c(1:2)] #remove the first three columns for the heatmap


data2_log <- log((data2*100)+1)



pdf("heatmap.ncRNA_genecounts.logcounts.pdf",height=12,width=8)
Heatmap(data2_log, name = "log(Counts)", 
  #col = colorRamp2(c(-3,0,4), c("cadetblue3","floralwhite", "maroon4"),space = "RGB"), 
    #cluster_rows = FALSE, 
    col = colorRamp2(c(0,1,2), c("gray","pink","red"),space = "RGB"),
    cluster_columns = FALSE,
    column_title = "Time-point", 
    column_title_gp = gpar(fontsize = 10, fontface = "bold"),
    column_names_gp = gpar(fontsize = 7, fontface = "bold"),
    #row_title = "ncRNAs", row_title_rot = 90,
    row_title_gp = gpar(fontsize = 8, fontface = "bold"),
    cluster_rows = TRUE,
    show_row_names = FALSE,
    row_names_gp = gpar(fontsize = 5), #row names size
    column_order = 1:dim(data3)[2],#Keep the column order, make clustering FALSE for this
    row_dend_side = "right", #Dendogram on the right side
    #row_order = 1:dim(data4)[1], #Keep the row order, make clustering FALSE for this
    #show_column_dend = TRUE, #
    #column_dend_side = "top",
    column_names_side = "bottom",
    row_split = data$Gene_Type, #Splitting by Class
    gap = unit(1, "mm"), #Gap
    )
dev.off()




















### Count pA on mtrRNAs

hpf2_mtrRNA <- subset(ribodep_hpf2.tail_class_count, Gene_Type=="Mt_rRNA")
hpf4_mtrRNA <- subset(ribodep_hpf4.tail_class_count, Gene_Type=="Mt_rRNA")
hpf6_mtrRNA <- subset(ribodep_hpf6.tail_class_count, Gene_Type=="Mt_rRNA")

hpf4_mtrRNA_pA <- subset(polyA_hpf4.tail_class_count, Gene_Type=="Mt_rRNA")



mtrRNA_zebrafish_ribodep <- rbind(hpf2_mtrRNA, hpf4_mtrRNA,hpf6_mtrRNA,hpf4_mtrRNA_pA)
mtrRNA_zebrafish_ribodep_uniq <-mtrRNA_zebrafish_ribodep[!duplicated(mtrRNA_zebrafish_ribodep[c("Gene_Name", "Sample")]),]












