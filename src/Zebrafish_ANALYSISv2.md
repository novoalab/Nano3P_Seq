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



library(EnvStats)


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










### COMPARE WITH LITERATURE 
bartel_hpf2.data <- read.delim("GSE52809_Dre_mock_2hpf.txt")
bartel_hpf4.data <- read.delim("GSE52809_Dre_mock_4hpf.txt")
bartel_hpf6.data <- read.delim("GSE52809_Dre_mock_6hpf.txt")


bartel_merging <- function(data_bartel,data_us, threshold) {
   id_convert <- read.delim("Transcript_ID_To_Gene_ID.txt")
   data_bartel2 <- data_bartel[,c("Transcript.ID", "Mean.TL", "Median.TL")] 
   data_bartel3 <- merge(data_bartel2, id_convert, by.x="Transcript.ID", by.y="Transcript.stable.ID")
   #Take unique gene names for our data
   data_us2 <-  data_us[!duplicated(data_us[c("Gene_Name", "Sample")]),]
   data_us3 <- subset(data_us2, Gene_Type =="protein_coding")
   data_us4 <- data_us3[,c("Gene_Name","Median_Length","Mean_Length","Gene_Count","Gene_Count_Norm","Sample")]
   hpf2_comparison <- merge(data_us4, data_bartel3, by.x="Gene_Name", by.y="Gene.stable.ID")
   hpf2_comparison_filtered <- subset(hpf2_comparison, Gene_Count > threshold)
   return(hpf2_comparison_filtered)
  }

bartel_hpf2.processed <- bartel_merging(bartel_hpf2.data, ribodep_hpf2_merged.reshape, 30)
bartel_hpf4.processed <- bartel_merging(bartel_hpf4.data, ribodep_hpf4_merged.reshape, 30)
bartel_hpf6.processed <- bartel_merging(bartel_hpf6.data, ribodep_hpf6_merged.reshape, 30)



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
  plot(my_x,my_y,col=dcols,cex=1, cex.lab=1,cex.main=3,lwd=5,pch=20,xlab=xlab,ylab=ylab)
  title(main=pdfname, col.main="black", font.main=4)
  #abline(v=0, lty=2)
  # Correlation
  test<-cor.test(my_x,my_y, method="pearson")
  print(test)
  cor222<-paste("Pearson's p =",round(as.numeric(test$estimate),3))
  #pval<-paste("Pval =",test$p.value)
  mtext(paste(cor222))
  #mtext(paste(cor222,pval,sep=" ; ")) #Print the subtitle with the dataset correlation
  dev.off()
}

plot_denscols_with_corr_pearson("PALSeq_vs_Nano3PSeq 2HPF", bartel_hpf2.processed$Median_Length , bartel_hpf2.processed$Median.TL, "NanoTailCapture (this work) Median Tail Length", "PAL-Seq (Subtelny et al 2014) Median Tail Length" )

plot_denscols_with_corr_pearson("PALSeq_vs_Nano3PSeq 4HPF", bartel_hpf4.processed$Median_Length , bartel_hpf4.processed$Median.TL, "NanoTailCapture (this work) Median Tail Length", "PAL-Seq (Subtelny et al 2014) Median Tail Length" )

plot_denscols_with_corr_pearson("PALSeq_vs_Nano3PSeq 6HPF", bartel_hpf6.processed$Median_Length , bartel_hpf6.processed$Median.TL, "NanoTailCapture (this work) Median Tail Length", "PAL-Seq (Subtelny et al 2014) Median Tail Length" )









### BOXPLOT OF THE MEDIAN TAIL PER GENE IN RIBODEPLETED RUN
  GENE_ID <- read.delim("Zebrafish_ID_Name_Conversion.txt")
  ribodep_all_unique <-ribodep_all_merged[!duplicated(ribodep_all_merged[c("Gene_Name", "Sample")]),]
  ribodep_all_unique_names <- merge(ribodep_all_unique, GENE_ID, by.x="Gene_Name", by.y="Gene.stable.ID")
  ribodep_all_unique_min20 <-  subset(ribodep_all_unique_names, Gene_Count > 20)
  ribodep_all_unique_min20_mRNA <- subset(ribodep_all_unique_min20, Gene_Type =="protein_coding")
  

  
  
    pdf(file= "Tails_Ribodepleted_Timepoints_ProteinCoding_Tails_Min20_merged.pdf",height=10,width=20,onefile=FALSE)
      print(ggplot(ribodep_all_unique_min20_mRNA, aes(x=Sample, y=Median_Length)) + 
        geom_quasirandom(varwidth = TRUE, aes(color=Sample))+
        geom_boxplot(aes(alpha=0), outlier.shape=NA)+
        geom_text_repel(data=subset(ribodep_all_unique_min20_mRNA,Median_Length==0), aes(label=Gene.name))+
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

polyA_hpf4.tail_class_count <- count_polyA_population(polyA_hpf4.reshape)




# LINE PLOT WITH THREE TIME POINTS
ribodep_hpf2.tail_class_count_unique <-ribodep_hpf2.tail_class_count[!duplicated(ribodep_hpf2.tail_class_count[c("Gene_Name", "Sample")]),]
ribodep_hpf4.tail_class_count_unique <-ribodep_hpf4.tail_class_count[!duplicated(ribodep_hpf4.tail_class_count[c("Gene_Name", "Sample")]),]
ribodep_hpf6.tail_class_count_unique <-ribodep_hpf6.tail_class_count[!duplicated(ribodep_hpf6.tail_class_count[c("Gene_Name", "Sample")]),]


columns <- c("Gene_Name", "Gene_Type", "Gene_Count", "Ratio_PolyA", "Ratio_NoPolyA")
merged_2_4_hours <-  merge(ribodep_hpf2.tail_class_count_unique[,columns], ribodep_hpf4.tail_class_count_unique[,columns], by.x=c("Gene_Name", "Gene_Type"),  by.y=c("Gene_Name", "Gene_Type"))

merged_2_4_6_hours <-  merge(merged_2_4_hours, ribodep_hpf6.tail_class_count_unique[,columns], by.x=c("Gene_Name", "Gene_Type"),  by.y=c("Gene_Name", "Gene_Type"))


colnames(merged_2_4_6_hours) <-  c("Gene_Name", "Gene_Type", "Gene_Count.2hpf","Ratio_PolyA.2hpf" , "Ratio_NoPolyA.2hpf","Gene_Count.4hpf" ,"Ratio_PolyA.4hpf" , "Ratio_NoPolyA.4hpf","Gene_Count.6hpf","Ratio_PolyA.6hpf" ,"Ratio_NoPolyA.6hpf" ) 

merged_2_4_6_hours_mRNA <- subset(merged_2_4_6_hours, Gene_Type=="protein_coding")

merged_2_4_6_hours_mRNA_min50 <- subset(merged_2_4_6_hours_mRNA, Gene_Count.2hpf> 10 & Gene_Count.4hpf> 10 & Gene_Count.6hpf> 10)

merged_2_4_6_hours_mRNA_min50_2 <- merged_2_4_6_hours_mRNA_min50[,c("Gene_Name", "Gene_Type","Ratio_PolyA.2hpf","Ratio_PolyA.4hpf" , "Ratio_PolyA.6hpf")]
merged_2_4_6_hours_mRNA_min50_2$Diff_2_4 <- abs(merged_2_4_6_hours_mRNA_min50_2$Ratio_PolyA.2hpf- merged_2_4_6_hours_mRNA_min50_2$Ratio_PolyA.4hpf)
merged_2_4_6_hours_mRNA_min50_2$Diff_4_6 <- abs(merged_2_4_6_hours_mRNA_min50_2$Ratio_PolyA.4hpf- merged_2_4_6_hours_mRNA_min50_2$Ratio_PolyA.6hpf)

merged_2_4_6_hours_mRNA_min50_changing <- subset(merged_2_4_6_hours_mRNA_min50_2,Diff_2_4>0.2 | Diff_4_6 > 0.2 )
merged_2_4_6_hours_mRNA_min50_changing$Type <- "Changing"

merged_2_4_6_hours_mRNA_min50_notchanging <- subset(merged_2_4_6_hours_mRNA_min50_2,Diff_2_4 < 0.2 & Diff_4_6 < 0.2)
merged_2_4_6_hours_mRNA_min50_notchanging$Type <- "NotChanging"


merged_2_4_6_hours_mRNA_min50_type <- rbind(merged_2_4_6_hours_mRNA_min50_changing,merged_2_4_6_hours_mRNA_min50_notchanging)

merged_2_4_6_hours_mRNA_min50_type2 <- merged_2_4_6_hours_mRNA_min50_type[,c("Gene_Name", "Gene_Type","Type","Ratio_PolyA.2hpf","Ratio_PolyA.4hpf" , "Ratio_PolyA.6hpf")]


library(reshape2)
merged_2_4_6_hours_mRNA_min50_melt <- melt(merged_2_4_6_hours_mRNA_min50_type2)


  pdf(file= "ZebrafishTimePoints_PolyA_Ratio_LinePlot_mRNAs_min10.pdf",height=5,width=10,onefile=FALSE)
    print(ggplot(data=merged_2_4_6_hours_mRNA_min50_melt, aes(x=variable, y=value, group=Gene_Name, colour=Type)) +
      geom_line(alpha=1/5)+
      scale_colour_manual(values=c("red","gray"))+
      theme_bw()+
      geom_text(data=subset(merged_2_4_6_hours_mRNA_min50_melt, value < 0.5), aes(label=Gene_Name))+
        geom_point())
    dev.off()





































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
zygotic <- read.delim("allTranscripts_riboZero_rnaSeq.zyfir.txt", header=FALSE)
zygotic$Group2 <- "Zygotic"

mir430 <- read.delim("allTranscripts_riboZero_rnaSeq.mir430.txt", header=FALSE)
mir430$Group2 <- "mir430"

id_convert <- read.delim("Transcript_ID_To_Gene_ID.txt")

groups <- rbind(maternal, zygotic, mir430)
groups_id <- merge(groups,id_convert, by.x="V1", by.y="Transcript.stable.ID")
groups_id$V1 <- NULL


ribodep_all_unique_246hpf_final3 <- merge(ribodep_all_unique_246hpf_final2, groups_id, all.x=TRUE, all.y=TRUE, by.x="Gene_Name", by.y="Gene.stable.ID")
ribodep_all_unique_246hpf_final3$Group2[is.na(ribodep_all_unique_246hpf_final3$Group2)] <- "Rest"





ribodep_all_unique_246hpf_genecounts <- ribodep_all_unique_246hpf_final3[,c("Gene_Name","Gene_Type","Group1","Group2", "Gene_Count_Norm.2hpf","Gene_Count_Norm.4hpf","Gene_Count_Norm.6hpf" )]


#Normalize by 2hpf
ribodep_all_unique_246hpf_genecounts$Gene_Count_Norm2.2hpf <- (ribodep_all_unique_246hpf_genecounts$Gene_Count_Norm.2hpf+1)/(ribodep_all_unique_246hpf_genecounts$Gene_Count_Norm.2hpf+1)
ribodep_all_unique_246hpf_genecounts$Gene_Count_Norm2.4hpf <- (ribodep_all_unique_246hpf_genecounts$Gene_Count_Norm.4hpf+1)/(ribodep_all_unique_246hpf_genecounts$Gene_Count_Norm.2hpf+1)
ribodep_all_unique_246hpf_genecounts$Gene_Count_Norm2.6hpf <- (ribodep_all_unique_246hpf_genecounts$Gene_Count_Norm.6hpf+1)/(ribodep_all_unique_246hpf_genecounts$Gene_Count_Norm.2hpf+1)

ribodep_all_unique_246hpf_genecounts_relativeto2hpf <- ribodep_all_unique_246hpf_genecounts[,c("Gene_Name","Gene_Type", "Group1","Group2", "Gene_Count_Norm2.2hpf" ,"Gene_Count_Norm2.4hpf", "Gene_Count_Norm2.6hpf")]


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





pdf(file="Zebrafish_Embryos_GeneCount_Line_Plot_Maternal_vs_Rest_Relative_to_2hpf_3Groups.pdf",height=10,width=10,onefile=FALSE)
print(ggplot(ribodep_all_unique_246hpf_genecounts_melted_mrna_relativeto2hpf, aes(x=variable, y=log(value), group=Gene_Name)) +
  geom_line(alpha=1, aes(color=Group2))+
  theme_bw()+
  facet_wrap(~Group2)+
  geom_point( aes(color=Group2)))
dev.off()










#### INVESTIGATE THE MATERNAL RNAs TAIL LENGTHS
GENE_ID <- read.delim("Zebrafish_ID_Name_Conversion.txt")
ribodep_all_unique <-ribodep_all_merged[!duplicated(ribodep_all_merged[c("Gene_Name", "Sample")]),]
ribodep_all_unique_2hpf <- subset(ribodep_all_unique, Sample=="Ribodep_2hpf_merged")
ribodep_all_unique_4hpf <- subset(ribodep_all_unique, Sample=="Ribodep_4hpf_merged")
ribodep_all_unique_6hpf <- subset(ribodep_all_unique, Sample=="Ribodep_6hpf_merged")

ribodep_all_unique_2hpf_mincov <- subset(ribodep_all_unique_2hpf, Gene_Count> 20)
ribodep_all_unique_4hpf_mincov <- subset(ribodep_all_unique_4hpf, Gene_Count> 20)
ribodep_all_unique_6hpf_mincov <- subset(ribodep_all_unique_6hpf, Gene_Count> 20)

ribodep_all_unique_2hpf_mincov <- ribodep_all_unique_2hpf_mincov[,c("Gene_Type", "Gene_Name", "Median_Length", "Gene_Count_Norm")]
ribodep_all_unique_4hpf_mincov <- ribodep_all_unique_4hpf_mincov[,c("Gene_Type", "Gene_Name", "Median_Length", "Gene_Count_Norm")]
ribodep_all_unique_6hpf_mincov <- ribodep_all_unique_6hpf_mincov[,c("Gene_Type", "Gene_Name", "Median_Length", "Gene_Count_Norm")]

columns <-c("Gene_Type", "Gene_Name")

ribodep_all_unique_2hpf_4hpf_mincov <- merge(ribodep_all_unique_2hpf_mincov, ribodep_all_unique_4hpf_mincov, all.x=TRUE, all.y=TRUE, by.x=columns, by.y=columns, suffixes = c(".2hpf", ".4hpf"))

ribodep_all_unique_246hpf_mincov <- merge(ribodep_all_unique_2hpf_4hpf_mincov, ribodep_all_unique_6hpf_mincov, all.x=TRUE, all.y=TRUE, by.x=columns, by.y=columns)

colnames(ribodep_all_unique_246hpf_mincov) <- c(colnames(ribodep_all_unique_2hpf_4hpf_mincov),"Median_Length.6hpf", "Gene_Count_Norm.6hpf")

ribodep_all_unique_246hpf_mincov_characters <- ribodep_all_unique_246hpf_mincov[,c(1,2)]
ribodep_all_unique_246hpf_mincov_numbers <- ribodep_all_unique_246hpf_mincov[,-c(1,2)]

ribodep_all_unique_246hpf_mincov_numbers[is.na(ribodep_all_unique_246hpf_mincov_numbers)] <- 0

ribodep_all_unique_246hpf_mincov_final <- cbind(ribodep_all_unique_246hpf_mincov_characters,ribodep_all_unique_246hpf_mincov_numbers )


maternal <- read.delim("264_top_maternal_decay.txt")
maternal$Group1 <- "Maternal"



ribodep_all_unique_246hpf_mincov_final2 <- merge(ribodep_all_unique_246hpf_mincov_final, maternal,all.x=TRUE, by.x="Gene_Name", by.y="Ensembl_Gene_ID")
ribodep_all_unique_246hpf_mincov_final2$Group1[is.na(ribodep_all_unique_246hpf_mincov_final2$Group1)] <- "Rest"
ribodep_all_unique_246hpf_mincov_final2$Ensembl_Transcript_ID <- NULL




#Gene profile
maternal <- read.delim("allTranscripts_riboZero_rnaSeq.maternal.txt", header=FALSE)
maternal$Group2 <- "Maternal"
zygotic <- read.delim("allTranscripts_riboZero_rnaSeq.zyfir.txt", header=FALSE)
zygotic$Group2 <- "Zygotic"

mir430 <- read.delim("allTranscripts_riboZero_rnaSeq.mir430.txt", header=FALSE)
mir430$Group2 <- "mir430"

id_convert <- read.delim("Transcript_ID_To_Gene_ID.txt")

groups <- rbind(maternal, zygotic, mir430)
groups_id <- merge(groups,id_convert, by.x="V1", by.y="Transcript.stable.ID")
groups_id$V1 <- NULL


ribodep_all_unique_246hpf_mincov_final3 <- merge(ribodep_all_unique_246hpf_mincov_final2, groups_id, all.x=TRUE, all.y=TRUE, by.x="Gene_Name", by.y="Gene.stable.ID")
ribodep_all_unique_246hpf_mincov_final3$Group2[is.na(ribodep_all_unique_246hpf_mincov_final3$Group2)] <- "Rest"





ribodep_all_unique_246hpf_mincov_taillengths <- ribodep_all_unique_246hpf_mincov_final3[,c("Gene_Name","Gene_Type","Group1","Group2", "Median_Length.2hpf","Median_Length.4hpf","Median_Length.6hpf" )]


#Normalize by 2hpf
ribodep_all_unique_246hpf_mincov_taillengths$Median_Length2.2hpf <- (ribodep_all_unique_246hpf_mincov_taillengths$Median_Length.2hpf+1)/(ribodep_all_unique_246hpf_mincov_taillengths$Median_Length.2hpf+1)
ribodep_all_unique_246hpf_mincov_taillengths$Median_Length2.4hpf <- (ribodep_all_unique_246hpf_mincov_taillengths$Median_Length.4hpf+1)/(ribodep_all_unique_246hpf_mincov_taillengths$Median_Length.2hpf+1)
ribodep_all_unique_246hpf_mincov_taillengths$Median_Length2.6hpf <- (ribodep_all_unique_246hpf_mincov_taillengths$Median_Length.6hpf+1)/(ribodep_all_unique_246hpf_mincov_taillengths$Median_Length.2hpf+1)

ribodep_all_unique_246hpf_mincov_taillengths_relativeto2hpf <- ribodep_all_unique_246hpf_mincov_taillengths[,c("Gene_Name","Gene_Type", "Group1","Group2", "Median_Length2.2hpf" ,"Median_Length2.4hpf", "Median_Length2.6hpf")]


ribodep_all_unique_246hpf_mincov_taillengths_melted_relativeto2hpf <- melt(ribodep_all_unique_246hpf_mincov_taillengths_relativeto2hpf)
ribodep_all_unique_246hpf_mincov_taillengths_melted_mrna_relativeto2hpf  <- subset(ribodep_all_unique_246hpf_mincov_taillengths_melted_relativeto2hpf, Gene_Type=="protein_coding")



pdf(file="Zebrafish_Embryos_TailLength_Line_Plot_Maternal_vs_Rest_Relative_to_2hpf.pdf",height=4,width=10,onefile=FALSE)
print(ggplot(ribodep_all_unique_246hpf_mincov_taillengths_melted_relativeto2hpf, aes(x=variable, y=log(value), group=Gene_Name)) +
  geom_line(data=subset(ribodep_all_unique_246hpf_mincov_taillengths_melted_relativeto2hpf, Group1=="Rest"), alpha=1/5, color="grey")+
  geom_line(data=subset(ribodep_all_unique_246hpf_mincov_taillengths_melted_relativeto2hpf, Group1=="Maternal"), alpha=1/2, color="red")+
  theme_bw()+
  geom_point(data=subset(ribodep_all_unique_246hpf_mincov_taillengths_melted_relativeto2hpf, Group1=="Rest"), color="grey")+
    geom_point(data=subset(ribodep_all_unique_246hpf_mincov_taillengths_melted_relativeto2hpf, Group1=="Maternal"), color="red"))

dev.off()





pdf(file="Zebrafish_Embryos_TailLength_Line_Plot_Maternal_vs_Rest_Relative_to_2hpf_3Groups.pdf",height=10,width=10,onefile=FALSE)
print(ggplot(ribodep_all_unique_246hpf_mincov_taillengths_melted_relativeto2hpf, aes(x=variable, y=log(value), group=Gene_Name)) +
  geom_line(alpha=1, aes(color=Group2))+
  theme_bw()+
  facet_wrap(~Group2)+
  geom_point( aes(color=Group2)))
dev.off()





#DONT Normalize by 2hpf

ribodep_all_unique_246hpf_mincov_taillengths_onlyabsolute <- ribodep_all_unique_246hpf_mincov_taillengths[,c("Gene_Name","Gene_Type", "Group1","Group2", "Median_Length.2hpf" ,"Median_Length.4hpf", "Median_Length.6hpf")]

ribodep_all_unique_246hpf_mincov_taillengths_melted <- melt(ribodep_all_unique_246hpf_mincov_taillengths_onlyabsolute)
ribodep_all_unique_246hpf_mincov_taillengths_melted_mrna  <- subset(ribodep_all_unique_246hpf_mincov_taillengths_melted, Gene_Type=="protein_coding")



pdf(file="Zebrafish_Embryos_TailLength_Line_Plot_Maternal_vs_Rest.pdf",height=4,width=10,onefile=FALSE)
print(ggplot(ribodep_all_unique_246hpf_mincov_taillengths_melted_mrna, aes(x=variable, y=log(value), group=Gene_Name)) +
  geom_line(data=subset(ribodep_all_unique_246hpf_mincov_taillengths_melted_mrna, Group1=="Rest"), alpha=1/5, color="grey")+
  geom_line(data=subset(ribodep_all_unique_246hpf_mincov_taillengths_melted_mrna, Group1=="Maternal"), alpha=1/2, color="red")+
  theme_bw()+
  geom_point(data=subset(ribodep_all_unique_246hpf_mincov_taillengths_melted_mrna, Group1=="Rest"), color="grey")+
    geom_point(data=subset(ribodep_all_unique_246hpf_mincov_taillengths_melted_mrna, Group1=="Maternal"), color="red"))

dev.off()





pdf(file="Zebrafish_Embryos_TailLength_Line_Plot_Maternal_vs_Rest_3Groups.pdf",height=10,width=10,onefile=FALSE)
print(ggplot(ribodep_all_unique_246hpf_mincov_taillengths_melted_mrna, aes(x=variable, y=log(value), group=Gene_Name)) +
  geom_line(alpha=1, aes(color=Group2))+
  theme_bw()+
  facet_wrap(~Group2)+
  geom_point( aes(color=Group2)))
dev.off()




pdf(file= "Zebrafish_Embryos_TailLength_Dotplot_Maternal_vs_Rest_3Groups.pdf",height=10,width=20,onefile=FALSE)
  print(ggplot(ribodep_all_unique_246hpf_mincov_taillengths_melted_mrna, aes(x=variable, y=log(value+1))) + 
    geom_quasirandom(varwidth = TRUE, aes(color=variable))+
    geom_boxplot(aes(alpha=0), outlier.shape=NA)+
    stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
            geom = "crossbar", width = 0.7, color="#c06c84")+
    theme_bw()+
    ggtitle("Zebrafish Embryo")+
    xlab("Time Points")+
    ylab("log(Median Tail Length)") +
    facet_wrap(~Group2, scales="free")+
    theme(axis.text=element_text(size=14),strip.text = element_text(size=13),
            axis.title=element_text(size=17,face="bold"),
            legend.title = element_text(size = 20),
            legend.text = element_text(color = "black", size=15)))
dev.off()








pdf(file="Zebrafish_Embryos_TailLength_Boxplot_Maternal_vs_Rest_3Groups.pdf",height=4,width=10,onefile=FALSE)
    print(ggplot(ribodep_all_unique_246hpf_mincov_taillengths_melted_mrna, aes(x = Group2, y = log(value+1) )) + 
      geom_boxplot(aes(fill = variable),position = position_dodge(0.9)) +
      #stat_n_text() + 
      #scale_fill_manual(values=colors)+
      #coord_cartesian(ylim = c(0,175))+
      theme_bw())
  dev.off()

            

 pdf(file="Zebrafish_Embryos_TailLength_Violinplot_Maternal_vs_Rest_3Groups.pdf",height=4,width=10,onefile=FALSE)
print(ggplot(ribodep_all_unique_246hpf_mincov_taillengths_melted_mrna, aes(x=Group2, y=log(value+1), fill=variable)) +
  geom_violin(position=position_dodge(1)))
dev.off()













































pdf(file= "Zebrafish_Embryos_GeneCount_Dot_Plot_Maternal_vs_Rest_Relative_to_2hpf.pdf",height=10,width=20,onefile=FALSE)
  print(ggplot(ribodep_all_unique_246hpf_genecounts_melted_mrna_relativeto2hpf, aes(x=variable, y=log(value+1))) + 
    geom_quasirandom(varwidth = TRUE, aes(color=variable))+
    geom_boxplot(aes(alpha=0), outlier.shape=NA)+
    stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
            geom = "crossbar", width = 0.7, color="#c06c84")+
    theme_bw()+
    ggtitle("Zebrafish Embryo")+
    xlab("Time Points")+
    ylab("log(Normalized Gene Count)") +
    facet_wrap(~Group, scales="free")+
    theme(axis.text=element_text(size=14),strip.text = element_text(size=13),
            axis.title=element_text(size=17,face="bold"),
            legend.title = element_text(size = 20),
            legend.text = element_text(color = "black", size=15)))
dev.off()








### DESEQ2
library(DESeq2)
ribodep_all_unique_246hpf_genecounts_mRNA <- subset(ribodep_all_unique_246hpf_genecounts, Gene_Type=="protein_coding")
deseq_count_data <- ribodep_all_unique_246hpf_genecounts_mRNA[,c("Gene_Name", "Gene_Count_Norm.2hpf", "Gene_Count_Norm.4hpf","Gene_Count_Norm.6hpf")]


dds <- DESeqDataSetFromMatrix(countData=deseq_count_data, 
                              design=~dex, tidy = TRUE)















pdf(file="MedianTails_Ribodepleted_Timepoints_Maternal.pdf",height=6,width=15,onefile=FALSE)
  print(ggplot(maternal_group, aes(x=Sample, y=Median_Length)) + 
    geom_quasirandom(varwidth = TRUE, aes(color=Sample))+
    geom_boxplot(aes(alpha=0), outlier.shape=NA)+
    stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
        geom = "crossbar", width = 0.7, color="#c06c84")+
    theme_bw()+
    xlab("Group")+
        ylab("Median Tail length") +
        facet_wrap(~Group)+
    theme(axis.text=element_text(size=14),strip.text = element_text(size=13),
        axis.title=element_text(size=17,face="bold"),
        legend.title = element_text(size = 20),
        legend.text = element_text(color = "black", size=15)))
dev.off()





#Only matching genes
maternal_group_complete <- vector()
for (gene in unique(maternal_group$Gene_Name)){
  subs <- subset(maternal_group, Gene_Name==gene)
  if (nrow(subs) ==3){
    maternal_group_complete <- rbind(maternal_group_complete, subs)
    } else {}
}





#Line Plot
pdf(file="Zebrafish_Embryos_MedianTail_Line_Plot_Maternal.pdf",height=5,width=15,onefile=FALSE)
print(ggplot(maternal_group_complete, aes(x=Sample, y=Median_Length, group=Gene_Name)) +
  geom_line(alpha=1, aes(color=Group))+
  facet_wrap(~Group)+
  geom_point( aes(color=Group)))
dev.off()



#Line Plot
pdf(file="Zebrafish_Embryos_GeneCount_Line_Plot_Maternal.pdf",height=5,width=15,onefile=FALSE)
print(ggplot(maternal_group_complete, aes(x=Sample, y=log(Gene_Count_Norm+1), group=Gene_Name)) +
  geom_line(alpha=1, aes(color=Group))+
  facet_wrap(~Group)+
  geom_point( aes(color=Group)))
dev.off()












### Count pA on mtrRNAs

hpf2_mtrRNA <- subset(ribodep_hpf2.tail_class_count, Gene_Type=="Mt_rRNA")
hpf4_mtrRNA <- subset(ribodep_hpf4.tail_class_count, Gene_Type=="Mt_rRNA")
hpf6_mtrRNA <- subset(ribodep_hpf6.tail_class_count, Gene_Type=="Mt_rRNA")

hpf4_mtrRNA_pA <- subset(polyA_hpf4.tail_class_count, Gene_Type=="Mt_rRNA")



mtrRNA_zebrafish_ribodep <- rbind(hpf2_mtrRNA, hpf4_mtrRNA,hpf6_mtrRNA,hpf4_mtrRNA_pA)
mtrRNA_zebrafish_ribodep_uniq <-mtrRNA_zebrafish_ribodep[!duplicated(mtrRNA_zebrafish_ribodep[c("Gene_Name", "Sample")]),]












