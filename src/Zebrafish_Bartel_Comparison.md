########################################################
######## ZEBRAFISH Bartel Comparison ####
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






