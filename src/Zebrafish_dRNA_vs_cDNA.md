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


# Import the data
#RIBODEPLETED DATA

#POLYA SELECTED
nano3p_pA_hpf4.data <- read.delim("4hpf_pAselected.genome11_sequin_ALLRNAs_Merged.bed", header=FALSE)

nano3p_ribo_hpf4.data <- read.delim("4hpf.genome11_sequin_ALLRNAs_Merged_Rep2.bed", header=FALSE)

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
  return(merged4)
}





nano3p_pA_hpf4.reshape <- reshape(nano3p_pA_hpf4.data,"Nano3Pseq_PolyA")
nano3p_ribo_hpf4.reshape <- reshape(nano3p_ribo_hpf4.data,"Nano3Pseq_Ribodep")
drna_hpf4.reshape <- reshape(drna_hpf4.data,"dRNA")

all <- rbind(nano3p_ribo_hpf4.reshape, nano3p_pA_hpf4.reshape,drna_hpf4.reshape )




#To look into rRNAs 
drna_hpf4.rRNA <- subset(drna_hpf4.reshape, Gene_Type=="rRNA")
  drna_hpf4.rRNA2 <-  drna_hpf4.rRNA [!duplicated(drna_hpf4.rRNA [c("Gene_Name", "Sample")]),]




#### Gene Type Count Comparison ###
simple_barplot_grouped <- function(data, label){
  #Remove duplicates
  data2 <-  data[!duplicated(data[c("Gene_Type", "Sample")]),]
  data2 <- subset(data2, Gene_Type != "synthetic" & Gene_Type != "ensembl_havana"  & Gene_Type != "ensembl" & Gene_Type != "havana" & Gene_Type != "polymorphic_pseudogene" & Gene_Type != "TEC" & Gene_Type != "unprocessed_pseudogene" & Gene_Type != "transcribed_unprocessed_pseudogene" & Gene_Type != "snRNA"& Gene_Type != "TR_J_gene")
  category_sum <- aggregate(.~Gene_Type+Sample, data2[,c("Gene_Type","Sample",  "Gene_Type_Count_Norm")], sum)
  category_sum <- category_sum[order(-category_sum$Gene_Type_Count_Norm),]
  category_sum$Gene_Type <- factor(category_sum$Gene_Type, levels = unique(category_sum$Gene_Type))
  pdf(file=paste(label, "Gene_Type_Normalized_Count.pdf",sep="_"),height=6,width=15,onefile=FALSE)
  print(ggplot(category_sum, aes(fill=Sample,  y=Gene_Type_Count_Norm, x=Gene_Type)) + 
    geom_bar(position="dodge", stat="identity")+
    #scale_y_continuous(trans = 'log2')+
    theme_bw())
    dev.off()
}



simple_barplot_grouped(all, "Nano3Pseq_vs_dRNA")


simple_barplot_grouped_log <- function(data, label){
  #Remove duplicates
  data2 <-  data[!duplicated(data[c("Gene_Type", "Sample")]),]
  data2 <- subset(data2, Gene_Type != "synthetic" & Gene_Type != "ensembl_havana"  & Gene_Type != "ensembl" & Gene_Type != "havana" & Gene_Type != "polymorphic_pseudogene" & Gene_Type != "TEC" & Gene_Type != "unprocessed_pseudogene" & Gene_Type != "transcribed_unprocessed_pseudogene" & Gene_Type != "snRNA" & Gene_Type != "TR_J_gene" & Gene_Type!= "Mt_tRNA" & Gene_Type!= "sense_intronic"& Gene_Type!= "snoRNA"& Gene_Type!= "ribozyme")
  category_sum <- aggregate(.~Gene_Type+Sample, data2[,c("Gene_Type","Sample",  "Gene_Type_Count_Norm")], sum)
  category_sum <- category_sum[order(-category_sum$Gene_Type_Count_Norm),]
  category_sum$Gene_Type <- factor(category_sum$Gene_Type, levels = unique(category_sum$Gene_Type))
  pdf(file=paste(label, "logGene_Type_Normalized_Count.pdf",sep="_"),height=6,width=15,onefile=FALSE)
  print(ggplot(category_sum, aes(fill=Sample,  y=log(Gene_Type_Count_Norm+1), x=Gene_Type)) + 
    geom_bar(position="dodge", stat="identity")+
    #scale_y_continuous(trans = 'log2')+
    theme_bw())
    dev.off()
}
simple_barplot_grouped_log(all, "Nano3Pseq_vs_dRNA")






  all_gene_types <-  all[!duplicated(all[c("Gene_Type", "Sample")]),]
  all_gene_types_coding <- subset(all_gene_types, Gene_Type=="protein_coding")
  all_gene_types_coding$Category_New <- "Coding"
  all_gene_types_noncoding <- subset(all_gene_types, Gene_Type!="protein_coding")
  all_gene_types_noncoding$Category_New <- "Noncoding"

  coding_sum <- aggregate(.~Category_New+Sample, all_gene_types_coding[,c("Category_New","Sample",  "Gene_Type_Count_Norm")], sum)
  noncoding_sum <- aggregate(.~Category_New+Sample, all_gene_types_noncoding[,c("Category_New","Sample",  "Gene_Type_Count_Norm")], sum)

  both_sum <- rbind(noncoding_sum,coding_sum)

both_sum$Category_New <- factor(both_sum$Category_New, levels = both_sum$Category_New)
  pdf(file="dRNA_cDNA_comparison_gene_type_barplot.pdf",height=6,width=6,onefile=FALSE)
print(ggplot(both_sum, aes(fill=Category_New, y=Gene_Type_Count_Norm, x=Sample)) + 
    geom_bar(position="stack", stat="identity"))
dev.off()









nano3p_pA_hpf4.unique <-  nano3p_pA_hpf4.reshape[!duplicated(nano3p_pA_hpf4.reshape[c("Gene_Name", "Sample")]),]
drna_hpf4.unique <-  drna_hpf4.reshape[!duplicated(drna_hpf4.reshape[c("Gene_Name", "Sample")]),]



columns1 <- c("Gene_Type", "Gene_Name","Gene_Count", "Gene_Count_Norm")
columns2 <- c("Gene_Type",  "Gene_Name")

nano3p_vs_dRNA_merged <- merge(nano3p_pA_hpf4.unique[,columns1] ,drna_hpf4.unique[,columns1], by.x=columns2, by.y=columns2  )



nano3p_vs_dRNA_merged_protein <- subset(nano3p_vs_dRNA_merged, Gene_Type=="protein_coding")
nano3p_vs_dRNA_merged_protein2 <- subset(nano3p_vs_dRNA_merged_protein, Gene_Count.x>20 & Gene_Count.y>20 )


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


plot_denscols_with_corr_pearson("Nano3P_vs_dRNA_Gene_Count", nano3p_vs_dRNA_merged_protein2$Gene_Count_Norm.x, nano3p_vs_dRNA_merged_protein2$Gene_Count_Norm.y, "Nano3P-seq", "dRNA" )



