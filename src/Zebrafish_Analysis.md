##################################################
######## ANALYSIS OF THE ZEBRAFISH RUNS #########
###################################################
########## OGUZHAN BEGIK APRIL 2020 ###############

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







# Import the data
#RIBODEPLETED DATA
ribodep_hpf2_rep1.data <- read.delim("2hpf.genome11_sequin_ALLRNAs_Merged_Rep1.bed", header=FALSE)
ribodep_hpf4_rep1.data <- read.delim("4hpf.genome11_sequin_ALLRNAs_Merged_Rep1.bed", header=FALSE)
ribodep_hpf6_rep1.data <- read.delim("6hpf.genome11_sequin_ALLRNAs_Merged_Rep1.bed", header=FALSE)

ribodep_hpf2_rep2.data <- read.delim("2hpf.genome11_sequin_ALLRNAs_Merged_Rep2.bed", header=FALSE)
ribodep_hpf4_rep2.data <- read.delim("4hpf.genome11_sequin_ALLRNAs_Merged_Rep2.bed", header=FALSE)
ribodep_hpf6_rep2.data <- read.delim("6hpf.genome11_sequin_ALLRNAs_Merged_Rep2.bed", header=FALSE)


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
	gene_type_major <- subset(gene_type_sum, Gene_Count_Norm > 5.5)
	gene_type_major$Category <- gene_type_major$Gene_Type
	gene_type_minor <- subset(gene_type_sum, Gene_Count_Norm < 5.5)
	gene_type_minor$Category <- "Other"
	gene_type <- rbind(gene_type_major, gene_type_minor)
	colnames(gene_type) <- c("Gene_Type", "Gene_Type_Count_Norm", "Category")
	#MErge
	merged4 <- merge(merged3, gene_type, by.x="Gene_Type", by.y="Gene_Type")
	return(merged4)
}





ribodep_hpf2_rep1.reshape <- reshape(ribodep_hpf2_rep1.data,ribodepleted_rep1.tails_processed,"Ribodep_2hpf_rep1")
ribodep_hpf4_rep1.reshape <- reshape(ribodep_hpf4_rep1.data,ribodepleted_rep1.tails_processed,"Ribodep_4hpf_rep1")
ribodep_hpf6_rep1.reshape <- reshape(ribodep_hpf6_rep1.data,ribodepleted_rep1.tails_processed,"Ribodep_6hpf_rep1")



ribodep_hpf2_rep2.reshape <- reshape(ribodep_hpf2_rep2.data,ribodepleted_rep2.tails_processed,"Ribodep_2hpf_rep2")
ribodep_hpf4_rep2.reshape <- reshape(ribodep_hpf4_rep2.data,ribodepleted_rep2.tails_processed,"Ribodep_4hpf_rep2")
ribodep_hpf6_rep2.reshape <- reshape(ribodep_hpf6_rep2.data,ribodepleted_rep2.tails_processed,"Ribodep_6hpf_rep2")


polyA_hpf4.reshape <- reshape(polyA_hpf4.data,polyA.tails_processed,"PolyA_4hpf")





polyA_vs_ribodep_4hpf <- rbind(ribodep_hpf4_rep1.reshape,polyA_hpf4.reshape)

ribodep_rep1_all <- rbind(ribodep_hpf2_rep1.reshape, ribodep_hpf4_rep1.reshape, ribodep_hpf6_rep1.reshape)
ribodep_rep2_all <- rbind(ribodep_hpf2_rep2.reshape, ribodep_hpf4_rep2.reshape, ribodep_hpf6_rep2.reshape)

ribodep_hpf2_replicability <- rbind(ribodep_hpf2_rep1.reshape,ribodep_hpf2_rep2.reshape) 
ribodep_hpf4_replicability <- rbind(ribodep_hpf4_rep1.reshape,ribodep_hpf4_rep2.reshape) 
ribodep_hpf6_replicability <- rbind(ribodep_hpf6_rep1.reshape,ribodep_hpf6_rep2.reshape) 


ribodep_rep1_all$Rep <- "Rep1"
ribodep_rep2_all$Rep <- "Rep2"

ribodep_rep1_rep2_all <- rbind(ribodep_rep1_all,ribodep_rep2_all)




#### Gene Type Count Comparison ###
simple_barplot_grouped <- function(data, label){
	#Remove duplicates
	data2 <-  data[!duplicated(data[c("Gene_Type", "Sample")]),]
	data2 <- subset(data2, Gene_Type != "synthetic" & Gene_Type != "ensembl_havana"  & Gene_Type != "ensembl" & Gene_Type != "havana" & Gene_Type != "polymorphic_pseudogene" & Gene_Type != "TEC" & Gene_Type != "unprocessed_pseudogene" & Gene_Type != "transcribed_unprocessed_pseudogene" & Gene_Type != "snRNA"& Gene_Type != "TR_J_gene")
	#Aggregate by Gene Type
	pdf(file=paste(label, "Gene_Type_Normalized_Count.pdf",sep="_"),height=6,width=20,onefile=FALSE)
	print(ggplot(data2, aes(fill=Sample, y=log(Gene_Type_Count_Norm+1), x=Gene_Type)) + 
    geom_bar(position="dodge", stat="identity")+
    theme_bw())
    dev.off()
}



simple_barplot_grouped(polyA_vs_ribodep_4hpf, "PolyA_vs_Ribodep_4HPF")
simple_barplot_grouped(ribodep_rep1_all, "Ribodepletion_All_Timepoints_rep1")
simple_barplot_grouped(ribodep_rep2_all, "Ribodepletion_All_Timepoints_rep2")

simple_barplot_grouped(ribodep_hpf2_replicability, "Ribodepletion_2hpf_Replicability")
simple_barplot_grouped(ribodep_hpf4_replicability, "Ribodepletion_4hpf_Replicability")
simple_barplot_grouped(ribodep_hpf6_replicability, "Ribodepletion_6hpf_Replicability")




##Overall Tail comparison (Single transcript)
tail_comparison_overall <- function(data, label) {
	data2 <- subset(data, Gene_Type=="lincRNA"|Gene_Type=="Mt_rRNA"|Gene_Type=="protein_coding"|Gene_Type=="rRNA" )
	pdf(file=paste(label, "Overall_Tail_Comparison_Single_Transcript.pdf",sep="_"),height=6,width=10,onefile=FALSE)
	print(ggplot(data2, aes(x=tail_length, color=Sample)) +
		geom_density()+
  		theme_bw()+
  		xlim(-10, 350)+
   		facet_wrap(~Gene_Type, scales="free"))
	dev.off()
}


tail_comparison_overall(polyA_vs_ribodep_4hpf, "PolyA_vs_Ribodep_4HPF")
tail_comparison_overall(ribodep_rep1_all, "Ribodepletion_All_Timepoints_rep1")
tail_comparison_overall(ribodep_rep2_all, "Ribodepletion_All_Timepoints_rep2")

tail_comparison_overall(ribodep_hpf2_replicability, "Ribodepletion_2hpf_Replicability")
tail_comparison_overall(ribodep_hpf4_replicability, "Ribodepletion_4hpf_Replicability")
tail_comparison_overall(ribodep_hpf6_replicability, "Ribodepletion_6hpf_Replicability")


##Overall Tail comparison (Single transcript) Both Rep1 and Rep2 
tail_comparison_overall_rep_merged<- function(data, label) {
	data2 <- subset(data, Gene_Type=="lincRNA"|Gene_Type=="Mt_rRNA"|Gene_Type=="protein_coding"|Gene_Type=="rRNA" )
	data2$Timepoint <- data2$Sample
	data2$Timepoint <- gsub("Ribodep_2hpf_rep1", "2hpf",data2$Timepoint )
	data2$Timepoint <- gsub("Ribodep_4hpf_rep1", "4hpf",data2$Timepoint )
	data2$Timepoint <- gsub("Ribodep_6hpf_rep1", "6hpf",data2$Timepoint )
	data2$Timepoint <- gsub("Ribodep_2hpf_rep2", "2hpf",data2$Timepoint )
	data2$Timepoint <- gsub("Ribodep_4hpf_rep2", "4hpf",data2$Timepoint )
	data2$Timepoint <- gsub("Ribodep_6hpf_rep2", "6hpf",data2$Timepoint )
	pdf(file=paste(label, "Overall_Tail_Comparison_Single_Transcript.pdf",sep="_"),height=6,width=10,onefile=FALSE)
	print(ggplot(data2, aes(x=tail_length, color=Timepoint)) +
		geom_density(data=subset(data2, Rep=="Rep1"), linetype="dashed")+
		geom_density(data=subset(data2, Rep=="Rep2"))+
  		theme_bw()+
  		xlim(-10, 350)+
   		facet_wrap(~Gene_Type, scales="free"))
	dev.off()
}
tail_comparison_overall_rep_merged(ribodep_rep1_rep2_all, "Ribodepletion_All_Timepoints_rep1_2")





## Tail comparison median per gene
tail_comparison_median_per_gene_protein<- function(data, label) {
	data2 <- subset(data, Gene_Type=="protein_coding")
	data2 <-  data2[!duplicated(data2[c("Gene_Name", "Sample")]),]
	pdf(file=paste(label, "Median_Tail_Per_Gene_Comparison_mRNA.pdf",sep="_"),height=6,width=10,onefile=FALSE)
	print(ggplot(data2, aes(x=Median_Length, color=Sample)) +
		geom_density()+
  		theme_bw()+
  		coord_cartesian(xlim=c(-10, 200))+
   		facet_wrap(~Gene_Type, scales="free"))
	dev.off()
}


tail_comparison_median_per_gene_protein(polyA_vs_ribodep_4hpf, "PolyA_vs_Ribodep_4HPF")
tail_comparison_median_per_gene_protein(ribodep_rep1_all, "Ribodepletion_All_Timepoints_rep1")
tail_comparison_median_per_gene_protein(ribodep_rep2_all, "Ribodepletion_All_Timepoints_rep2")

tail_comparison_median_per_gene_protein(ribodep_hpf2_replicability, "Ribodepletion_2hpf_Replicability")
tail_comparison_median_per_gene_protein(ribodep_hpf4_replicability, "Ribodepletion_4hpf_Replicability")
tail_comparison_median_per_gene_protein(ribodep_hpf6_replicability, "Ribodepletion_6hpf_Replicability")




tail_comparison_median_per_gene_protein_rep_merged <- function(data, label) {
	data2 <- subset(data, Gene_Type=="protein_coding")
	data2 <-  data2[!duplicated(data2[c("Gene_Name", "Sample")]),]
	data2$Timepoint <- data2$Sample
	data2$Timepoint <- gsub("Ribodep_2hpf_rep1", "2hpf",data2$Timepoint )
	data2$Timepoint <- gsub("Ribodep_4hpf_rep1", "4hpf",data2$Timepoint )
	data2$Timepoint <- gsub("Ribodep_6hpf_rep1", "6hpf",data2$Timepoint )
	data2$Timepoint <- gsub("Ribodep_2hpf_rep2", "2hpf",data2$Timepoint )
	data2$Timepoint <- gsub("Ribodep_4hpf_rep2", "4hpf",data2$Timepoint )
	data2$Timepoint <- gsub("Ribodep_6hpf_rep2", "6hpf",data2$Timepoint )
	pdf(file=paste(label, "Median_Tail_Per_Gene_Comparison_mRNA.pdf",sep="_"),height=6,width=10,onefile=FALSE)
	print(ggplot(data2, aes(x=Median_Length, color=Timepoint)) +
		geom_density(data=subset(data2, Rep=="Rep1"), linetype="dashed")+
		geom_density(data=subset(data2, Rep=="Rep2"))+
		theme_bw()+
  		coord_cartesian(xlim=c(-10, 200))+
   		facet_wrap(~Gene_Type, scales="free"))
	dev.off()
}

tail_comparison_median_per_gene_protein_rep_merged(ribodep_rep1_rep2_all, "Ribodepletion_All_Timepoints_rep1_2")




### Rep1 vs Rep2


merge_two_rep_tail <- function(rep1, rep2,label) {
	rep1_unique <-   rep1[!duplicated(rep1[c("Gene_Name", "Sample")]),]
	rep2_unique <-   rep2[!duplicated(rep2[c("Gene_Name", "Sample")]),]
	data_merged <- merge(rep1_unique,rep2_unique, by.x="Gene_Name", by.y="Gene_Name" )
	data_merged_mRNA <-subset(data_merged, Gene_Type.x =="protein_coding")
	data_merged_mRNA_min20 <- subset(data_merged_mRNA, Gene_Count.x > 30 & Gene_Count.y >30)
	 corr <- cor.test(data_merged_mRNA_min20$Median_Length.x, data_merged_mRNA_min20$Median_Length.y, method = "pearson", conf.level = 0.95)
	 value <- as.numeric(corr$estimate)
      pdf(file=paste(label,"Rep1vsRep2_TailMedian_Min30.pdf",sep="_"),height=4,width=4,onefile=FALSE)
		print(ggplot(data_merged_mRNA_min20, aes(x=Median_Length.x, y=Median_Length.y)) + 
		theme_bw()+ 
       ylim(0,160)+
       xlim(0,160)+
        ggtitle(paste("Rep1 vs Rep2", label))+
        annotate(geom="text", x=50, y=150, label=paste("Pearson Correlation =", value),
              color="red", size=2)+
		xlab("Rep1")+
        ylab("Rep2")+
  		geom_point())
  		dev.off()
}
merge_two_rep_tail(ribodep_hpf2_rep1.reshape,ribodep_hpf2_rep2.reshape, "2HPF")
merge_two_rep_tail(ribodep_hpf4_rep1.reshape,ribodep_hpf4_rep2.reshape, "4HPF")
merge_two_rep_tail(ribodep_hpf6_rep1.reshape,ribodep_hpf6_rep2.reshape, "6HPF")




merge_two_rep_genecount <- function(rep1, rep2,label) {
	rep1_unique <-   rep1[!duplicated(rep1[c("Gene_Name", "Sample")]),]
	rep2_unique <-   rep2[!duplicated(rep2[c("Gene_Name", "Sample")]),]
	data_merged <- merge(rep1_unique,rep2_unique, by.x="Gene_Name", by.y="Gene_Name" )
	data_merged_mRNA <-subset(data_merged, Gene_Type.x =="protein_coding")
	data_merged_mRNA_min20 <- subset(data_merged_mRNA, Gene_Count.x > 30 & Gene_Count.y >30)
	 corr <- cor.test(data_merged_mRNA_min20$Gene_Count.x, data_merged_mRNA_min20$Gene_Count.y, method = "pearson", conf.level = 0.95)
	 value <- as.numeric(corr$estimate)
      pdf(file=paste(label,"Rep1vsRep2_Genecount_Min30.pdf",sep="_"),height=4,width=4,onefile=FALSE)
		print(ggplot(data_merged_mRNA_min20, aes(x=log(Gene_Count.x), y=log(Gene_Count.y))) + 
		theme_bw()+ 
        ylim(3,7.5)+
        xlim(3,7.5)+
        ggtitle(paste("Rep1 vs Rep2", label))+
        annotate(geom="text", x=4, y=7, label=paste("Pearson Correlation =", value),
              color="red", size=2)+
		xlab("log Gene Count Rep1")+
        ylab("log Gene Count Rep2")+
  		geom_point())
  		dev.off()
}
merge_two_rep_genecount(ribodep_hpf2_rep1.reshape,ribodep_hpf2_rep2.reshape, "2HPF")
merge_two_rep_genecount(ribodep_hpf4_rep1.reshape,ribodep_hpf4_rep2.reshape, "4HPF")
merge_two_rep_genecount(ribodep_hpf6_rep1.reshape,ribodep_hpf6_rep2.reshape, "6HPF")






##### DENS COLS

merge_two_rep <- function(rep1, rep2) {
	rep1_unique <-   rep1[!duplicated(rep1[c("Gene_Name", "Sample")]),]
	rep2_unique <-   rep2[!duplicated(rep2[c("Gene_Name", "Sample")]),]
	data_merged <- merge(rep1_unique,rep2_unique, by.x="Gene_Name", by.y="Gene_Name" )
	data_merged_mRNA <-subset(data_merged, Gene_Type.x =="protein_coding")
	data_merged_mRNA_min30 <- subset(data_merged_mRNA, Gene_Count.x > 30 & Gene_Count.y >30)
	return(data_merged_mRNA_min30)
}

data_merged_mRNA_min30_2hpf <- merge_two_rep(ribodep_hpf2_rep1.reshape,ribodep_hpf2_rep2.reshape)
data_merged_mRNA_min30_4hpf <- merge_two_rep(ribodep_hpf4_rep1.reshape,ribodep_hpf4_rep2.reshape)
data_merged_mRNA_min30_6hpf <- merge_two_rep(ribodep_hpf6_rep1.reshape,ribodep_hpf6_rep2.reshape)





plot_denscols_with_corr_pearson<-function(pdfname,my_x,my_y,xlab,ylab) {
	pdf(file=paste(pdfname, "_pearson.pdf",sep=""), height=6, width=6)
	dcols<-densCols(my_x,my_y, colramp=colorRampPalette(blues9[-(1:3)]))
	plot(my_x,my_y,col=dcols,cex=1, cex.lab=1,cex.main=3,lwd=5,pch=20,xlab=xlab,ylab=ylab)
	title(main=pdfname, col.main="black", font.main=4)
	#abline(v=0, lty=2)
	# Correlation
	test<-cor.test(my_x,my_y, method="pearson")
	print(test)
	cor222<-paste("Pearson's rho =",round(as.numeric(test$estimate),3))
	#pval<-paste("Pval =",test$p.value)
	mtext(paste(cor222))
	#mtext(paste(cor222,pval,sep=" ; ")) #Print the subtitle with the dataset correlation
	dev.off()
}

plot_denscols_with_corr_pearson("Rep1_Rep2_2hpf_genecount", log(data_merged_mRNA_min30_2hpf$Gene_Count.x) , log(data_merged_mRNA_min30_2hpf$Gene_Count.y), "2HPF_Rep1_logGeneCount", "2HPF_Rep2_logGeneCount" )
plot_denscols_with_corr_pearson("Rep1_Rep2_4hpf_genecount", log(data_merged_mRNA_min30_4hpf$Gene_Count.x) , log(data_merged_mRNA_min30_4hpf$Gene_Count.y), "4HPF_Rep1_logGeneCount", "4HPF_Rep2_logGeneCount" )

plot_denscols_with_corr_pearson("Rep1_Rep2_6hpf_genecount", log(data_merged_mRNA_min30_6hpf$Gene_Count.x) , log(data_merged_mRNA_min30_6hpf$Gene_Count.y), "6HPF_Rep1_logGeneCount", "6HPF_Rep2_logGeneCount" )





plot_denscols_with_corr_pearson_tail<-function(pdfname,my_x,my_y,xlab,ylab) {
	pdf(file=paste(pdfname, "_pearson.pdf",sep=""), height=6, width=6)
	dcols<-densCols(my_x,my_y, colramp=colorRampPalette(blues9[-(1:3)]))
	plot(my_x,my_y,col=dcols,cex=1, cex.lab=1,cex.main=3,lwd=5,pch=20,xlab=xlab,ylab=ylab,xlim = c(0,160),ylim = c(0,160) )
	title(main=pdfname, col.main="black", font.main=4)
	#abline(v=0, lty=2)
	# Correlation
	test<-cor.test(my_x,my_y, method="pearson")
	print(test)
	cor222<-paste("Pearson's rho =",round(as.numeric(test$estimate),3))
	#pval<-paste("Pval =",test$p.value)
	mtext(paste(cor222))
	#mtext(paste(cor222,pval,sep=" ; ")) #Print the subtitle with the dataset correlation
	dev.off()
}




plot_denscols_with_corr_pearson_tail("Rep1_Rep2_2hpf_tail", data_merged_mRNA_min30_2hpf$Median_Length.x , data_merged_mRNA_min30_2hpf$Median_Length.y, "2HPF_Rep1_Median_Length", "2HPF_Rep2_Median_Length" )


plot_denscols_with_corr_pearson_tail("Rep1_Rep2_4hpf_tail", data_merged_mRNA_min30_4hpf$Median_Length.x , data_merged_mRNA_min30_4hpf$Median_Length.y, "4HPF_Rep1_Median_Length", "4HPF_Rep2_Median_Length" )


plot_denscols_with_corr_pearson_tail("Rep1_Rep2_6hpf_tail", data_merged_mRNA_min30_6hpf$Median_Length.x , data_merged_mRNA_min30_6hpf$Median_Length.y, "6HPF_Rep1_Median_Length", "6HPF_Rep2_Median_Length" )











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

bartel_hpf2.processed <- bartel_merging(bartel_hpf2.data, ribodep_hpf2.reshape, 20)
bartel_hpf4.processed <- bartel_merging(bartel_hpf4.data, ribodep_hpf4.reshape, 20)
bartel_hpf6.processed <- bartel_merging(bartel_hpf6.data, ribodep_hpf6.reshape, 20)



#SCATTER PLOT 
scatter_bartel_comparison <- function(data, label){
	 corr <- cor.test(data$Median_Length, data$Median.T, method = "pearson", conf.level = 0.95)
	 value <- as.numeric(corr$estimate)
      pdf(file=paste(label,"Bartel_vs_OurStudy_TailMedian_Min20.pdf", sep="_"),height=4,width=4,onefile=FALSE)
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




# COMPARE POLYA VS RIBODEPLETED IN 


ribodep_hpf4_unique <-   ribodep_hpf4.reshape[!duplicated(ribodep_hpf4.reshape[c("Gene_Name", "Sample")]),]
polyA_hpf4_unique <-   polyA_hpf4.reshape[!duplicated(polyA_hpf4.reshape[c("Gene_Name", "Sample")]),]
hpf4_merged <- merge(ribodep_hpf4_unique,polyA_hpf4_unique, by.x="Gene_Name", by.y="Gene_Name" )
hpf4_merged_mRNA <-subset(hpf4_merged, Gene_Type.x =="protein_coding")
hpf4_merged_mRNA_min20 <- subset(hpf4_merged_mRNA, Gene_Count.x > 20 & Gene_Count.y >20)

	 corr <- cor.test(hpf4_merged_mRNA_min20$Median_Length.x, hpf4_merged_mRNA_min20$Median_Length.y, method = "pearson", conf.level = 0.95)
	 value <- as.numeric(corr$estimate)
      pdf(file="PolyA_vs_Ribodep_4hpf_TailMedian_Min20.pdf",height=4,width=4,onefile=FALSE)
		print(ggplot(hpf4_merged_mRNA_min20, aes(x=Median_Length.x, y=Median_Length.y)) + 
		theme_bw()+ 
        ylim(0,150)+
        xlim(0,150)+
        ggtitle("PolyA vs Ribodepletion 4hpf")+
        annotate(geom="text", x=50, y=150, label=paste("Pearson Correlation =", value),
              color="red", size=2)+
		xlab("Ribodepletion")+
        ylab("PolyA Selection")+
  		geom_point())
  		dev.off()



hpf4_merged_mRNA_min20$Abs_Difference <- abs(hpf4_merged_mRNA_min20$Median_Length.x- hpf4_merged_mRNA_min20$Median_Length.y)
hpf4_merged_mRNA_min20_significant_changing <- subset(hpf4_merged_mRNA_min20, Abs_Difference > 50)




### BOXPLOT OF THE MEDIAN TAIL PER GENE IN RIBODEPLETED RUN
	GENE_ID <- read.delim("Zebrafish_ID_Name_Conversion.txt")
	ribodep_all_unique <-ribodep_rep2_all[!duplicated(ribodep_rep2_all[c("Gene_Name", "Sample")]),]
	ribodep_all_unique_names <- merge(ribodep_all_unique, GENE_ID, by.x="Gene_Name", by.y="Gene.stable.ID")
	ribodep_all_unique_min20 <-  subset(ribodep_all_unique_names, Gene_Count > 20)
	ribodep_all_unique_min20_mRNA <- subset(ribodep_all_unique_min20, Gene_Type =="protein_coding")
	

	
	
		pdf(file= "Tails_Ribodepleted_Timepoints_ProteinCoding_Tails_Min20_rep2.pdf",height=10,width=20,onefile=FALSE)
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
	








#### MERGE THE RIBODEPLETED RUNS WITH THEIR GROUP OF GENES

#Gene profile
maternal <- read.delim("allTranscripts_riboZero_rnaSeq.maternal.txt", header=FALSE)
maternal$Group <- "Maternal"
zygotic <- read.delim("allTranscripts_riboZero_rnaSeq.zyfir.txt", header=FALSE)
zygotic$Group <- "Zygotic"

mir430 <- read.delim("allTranscripts_riboZero_rnaSeq.mir430.txt", header=FALSE)
mir430$Group <- "mir430"

	id_convert <- read.delim("Transcript_ID_To_Gene_ID.txt")

groups <- rbind(maternal, zygotic, mir430)
groups_id <- merge(groups,id_convert, by.x="V1", by.y="Transcript.stable.ID")



merged_protein_group <- merge(ribodep_all_unique_min20_mRNA, groups_id, by.x="Gene_Name", by.y="Gene.stable.ID")








pdf(file="MedianTails_Ribodepleted_Timepoints_ProteinCoding_byGroup.pdf",height=6,width=15,onefile=FALSE)
	print(ggplot(merged_protein_group, aes(x=Sample, y=Median_Length)) + 
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
merged_protein_group_complete <- vector()
for (gene in unique(merged_protein_group$Gene_Name)){
	subs <- subset(merged_protein_group, Gene_Name==gene)
	if (nrow(subs) ==3){
		merged_protein_group_complete <- rbind(merged_protein_group_complete, subs)
		} else {}
}





#Line Plot
pdf(file="Zebrafish_Embryos_MedianTail_Line_Plot_Grouped.pdf",height=5,width=15,onefile=FALSE)
print(ggplot(merged_protein_group_complete, aes(x=Sample, y=Median_Length, group=Gene_Name)) +
  geom_line(alpha=1, aes(color=Group))+
  facet_wrap(~Group)+
  geom_point( aes(color=Group)))
dev.off()



#Line Plot
pdf(file="Zebrafish_Embryos_GeneCount_Line_Plot_Grouped.pdf",height=5,width=15,onefile=FALSE)
print(ggplot(merged_protein_group_complete, aes(x=Sample, y=log(Gene_Count_Norm+1), group=Gene_Name)) +
  geom_line(alpha=1, aes(color=Group))+
  facet_wrap(~Group)+
  geom_point( aes(color=Group)))
dev.off()






      pdf(file="Zebrafish_Embryos_MedianTail_Count_Normalized_Scatter.pdf",height=4,width=12,onefile=FALSE)
		print(ggplot(merged_protein_group, aes(x=log(Gene_Count_Norm+1), y=Median_Length,color=Group)) + 
		theme_bw()+ 
		facet_wrap(~Sample)+
		xlab("Log(Normalized Count)")+
        ylab("Median Tail Length")+
  		geom_point())
  		dev.off()





#### CHECKING THE RIBOSOMAL RNAS
##1 . OVERALL TAIL LENGTH DISTRIBUTION 
	ribodep_all_rRNA <- subset(ribodep_all, Gene_Name=="18s_Maternal" | Gene_Name=="28s_Maternal"| Gene_Name=="ENSDARG00000080337" | Gene_Name=="ENSDARG00000082753" )
pdf(file="Tails_Ribodepleted_Timepoints_rRNA.pdf",height=20,width=10,onefile=FALSE)
	print(ggplot(ribodep_all_rRNA, aes(x=Sample, y=tail_length)) + 
		geom_quasirandom(varwidth = TRUE, aes(color=Sample))+
		geom_boxplot(aes(alpha=0), outlier.shape=NA)+
		stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
    		geom = "crossbar", width = 0.7, color="#c06c84")+
		theme_bw()+
		xlab("Group")+
      	ylab("Tail length") +
      	facet_wrap(~Gene_Name)+
		theme(axis.text=element_text(size=14),strip.text = element_text(size=13),
    		axis.title=element_text(size=17,face="bold"),
    		legend.title = element_text(size = 20),
    		legend.text = element_text(color = "black", size=15)))
dev.off()

























###### CHECK THE SHORT/LONG RATIO PER GENE




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

ribodep_hpf2.tail_class_count <- count_polyA_population(ribodep_hpf2.reshape)
ribodep_hpf4.tail_class_count <- count_polyA_population(ribodep_hpf4.reshape)
ribodep_hpf6.tail_class_count <- count_polyA_population(ribodep_hpf6.reshape)







#Scatter plot with different hours 
ribodep_hpf2.tail_class_count_unique <-ribodep_hpf2.tail_class_count[!duplicated(ribodep_hpf2.tail_class_count[c("Gene_Name", "Sample")]),]
ribodep_hpf4.tail_class_count_unique <-ribodep_hpf4.tail_class_count[!duplicated(ribodep_hpf4.tail_class_count[c("Gene_Name", "Sample")]),]


columns <- c("Gene_Name", "Gene_Type", "Gene_Count", "Ratio_PolyA", "Ratio_NoPolyA")
merged_2_4_hours <-  merge(ribodep_hpf2.tail_class_count_unique[,columns], ribodep_hpf4.tail_class_count_unique[,columns], by.x=c("Gene_Name", "Gene_Type"),  by.y=c("Gene_Name", "Gene_Type"))

merged_2_4_hours_mRNA <- subset(merged_2_4_hours, Gene_Type=="protein_coding")
merged_2_4_hours_mRNA$Change <- abs(merged_2_4_hours_mRNA$Ratio_PolyA.x-merged_2_4_hours_mRNA$Ratio_PolyA.y)

scatter_plot <- function(data,label1, label2){
	pdf(file=paste(label1,label2, "PolyA_Ratio_Scatter_Plot_mRNAs.pdf", sep="_"),height=5,width=8,onefile=FALSE)
	print(ggplot(data, aes(x=Ratio_PolyA.x, y=Ratio_PolyA.y, colour=Change)) +
			geom_point(size=1)+
			geom_abline(slope=1, intercept=0,linetype="dashed", size=0.2, color= "black")+
			#geom_text_repel(data=subset(data, Change>0.5), aes(label=Gene_Name), colour="black",segment.size  = 0.4,segment.color = "grey50",size=5)+
			xlim(0,1)+
			ylim(0,1)+
			xlab(paste(label1, "Long PolyA Fraction"))+
			ylab(paste(label2,  "Long PolyA Fraction")) +
			theme_bw()+
			theme(axis.text.x = element_text(face="bold", color="black",size=11),
				 axis.text.y = element_text(face="bold", color="black", size=11),
			plot.title = element_text(color="black", size=24, face="bold.italic",hjust = 0.5),
			axis.title.x = element_text(color="black", size=15, face="bold"),
			axis.title.y = element_text(color="black", size=15, face="bold"),
			panel.background = element_blank(),
			axis.line = element_line(colour = "black", size=0.5),
			legend.title = element_text(color = "black", size = 20,face="bold"),
			legend.text = element_text(color = "black", size=20)))
dev.off()
}

scatter_plot(merged_2_4_hours_mRNA, "2HPF", "4HPF")


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










#### CHECK BINOMIALITY

library(diptest)
bimodality_test <- function(data) {
	data_cov <- subset(data, Gene_Count >50)
	data_mRNA <- subset(data_cov, Gene_Type=="protein_coding")
	bimodal_table <- vector()
	for (gene in unique(data_mRNA$Gene_Name)) {
		subs <- subset(data_mRNA, Gene_Name==gene)
		test <- (dip.test(subs$tail_length, simulate.p.value = FALSE, B = 2000))
		statistic <- as.numeric(as.character(test$statistic))
		p_value <- as.numeric(as.character(test$p.value))
		gene_name <- gene
		data_table <- as.data.frame(cbind(gene_name,statistic, p_value))
		data_table$statistic <- as.numeric(as.character(data_table$statistic))
		data_table$p_value <- as.numeric(as.character(data_table$p_value))
		bimodal_table <- rbind(bimodal_table,data_table)

	}
	return(bimodal_table)
}




ribodep_hpf2.bimodality <- bimodality_test(ribodep_hpf2.reshape)






