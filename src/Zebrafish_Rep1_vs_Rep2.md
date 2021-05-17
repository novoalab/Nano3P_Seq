########################################################
##### ANALYSIS OF THE ZEBRAFISH RUNS REP COMPARISON ####
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


# Import the data
#RIBODEPLETED DATA
ribodep_hpf2_rep1.data <- read.delim("2hpf.genome11_sequin_ALLRNAs_Merged_Rep1.bed", header=FALSE)
ribodep_hpf4_rep1.data <- read.delim("4hpf.genome11_sequin_ALLRNAs_Merged_Rep1.bed", header=FALSE)
ribodep_hpf6_rep1.data <- read.delim("6hpf.genome11_sequin_ALLRNAs_Merged_Rep1.bed", header=FALSE)

ribodep_hpf2_rep2.data <- read.delim("2hpf.genome11_sequin_ALLRNAs_Merged_Rep2.bed", header=FALSE)
ribodep_hpf4_rep2.data <- read.delim("4hpf.genome11_sequin_ALLRNAs_Merged_Rep2.bed", header=FALSE)
ribodep_hpf6_rep2.data <- read.delim("6hpf.genome11_sequin_ALLRNAs_Merged_Rep2.bed", header=FALSE)



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


ribodep_hpf2_rep1.reshape <- reshape(ribodep_hpf2_rep1.data,ribodepleted_rep1.tails_processed,"Ribodep_2hpf_rep1")
ribodep_hpf4_rep1.reshape <- reshape(ribodep_hpf4_rep1.data,ribodepleted_rep1.tails_processed,"Ribodep_4hpf_rep1")
ribodep_hpf6_rep1.reshape <- reshape(ribodep_hpf6_rep1.data,ribodepleted_rep1.tails_processed,"Ribodep_6hpf_rep1")



ribodep_hpf2_rep2.reshape <- reshape(ribodep_hpf2_rep2.data,ribodepleted_rep2.tails_processed,"Ribodep_2hpf_rep2")
ribodep_hpf4_rep2.reshape <- reshape(ribodep_hpf4_rep2.data,ribodepleted_rep2.tails_processed,"Ribodep_4hpf_rep2")
ribodep_hpf6_rep2.reshape <- reshape(ribodep_hpf6_rep2.data,ribodepleted_rep2.tails_processed,"Ribodep_6hpf_rep2")



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
	category_sum <- aggregate(.~Gene_Type+Sample, data2[,c("Gene_Type","Sample",  "Gene_Type_Count_Norm")], sum)
  category_sum <- category_sum[order(-category_sum$Gene_Type_Count_Norm),]
  category_sum$Gene_Type <- factor(category_sum$Gene_Type, levels = unique(category_sum$Gene_Type))
	#Aggregate by Gene Type
	pdf(file=paste(label, "Gene_Type_Normalized_Count.pdf",sep="_"),height=6,width=20,onefile=FALSE)
	print(ggplot(category_sum, aes(fill=Sample, y=log(Gene_Type_Count_Norm+1), x=Gene_Type)) + 
    geom_bar(position="dodge", stat="identity")+
    theme_bw())
    dev.off()
}

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


tail_comparison_overall(ribodep_rep1_all, "Ribodepletion_All_Timepoints_rep1")
tail_comparison_overall(ribodep_rep2_all, "Ribodepletion_All_Timepoints_rep2")

tail_comparison_overall(ribodep_hpf2_replicability, "Ribodepletion_2hpf_Replicability")
tail_comparison_overall(ribodep_hpf4_replicability, "Ribodepletion_4hpf_Replicability")
tail_comparison_overall(ribodep_hpf6_replicability, "Ribodepletion_6hpf_Replicability")


##Overall Tail comparison (Single transcript) Both Rep1 and Rep2 
tail_comparison_overall_rep_merged<- function(data, label) {
	data2 <- subset(data,Gene_Type=="protein_coding")
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
	data3 <- subset(data2, Gene_Count > 20)
	pdf(file=paste(label, "Median_Tail_Per_Gene_Comparison_mRNA.pdf",sep="_"),height=6,width=10,onefile=FALSE)
	print(ggplot(data3, aes(x=Median_Length, color=Timepoint)) +
		geom_density(data=subset(data3, Rep=="Rep1"), linetype="dashed")+
		geom_density(data=subset(data3, Rep=="Rep2"))+
		theme_bw()+
  		coord_cartesian(xlim=c(-10, 175))+
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
	cor222<-paste("Pearson's p =",round(as.numeric(test$estimate),3))
	#pval<-paste("Pval =",test$p.value)
	mtext(paste(cor222))
	#mtext(paste(cor222,pval,sep=" ; ")) #Print the subtitle with the dataset correlation
	dev.off()
}

plot_denscols_with_corr_pearson("Rep1_Rep2_2hpf_genecount", log(data_merged_mRNA_min30_2hpf$Gene_Count.x) , log(data_merged_mRNA_min30_2hpf$Gene_Count.y), "2HPF_Rep1_logGeneCount", "2HPF_Rep2_logGeneCount" )
plot_denscols_with_corr_pearson("Rep1_Rep2_4hpf_genecount", log(data_merged_mRNA_min30_4hpf$Gene_Count.x) , log(data_merged_mRNA_min30_4hpf$Gene_Count.y), "4HPF_Rep1_logGeneCount", "4HPF_Rep2_logGeneCount" )

plot_denscols_with_corr_pearson("Rep1_Rep2_6hpf_genecount", log(data_merged_mRNA_min30_6hpf$Gene_Count.x) , log(data_merged_mRNA_min30_6hpf$Gene_Count.y), "6HPF_Rep1_logGeneCount", "6HPF_Rep2_logGeneCount" )





plot_denscols_with_corr_pearson("Rep1_Rep2_2hpf_tail", data_merged_mRNA_min30_2hpf$Median_Length.x , data_merged_mRNA_min30_2hpf$Median_Length.y, "2HPF_Rep1_Median_Length", "2HPF_Rep2_Median_Length" )


plot_denscols_with_corr_pearson("Rep1_Rep2_4hpf_tail", data_merged_mRNA_min30_4hpf$Median_Length.x , data_merged_mRNA_min30_4hpf$Median_Length.y, "4HPF_Rep1_Median_Length", "4HPF_Rep2_Median_Length" )


plot_denscols_with_corr_pearson("Rep1_Rep2_6hpf_tail", data_merged_mRNA_min30_6hpf$Median_Length.x , data_merged_mRNA_min30_6hpf$Median_Length.y, "6HPF_Rep1_Median_Length", "6HPF_Rep2_Median_Length" )






### CHECK THE REPLICABILITY OF THE SEQUINS

ribodep_rep1_rep2_all_sequin <- subset(ribodep_rep1_rep2_all, Gene_Type == "synthetic")
ribodep_rep1_rep2_all_sequin$Group <- gsub("_.*","",ribodep_rep1_rep2_all_sequin$Gene_Name)


      pdf(file="Sequin_Rep1_Rep2_Zebrafish_Tails.pdf",height=4,width=8,onefile=FALSE)
            print(ggplot(ribodep_rep1_rep2_all_sequin, aes(x = Group, y = tail_length, fill=Sample )) + 
              geom_boxplot(aes(fill = Sample),position = position_dodge(0.9), outlier.shape = NA) +
              theme_bw()+
              coord_cartesian(ylim = c(0, 100))+
              geom_hline(yintercept=30, linetype="dashed" ,size=0.3)+
              geom_hline(yintercept=60, linetype="dashed" ,size=0.3))
         dev.off()









ribodep_rep1_rep2_all_sequin_uniq <- ribodep_rep1_rep2_all_sequin[!duplicated(ribodep_rep1_rep2_all_sequin[c("Gene_Name", "Sample")]),]


#### Gene Type Count Comparison ###
simple_barplot_grouped <- function(data, label){
	pdf(file=paste(label, "Median_Tail_Length_Grouped.pdf",sep="_"),height=6,width=8,onefile=FALSE)
	print(ggplot(data, aes(fill=Sample, y=Median_Length, x=Group)) + 
    geom_bar(aes(Group, Median_Length, fill = Sample), position="dodge", stat = "summary", fun.y = "mean")+
    theme_bw())
    dev.off()
}


simple_barplot_grouped(ribodep_rep1_rep2_all_sequin_uniq, "Ribodepletion_BothReps_Sequins")



