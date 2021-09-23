### Tail Content v2 ##
library(ggplot2)
library(reshape2)
library(stringr)
library(EnvStats)
library(ggpubr)

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



tail_content_2hpf.rep1 <- read.delim("cDNA786327_2hpf_tailcontent.csv")
tail_content_4hpf.rep1 <- read.delim("cDNA786327_4hpf_tailcontent.csv")
tail_content_6hpf.rep1 <- read.delim("cDNA786327_6hpf_tailcontent.csv")


tail_content_2hpf.rep2 <- read.delim("cDNA123791_2hpf_tailcontent.csv")
tail_content_4hpf.rep2 <- read.delim("cDNA123791_4hpf_tailcontent.csv")
tail_content_6hpf.rep2 <- read.delim("cDNA123791_6hpf_tailcontent.csv")


tail_content_2hpf <- rbind(tail_content_2hpf.rep1,tail_content_2hpf.rep2) 
tail_content_4hpf <- rbind(tail_content_4hpf.rep1,tail_content_4hpf.rep2) 
tail_content_6hpf <- rbind(tail_content_6hpf.rep1,tail_content_6hpf.rep2) 

# Reshape the tables and remove low quality reads
reshape<- function(data,tails,label, content) {
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
  ###CONTENT
  content2 <- subset(content, score > 40)
  content2 <- content2[-1,]
  content3 <- content2[!duplicated(content2[c("readid")]),]
  base_columns <- str_split_fixed(content3$potential_ttt, "", 20)
  base_freq_columns <-str_split_fixed(content3$polyT_ACGTN, ":", 5)
  colnames(base_columns) <- c("p20", "p19","p18","p17", "p16", "p15", "p14","p13","p12", "p11", "p10","p9","p8", "p7", "p6","p5","p4", "p3", "p2","p1")
  base_column_name <- colnames(base_columns)
  colnames(base_freq_columns) <- c("A","C", "G", "T", "N")
  base_freq_column_name <- colnames(base_freq_columns)
  content4 <- cbind(content3[,c("readid", "score", "most_common_base")], base_columns, base_freq_columns)
  content4$A <- as.numeric(as.character(content4$A))
  content4$C <- as.numeric(as.character(content4$C))
  content4$G <- as.numeric(as.character(content4$G))
  content4$T <- as.numeric(as.character(content4$T))
  content4$N <- as.numeric(as.character(content4$N))
  #Merge
  merged5 <- merge(merged4,content4,by.x="Read_ID", by.y="readid")
  merged5 <- subset(merged5 , Gene_Type=="protein_coding")
  merged6 <- merged5[,c("Read_ID", "Gene_Name", "tail_length", "Sample","most_common_base", base_column_name,  base_freq_column_name)]
  return(merged6)
}


ribodep_hpf2_rep1.reshape <- reshape(ribodep_hpf2_merged.data,ribodepleted_merged.tails_processed,"Ribodep_2hpf_rep1",tail_content_2hpf.rep1)
ribodep_hpf4_rep1.reshape <- reshape(ribodep_hpf4_merged.data,ribodepleted_merged.tails_processed,"Ribodep_4hpf_rep1",tail_content_4hpf.rep1)
ribodep_hpf6_rep1.reshape <- reshape(ribodep_hpf6_merged.data,ribodepleted_merged.tails_processed,"Ribodep_6hpf_rep1",tail_content_6hpf.rep1)

ribodep_hpf2_rep2.reshape <- reshape(ribodep_hpf2_merged.data,ribodepleted_merged.tails_processed,"Ribodep_2hpf_rep2",tail_content_2hpf.rep2)
ribodep_hpf4_rep2.reshape <- reshape(ribodep_hpf4_merged.data,ribodepleted_merged.tails_processed,"Ribodep_4hpf_rep2",tail_content_4hpf.rep2)
ribodep_hpf6_rep2.reshape <- reshape(ribodep_hpf6_merged.data,ribodepleted_merged.tails_processed,"Ribodep_6hpf_rep2",tail_content_6hpf.rep2)



ribodep_all_rep1 <- rbind(ribodep_hpf2_rep1.reshape, ribodep_hpf4_rep1.reshape, ribodep_hpf6_rep1.reshape)
ribodep_all_rep2 <- rbind(ribodep_hpf2_rep2.reshape, ribodep_hpf4_rep2.reshape, ribodep_hpf6_rep2.reshape)

ribodep_hpf2_both.reshape <- reshape(ribodep_hpf2_merged.data,ribodepleted_merged.tails_processed,"Ribodep_2hpf_both",tail_content_2hpf)
ribodep_hpf4_both.reshape <- reshape(ribodep_hpf4_merged.data,ribodepleted_merged.tails_processed,"Ribodep_4hpf_both",tail_content_4hpf)
ribodep_hpf6_both.reshape <- reshape(ribodep_hpf6_merged.data,ribodepleted_merged.tails_processed,"Ribodep_6hpf_both",tail_content_6hpf)


ribodep_all_both <- rbind(ribodep_hpf2_both.reshape,ribodep_hpf4_both.reshape,ribodep_hpf6_both.reshape)





#filter for the last 3 bases 
filter_last_3_bases <- function(data) {
	data_U <- subset(data, p1=="A" & p2=="A" & p3=="A" )
	data_U$Group <- "Last3U"
	#data_U_MostCommonA <- subset(data_U, most_common_base=="T")
	#data_U_MostCommonA$Group <- "Last3U_CommonbaseA"
	data_rest <- subset(data, p1!="A" |p2!="A" |p3!="A")
	data_rest$Group <- "Rest"
	data_all <- rbind(data_U,data_rest)
	return(data_all)
}


filtered_3bases_2hpf_rep1 <- filter_last_3_bases(ribodep_hpf2_rep1.reshape)
filtered_3bases_4hpf_rep1<- filter_last_3_bases(ribodep_hpf4_rep1.reshape)
filtered_3bases_6hpf_rep1 <- filter_last_3_bases(ribodep_hpf6_rep1.reshape)

filtered_3bases_2hpf_rep2 <- filter_last_3_bases(ribodep_hpf2_rep2.reshape)
filtered_3bases_4hpf_rep2<- filter_last_3_bases(ribodep_hpf4_rep2.reshape)
filtered_3bases_6hpf_rep2 <- filter_last_3_bases(ribodep_hpf6_rep2.reshape)


filtered_3bases_2hpf_both <- filter_last_3_bases(ribodep_hpf2_both.reshape)
filtered_3bases_4hpf_both<- filter_last_3_bases(ribodep_hpf4_both.reshape)
filtered_3bases_6hpf_both <- filter_last_3_bases(ribodep_hpf6_both.reshape)


filtered_3bases_rep1_alltimepoints <- rbind(filtered_3bases_2hpf_rep1,filtered_3bases_4hpf_rep1,filtered_3bases_6hpf_rep1)
filtered_3bases_rep2_alltimepoints <- rbind(filtered_3bases_2hpf_rep2,filtered_3bases_4hpf_rep2,filtered_3bases_6hpf_rep2)

filtered_3bases_both_alltimepoints <- rbind(filtered_3bases_2hpf_both,filtered_3bases_4hpf_both,filtered_3bases_6hpf_both)


boxplot <- function(data, label){
pdf(paste(label, "Tail_Length_Tail_Content_Boxplot_Per_Read.pdf",sep="_"),height=4,width=10,onefile=FALSE)
print(ggplot(data, aes(x=Sample, y=tail_length,colour=Group)) + 
  stat_n_text(aes(color = Group)) + 
  geom_boxplot())
dev.off()
}

boxplot(filtered_3bases_rep1_alltimepoints, "Rep1")
boxplot(filtered_3bases_rep2_alltimepoints, "Rep2")
boxplot(filtered_3bases_both_alltimepoints, "Both")



boxplot_stats <- function(data, label){
pdf(paste(label, "Tail_Length_Tail_Content_Boxplot_Per_Read.pdf",sep="_"),height=4,width=10,onefile=FALSE)
print(ggplot(data, aes(x=Sample, y=tail_length,colour=Group)) + 
  stat_n_text(aes(color = Group)) + 
  stat_compare_means(comparisons = my_comparisons, label.y = c(5, 5.5, 6))+
  geom_boxplot())
dev.off()
}

boxplot(filtered_3bases_rep1_alltimepoints, "Rep1")
boxplot(filtered_3bases_rep2_alltimepoints, "Rep2")
boxplot(filtered_3bases_both_alltimepoints, "Both")


my_comparisons <- list( c("Last3U","Rest"))

pdf(paste("Both_Rep_Tail_Length_Tail_Content_Boxplot_Per_Read.pdf",sep="_"),height=4,width=10,onefile=FALSE)
    print(ggplot(filtered_3bases_both_alltimepoints, aes(x = Group, y = tail_length)) + 
      geom_boxplot(aes(fill = Group),position = position_dodge(0.9)) +
      ylab("Tail Length (nt)")+
      stat_compare_means(comparisons = my_comparisons, label.y = c(310))+
      facet_wrap(~Sample,nrow=1)+
      stat_n_text() + 
      #scale_fill_manual(values=colors)+
      #coord_cartesian(ylim = c(0,175))+
      theme_bw())
  dev.off()



  maternal <- read.delim("264_top_maternal_decay.txt")
  maternal$Group <- "Maternal"



### check the reads 
check_reads <- function(data, label) {
  last3U <- subset(data, Group=="Last3U")
  rest <- subset(data, Group=="Rest")
  last3U_count <- as.data.frame(table(last3U$Gene_Name))
  colnames(last3U_count) <- c("Gene_Name", "Count_Last3U")
  rest_count <- as.data.frame(table(rest$Gene_Name))
  colnames(rest_count) <- c("Gene_Name", "Count_Rest")
  merged <- merge(last3U_count, rest_count, by.x="Gene_Name", by.y="Gene_Name")
  merged$count <- merged$Count_Last3U + merged$Count_Rest
  merged$Ratio_Last3U <- merged$Count_Last3U/ merged$count
  merged$Ratio_Rest<- merged$Count_Rest/ merged$count
  merged2 <- subset(merged,count>10)
  merged2$Sample <- label
  merged3 <- merge(merged2, maternal,all.x=TRUE, by.x="Gene_Name", by.y="Ensembl_Gene_ID")
  merged3$Group[is.na(merged3$Group)] <- "Rest"
  merged3$Ensembl_Transcript_ID <- NULL
  merged3$Percentage <-  merged3$Ratio_Last3U*100
  merged3<- merged3[!duplicated(merged3$Gene_Name), ]
  return(merged3)
}

read_count_2hpf <- check_reads(filtered_3bases_2hpf_both,"2HPF" )
read_count_4hpf <- check_reads(filtered_3bases_4hpf_both,"4HPF" )
read_count_6hpf <- check_reads(filtered_3bases_6hpf_both,"6HPF" )

read_count_all <- rbind(read_count_2hpf, read_count_4hpf, read_count_6hpf)


read_count_2hpf_rep1 <- check_reads(filtered_3bases_2hpf_rep1,"2HPF_rep1" )
read_count_2hpf_rep2 <- check_reads(filtered_3bases_2hpf_rep2,"2HPF_rep2" )



read_count_4hpf_rep1 <- check_reads(filtered_3bases_4hpf_rep1,"4HPF_rep1" )
read_count_4hpf_rep2 <- check_reads(filtered_3bases_4hpf_rep2,"4HPF_rep2" )


read_count_6hpf_rep1 <- check_reads(filtered_3bases_6hpf_rep1,"6HPF_rep1" )
read_count_6hpf_rep2 <- check_reads(filtered_3bases_6hpf_rep2,"6HPF_rep2" )

library(ggrepel)
scatter_plot <- function(data, label) {
  pdf(paste(label, "Tail_Content_Ratios_Min10Cov.pdf",sep="_"),height=5,width=5,onefile=FALSE)
  print(ggplot(data, aes(Ratio_Rest,Ratio_Last3U,colour = factor(Group)))+
  #geom_text_repel(data=subset(data, Ratio_Last3U>0.1),aes(y = Ratio_Last3U, label = Gene_Name), color = "black", size=6)+
  geom_point())
  dev.off()
}

scatter_plot(read_count_2hpf,"2HPF")
scatter_plot(read_count_4hpf,"4HPF")
scatter_plot(read_count_6hpf,"6HPF")


##Overall Tail comparison (Single transcript)
density_plot <- function(data, label) {
  pdf(file=paste(label, "Tail_Content_Ratios_Min10Cov_Density.pdf",sep="_"),height=3,width=5, onefile=FALSE)
  print(ggplot(data, aes(x=Ratio_Last3U, color=Group)) +
  geom_density()+
  theme_bw()+
  xlab("Last3U Ratio")+
  ylab("Density"))
  dev.off()
}


density_plot(read_count_2hpf,"2HPF")
density_plot(read_count_4hpf,"4HPF")
density_plot(read_count_6hpf,"6HPF")

density_plot(read_count_2hpf_rep1,"2HPF_rep1")
density_plot(read_count_2hpf_rep2,"2HPF_rep2")

density_plot(read_count_4hpf_rep1,"4HPF_rep1")
density_plot(read_count_4hpf_rep2,"4HPF_rep2")


density_plot(read_count_6hpf_rep1,"6HPF_rep1")
density_plot(read_count_6hpf_rep2,"6HPF_rep2")

  
#my_comparisons <- list( c("2HPF", "4HPF"), c("4HPF", "6HPF"), c("2HPF", "6HPF") )
boxplot <- function(data,label) {
  pdf(file=paste(label, "Tail_Content_Ratios_Min10Cov_boxplot.pdf",sep="_"),height=3,width=10, onefile=FALSE)
    print(ggplot(data, aes(x = Group, y = Percentage+0.1 )) + 
      geom_boxplot(aes(fill = Group),position = position_dodge(0.9)) +
      ylab("Last3U Percentage")+
      stat_compare_means()+
      facet_wrap(~Sample,nrow=1)+
      #facet_zoom(ylim = c(0, 0.1))+
      stat_n_text() + 
      scale_y_continuous(trans = 'log2')+
      #scale_fill_manual(values=colors)+
      #coord_cartesian(ylim = c(0,175))+
      theme_bw())
  dev.off()
}


boxplot(read_count_all, "AllTimepoints")


library(ggbeeswarm)


dotplot <- function(data,label) {
pdf(file=paste(label, "Tail_Content_Ratios_Min10Cov_dotplot.pdf",sep="_"),height=6,width=14,onefile=FALSE)
      print(ggplot(data, aes(x=Group, y=Percentage+0.1, fill=Group)) + 
        geom_quasirandom(varwidth = TRUE, aes())+
        geom_boxplot(aes(alpha=0), outlier.shape=NA)+
        stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                geom = "crossbar", width = 0.7, color="#c06c84")+
        stat_n_text()+
        theme_bw()+
        ggtitle(label)+
        facet_wrap(~Sample,nrow=1)+
        xlab("Group")+
        scale_y_continuous(trans = 'log2')+
        ylab("log scaled (pU Percentage)") +
        theme(axis.text=element_text(size=14),strip.text = element_text(size=13),
                axis.title=element_text(size=17,face="bold"),
                legend.title = element_text(size = 20),
                legend.text = element_text(color = "black", size=15)))
    dev.off()
}
dotplot(read_count_all, "AllTimepoints")



  
#my_comparisons <- list( c("2HPF", "4HPF"), c("4HPF", "6HPF"), c("2HPF", "6HPF") )
violinplot <- function(data,label) {
  pdf(file=paste(label, "Tail_Content_Ratios_Min10Cov_violinplot.pdf",sep="_"),height=3,width=10, onefile=FALSE)
    print(ggplot(data, aes(x = Group, y = Percentage+0.1 )) + 
      geom_violin(aes(fill = Group),position = position_dodge(0.9)) +
      ylab("Last3U Percentage")+
      stat_compare_means()+
      facet_wrap(~Sample,nrow=1)+
      #facet_zoom(ylim = c(0, 0.1))+
      stat_n_text() + 
      scale_y_continuous(trans = 'log2')+
      #scale_fill_manual(values=colors)+
      #coord_cartesian(ylim = c(0,175))+
      theme_bw())
  dev.off()
}


violinplot(read_count_all, "AllTimepoints")






##Overall Tail comparison (Single transcript)
density_plot_facet <- function(data, label) {
  pdf(file=paste(label, "Tail_Content_Ratios_Min10Cov_Density.pdf",sep="_"),height=3,width=7, onefile=FALSE)
  print(ggplot(data, aes(x=Percentage+0.1, color=Group)) +
  geom_density()+
  facet_wrap(~Sample,nrow=1)+
  scale_x_continuous(trans = 'log2')+
  theme_bw()+
  xlab("log scaled (pU Percentage)")+
  ylab("Density"))
  dev.off()
}
density_plot_facet(read_count_all, "AllTimepoints")















#Rep1 vs Rep2


rep1_vs_rep2_2hpf <- merge(read_count_2hpf_rep1[,c("Gene_Name","Ratio_Last3U","Ratio_Rest")],read_count_2hpf_rep2[,c("Gene_Name","Ratio_Last3U","Ratio_Rest")],by.x="Gene_Name", by.y="Gene_Name")



scatter_plot_rep1_vs_rep2 <- function(data, label) {
  pdf(paste(label, "Tail_Content_Ratios_Min10Cov_Rep1vsRep2.pdf",sep="_"),height=5,width=5,onefile=FALSE)
  print(ggplot(data, aes(Ratio_Last3U.x,Ratio_Last3U.y))+
  geom_point())
  dev.off()
}

scatter_plot_rep1_vs_rep2(rep1_vs_rep2_2hpf,"2HPF_Rep1_vs_Rep2")


#VENN DIAGRAM

read_count_2hpf_rep1_min20 <- subset(read_count_2hpf_rep1,count >20 )
read_count_2hpf_rep2_min20 <- subset(read_count_2hpf_rep2,count >20 )

read_count_4hpf_rep1_min20 <- subset(read_count_4hpf_rep1,count >20 )
read_count_4hpf_rep2_min20 <- subset(read_count_4hpf_rep2,count >20 )

read_count_6hpf_rep1_min20 <- subset(read_count_6hpf_rep1,count >20 )
read_count_6hpf_rep2_min20 <- subset(read_count_6hpf_rep2,count >20 )


read_count_2hpf_rep1_min20_2percent <- subset(read_count_2hpf_rep1_min20, Ratio_Last3U > 0.02)
read_count_2hpf_rep2_min20_2percent <- subset(read_count_2hpf_rep2_min20, Ratio_Last3U > 0.02)

read_count_4hpf_rep1_min20_2percent <- subset(read_count_4hpf_rep1_min20, Ratio_Last3U > 0.02)
read_count_4hpf_rep2_min20_2percent <- subset(read_count_4hpf_rep2_min20, Ratio_Last3U > 0.02)

read_count_6hpf_rep1_min20_2percent <- subset(read_count_6hpf_rep1_min20, Ratio_Last3U > 0.02)
read_count_6hpf_rep2_min20_2percent <- subset(read_count_6hpf_rep2_min20, Ratio_Last3U > 0.02)

library(VennDiagram)

 list_2hpf_rep1_rep2 <- list(read_count_2hpf_rep1_min20_2percent[,1],read_count_2hpf_rep2_min20_2percent[,1] )
 list_4hpf_rep1_rep2 <- list(read_count_4hpf_rep1_min20_2percent[,1],read_count_4hpf_rep2_min20_2percent[,1] )
 list_6hpf_rep1_rep2 <- list(read_count_6hpf_rep1_min20_2percent[,1],read_count_6hpf_rep2_min20_2percent[,1] )




venn_reps <- function(data, label){
venn.diagram(x =data,
     category.names = c("Rep1" , "Rep2"),
     filename = paste(label,'Rep1_vs_Rep2_Last3U_Genes.tiff',sep="_"),
     lwd = 2,
     lty = 'blank',
     fill = c("red", "blue"),
     output=TRUE,
       # Output features
      imagetype="tiff" ,
      height = 500 , 
      width = 500 , 
      resolution = 500,
      compression = "lzw",
     # Numbers
      cex = .3,
      fontface = "bold",
      fontfamily = "sans",
      # Set names
      cat.cex = 0.2,
      cat.fontface = "bold",
      cat.default.pos = "outer",
      #cat.pos = c(-27, 27),
      #cat.dist = c(0.055, 0.055),
      #cat.fontfamily = "sans",
      #rotation = 1
  )
}

venn_reps(list_2hpf_rep1_rep2, "2hpf")
venn_reps(list_4hpf_rep1_rep2, "4hpf")
venn_reps(list_6hpf_rep1_rep2, "6hpf")








# 2 ,4 and 6 hpf gene group overlap

read_count_2hpf_min20 <- subset(read_count_2hpf,count >20 )
read_count_4hpf_min20 <- subset(read_count_4hpf,count >20 )
read_count_6hpf_min20 <- subset(read_count_6hpf,count >20 )


read_count_2hpf_min20_2percent <- subset(read_count_2hpf_min20, Ratio_Last3U > 0.02)
read_count_4hpf_min20_2percent <- subset(read_count_4hpf_min20, Ratio_Last3U > 0.02)
read_count_6hpf_min20_2percent <- subset(read_count_6hpf_min20, Ratio_Last3U > 0.02)

library(VennDiagram)

 List_all <- list(read_count_2hpf_min20_2percent[,1],read_count_4hpf_min20_2percent[,1],read_count_6hpf_min20_2percent[,1] )


significant_pU_genes <- as.data.frame(Reduce(union, List_all))
colnames(significant_pU_genes) <- "Gene_Name"



 Reduce(intersect, List_all)



venn.diagram(x =List_all,
     category.names = c("2HPF" , "4HPF", "6HPF"),
     filename = 'Venn_All_Timepoints_PolyU_2percent_min20cov.tiff',
     lwd = 2,
     lty = 'blank',
     fill = c("#f5a25d", "#bedbbb", "#776d8a"),
     output=TRUE,
     # Output features
    imagetype="tiff" ,
    height = 6000 , 
    width = 6000 , 
    resolution = 6000,
    compression = "lzw",
   # Numbers
    cex = .5,
    fontface = "bold",
    fontfamily = "sans",
    # Set names
    cat.cex = 0.3,
    cat.fontface = "bold",
    cat.default.pos = "outer",
    cat.pos = c(-25, 0, 25 ),
    cat.dist = c(0.15, 0.15 ,-0.5),
    cat.fontfamily = "sans",
    #rotation = 1
)




##### LINE PLOT FOR THESE GENES

rep1_vs_rep2_2hpf <- merge(read_count_2hpf_rep1[,c("Gene_Name","Percentage")],read_count_2hpf_rep2[,c("Gene_Name","Percentage")],by.x="Gene_Name", by.y="Gene_Name")
rep1_vs_rep2_4hpf <- merge(read_count_4hpf_rep1[,c("Gene_Name","Percentage")],read_count_4hpf_rep2[,c("Gene_Name","Percentage")],by.x="Gene_Name", by.y="Gene_Name")
rep1_vs_rep2_6hpf <- merge(read_count_6hpf_rep1[,c("Gene_Name","Percentage")],read_count_6hpf_rep2[,c("Gene_Name","Percentage")],by.x="Gene_Name", by.y="Gene_Name")



rep1_vs_rep2_2hpf$Mean_Percentage <- (rep1_vs_rep2_2hpf$Percentage.x + rep1_vs_rep2_2hpf$Percentage.y)/2
rep1_vs_rep2_4hpf$Mean_Percentage <- (rep1_vs_rep2_4hpf$Percentage.x + rep1_vs_rep2_4hpf$Percentage.y)/2
rep1_vs_rep2_6hpf$Mean_Percentage <- (rep1_vs_rep2_6hpf$Percentage.x + rep1_vs_rep2_6hpf$Percentage.y)/2


significant_pU_genes_2hpf <- merge(significant_pU_genes, rep1_vs_rep2_2hpf, by.x="Gene_Name", by.y="Gene_Name")
significant_pU_genes_24hpf <- merge(significant_pU_genes_2hpf, rep1_vs_rep2_4hpf, by.x="Gene_Name", by.y="Gene_Name")
significant_pU_genes_246hpf <- merge(significant_pU_genes_24hpf, rep1_vs_rep2_6hpf, by.x="Gene_Name", by.y="Gene_Name")

significant_pU_genes_246hpf_final <- significant_pU_genes_246hpf[,c("Gene_Name", "Mean_Percentage.x", "Mean_Percentage.y", "Mean_Percentage")]

colnames(significant_pU_genes_246hpf_final) <- c("Gene_Name", "Percentage_2hpf", "Percentage_4hpf", "Percentage_6hpf")


  maternal <- read.delim("264_top_maternal_decay.txt")
  maternal$Group <- "Maternal"

  significant_pU_genes_246hpf_maternal <- merge(significant_pU_genes_246hpf_final, maternal,all.x=TRUE, by.x="Gene_Name", by.y="Ensembl_Gene_ID")
  significant_pU_genes_246hpf_maternal$Group[is.na(significant_pU_genes_246hpf_maternal$Group)] <- "Rest"
  significant_pU_genes_246hpf_maternal$Ensembl_Transcript_ID <- NULL

significant_pU_genes_246hpf_maternal_melted <- melt(significant_pU_genes_246hpf_maternal)



line_plot_pU <- function(data, label) {
pdf(file=paste(label,"pU_percentage_significant_genes.pdf",sep="_"),height=3,width=6,onefile=FALSE)
  print(ggplot(data=data , aes(x=variable, y=value, group=Gene_Name, color=Group)) +
    geom_line()+
   # geom_line(data=subset(data, Replicate=="Rep2"))+
    theme_bw()+
   # scale_y_continuous(trans = 'log2')+
    ylab("pU Percentage"))
    dev.off()
}


line_plot_pU(significant_pU_genes_246hpf_maternal_melted, "Mean_all_time")



  print(ggplot(data=data , aes(x=variable, y=count_percentage, group=value, color=value)) +





test <- subset(subs_last3U, Read_ID == "88012a53-f012-4318-90a7-cdd6213066fe")


88012a53-f012-4318-90a7-cdd6213066fe


#### EXTRACT READS FOR THE LAST3U 

extract_reads_last3U <- function(data,gene,label) {
  subs <- subset(data, Gene_Name==gene)
  subs_last3U <-  subset(data, Group=="Last3U")
  subs_toextract <- subs_last3U[,c("Read_ID")]
  write.table(subs_toextract, file=paste(gene,label,"reads.tsv", sep="_"), quote=FALSE, col.names=FALSE, row.names=FALSE)
}


extract_reads_last3U(filtered_3bases_2hpf_both,"ENSDARG00000006580","2HPF")




java -jar /users/enovoa/boguzhan/Software/picard/build/libs/picard.jar FilterSamReads \
       I=2hpf_nonrRNA_both_trimmed_porechop.genome11_sequin.sorted.bam \
       O=ENSDARG00000006580_2HPF.bam\
       READ_LIST_FILE=ENSDARG00000006580_2HPF_reads.tsv\
       FILTER=includeReadList

samtools index ENSDARG00000006580_2HPF.bam






java -jar /users/enovoa/boguzhan/Software/picard/build/libs/picard.jar FilterSamReads \
       I=2hpf.genome11_sequin.restRNAs.bam \
       O=ENSDARG00000006580_2HPF_untrimmed.bam\
       READ_LIST_FILE=ENSDARG00000006580_2HPF_reads.tsv\
       FILTER=includeReadList

samtools index ENSDARG00000006580_2HPF_untrimmed.bam

