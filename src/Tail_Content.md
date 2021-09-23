# Tail Content
```bash
python /users/enovoa/hliu/analysis/mygit/projects_notes/oz_ribo_switch_find_adapter_in_reads/adapterFinder.py tails.csv cDNA123791_2hpf.fastq 0.7 > cDNA123791_2hpf_tailcontent.csv

python /users/enovoa/hliu/analysis/mygit/projects_notes/oz_ribo_switch_find_adapter_in_reads/adapterFinder.py tails.csv cDNA123791_4hpf.fastq 0.7 > cDNA123791_4hpf_tailcontent.csv

python /users/enovoa/hliu/analysis/mygit/projects_notes/oz_ribo_switch_find_adapter_in_reads/adapterFinder.py tails.csv cDNA123791_6hpf.fastq 0.7 > cDNA123791_6hpf_tailcontent.csv

python /users/enovoa/hliu/analysis/mygit/projects_notes/oz_ribo_switch_find_adapter_in_reads/adapterFinder.py cDNA786327_tails.csv cDNA786327_2hpf.fastq 0.7 > cDNA786327_2hpf_tailcontent.csv

python /users/enovoa/hliu/analysis/mygit/projects_notes/oz_ribo_switch_find_adapter_in_reads/adapterFinder.py cDNA786327_tails.csv cDNA786327_4hpf.fastq 0.7 > cDNA786327_4hpf_tailcontent.csv

python /users/enovoa/hliu/analysis/mygit/projects_notes/oz_ribo_switch_find_adapter_in_reads/adapterFinder.py cDNA786327_tails.csv cDNA786327_6hpf.fastq 0.7 > cDNA786327_6hpf_tailcontent.csv


qsub -cwd -q long-sl7 -l virtual_free=40G tail.sh





sed '$d'  cDNA786327_6hpf_tailcontent.csv >  cDNA786327_6hpf_tailcontent2.csv
sed '$d'  cDNA123791_6hpf_tailcontent.csv >  cDNA123791_6hpf_tailcontent2.csv

```


```R
library(ggplot2)
library(reshape2)
library(stringr)



rep2.2hpf <- read.delim("cDNA123791_2hpf_tailcontent.csv")
rep2.4hpf <- read.delim("cDNA123791_4hpf_tailcontent.csv")
rep2.6hpf <- read.delim("cDNA123791_6hpf_tailcontent2.csv")


rep1.2hpf <- read.delim("cDNA786327_2hpf_tailcontent.csv")
rep1.4hpf <- read.delim("cDNA786327_4hpf_tailcontent.csv")
rep1.6hpf <- read.delim("cDNA786327_6hpf_tailcontent2.csv")


manipulate <- function(data) {
data<- subset(data, score > 40)
data <- data[-1,]
data2 <- data[!duplicated(data[c("readid")]),]
base_columns <- str_split_fixed(data2$potential_ttt, "", 20)
base_freq_columns <-str_split_fixed(data2$polyT_ACGTN, ":", 5)
colnames(base_columns) <- c("p20", "p19","p18","p17", "p16", "p15", "p14","p13","p12", "p11", "p10","p9","p8", "p7", "p6","p5","p4", "p3", "p2","p1")
colnames(base_freq_columns) <- c("A","C", "G", "T", "N")
data3 <- cbind(data2[,c("readid", "score", "most_common_base")], base_columns, base_freq_columns)
data3$A <- as.numeric(as.character(data3$A))
data3$C <- as.numeric(as.character(data3$C))
data3$G <- as.numeric(as.character(data3$G))
data3$T <- as.numeric(as.character(data3$T))
data3$N <- as.numeric(as.character(data3$N))
return(data3)
}

processed.rep2.2hpf <- manipulate(rep2.2hpf)
processed.rep1.2hpf <- manipulate(rep1.2hpf)

processed.rep2.4hpf <- manipulate(rep2.4hpf)
processed.rep1.4hpf <- manipulate(rep1.4hpf)

processed.rep2.6hpf <- manipulate(rep2.6hpf)
processed.rep1.6hpf <- manipulate(rep1.6hpf)



## First making a barplot for the percentage of the bases in 20 nt
summing <- function(data) {
A_total <- sum(data$A)
G_total <- sum(data$G)
C_total <- sum(data$C)
T_total <- sum(data$T)
Total <- A_total+G_total+C_total+T_total
T_freq <- T_total/Total
C_freq <- C_total/Total
G_freq <- G_total/Total
A_freq <- A_total/Total
Frequencies <- as.data.frame(rbind(A_freq,G_freq,C_freq,T_freq))
Frequencies$Base <- rownames(Frequencies)
Frequencies$V1 <- as.numeric(as.character(Frequencies$V1))
Frequencies$V1 <- Frequencies$V1*100
return(Frequencies)
}

frequency.rep2.2hpf <- summing(processed.rep2.2hpf)
frequency.rep1.2hpf <- summing(processed.rep1.2hpf)

frequency.rep2.4hpf <- summing(processed.rep2.4hpf)
frequency.rep1.4hpf <- summing(processed.rep1.4hpf)

frequency.rep2.6hpf <- summing(processed.rep2.6hpf)
frequency.rep1.6hpf <- summing(processed.rep1.6hpf)




plot_percentage <- function(data, label){
pdf(file=paste(label,"Tail_Content_20nt_base_percentage.pdf",sep="_"),height=6,width=10,onefile=FALSE)
	print(ggplot(data=data, aes(x=Base, y=V1, fill=Base)) +
		theme_bw()+
		ggtitle(label)+
		scale_y_continuous(trans = 'log2')+
		ylab("Log Scaled Percentage")+
  		geom_bar(stat="identity"))
	dev.off()
}
plot_percentage(frequency.rep2.2hpf, "Rep2_2HPF")
plot_percentage(frequency.rep1.2hpf, "Rep1_2HPF")


plot_percentage(frequency.rep2.4hpf, "Rep2_4HPF")
plot_percentage(frequency.rep1.4hpf, "Rep1_4HPF")


plot_percentage(frequency.rep2.6hpf, "Rep2_6HPF")
plot_percentage(frequency.rep1.6hpf, "Rep1_6HPF")


## Looking at each position
per_position_15 <- function(data, label){
data2 <- data[,c( "readid", "score", "most_common_base", "p15", "p14","p13","p12", "p11", "p10","p9","p8", "p7", "p6","p5","p4", "p3", "p2","p1")]
data3 <- melt(data2, id.vars=c("readid", "score", "most_common_base"))
data4 <- subset(data3, value =="A" | value =="C" |value =="G" |value =="T")
data4$count <- 1
position_sum <- aggregate(.~variable+value, data4[,c("variable", "value","count")], sum)
final <- vector()
for (pos in unique(position_sum$variable)) {
	subs <- subset(position_sum, variable==pos)
	total <- sum(subs$count)
	subs$count_percentage <- subs$count/total
	final <- rbind(final,subs)
}
final$Replicate <- label
return(final)
}


per_pos_15.rep2.2hpf <- per_position_15(processed.rep2.2hpf, "Rep2")
per_pos_15.rep1.2hpf <- per_position_15(processed.rep1.2hpf, "Rep1")

per_pos_15.rep2.4hpf <- per_position_15(processed.rep2.4hpf, "Rep2")
per_pos_15.rep1.4hpf <- per_position_15(processed.rep1.4hpf, "Rep1")

per_pos_15.rep2.6hpf <- per_position_15(processed.rep2.6hpf, "Rep2")
per_pos_15.rep1.6hpf <- per_position_15(processed.rep1.6hpf, "Rep1")




per_pos_15.2hpf <- rbind(per_pos_15.rep1.2hpf, per_pos_15.rep2.2hpf)
per_pos_15.4hpf <- rbind(per_pos_15.rep1.4hpf, per_pos_15.rep2.4hpf)
per_pos_15.6hpf <- rbind(per_pos_15.rep1.6hpf, per_pos_15.rep2.6hpf)


per_pos_15.2hpf$Timepoint <- "2HPF"
per_pos_15.4hpf$Timepoint <- "4HPF"
per_pos_15.6hpf$Timepoint <- "6HPF"


per_pos_15.all <- rbind(per_pos_15.2hpf,per_pos_15.4hpf,per_pos_15.6hpf)



plot_per_pos <- function(data, label) {
pdf(file=paste(label,"Tail_Content_15nt_Distribution_LinePlot.pdf",sep="_"),height=3,width=6,onefile=FALSE)
	print(ggplot(data=data , aes(x=variable, y=count_percentage, group=value, color=value)) +
    geom_line()+
    theme_bw()+
    scale_y_continuous(trans = 'log2')+
    ylab("Percentage"))
    dev.off()
}


plot_per_pos(per_pos_15.rep2.2hpf , "Rep2_2HPF")
plot_per_pos(per_pos_15.rep1.2hpf , "Rep1_2HPF")




plot_per_pos_both <- function(data, label) {
pdf(file=paste(label,"Tail_Content_15nt_Distribution_LinePlot.pdf",sep="_"),height=3,width=6,onefile=FALSE)
	print(ggplot(data=data , aes(x=variable, y=count_percentage, group=value, color=value)) +
    geom_line(data=subset(data, Replicate=="Rep1"), linetype="dashed")+
    geom_line(data=subset(data, Replicate=="Rep2"))+
    theme_bw()+
    scale_y_continuous(trans = 'log2')+
    ylab("Percentage"))
    dev.off()
}


plot_per_pos_both(per_pos_15.2hpf , "2HPF")
plot_per_pos_both(per_pos_15.4hpf , "4HPF")
plot_per_pos_both(per_pos_15.6hpf , "6HPF")











merge_bed <- function(data, bedfile, label) {
bed_file <- read.delim(bedfile,header=FALSE)
bed_file2 <- bed_file[,c("V1", "V4", "V5", "V6", "V16", "V17")]
colnames(bed_file2) <- c("Chr","Read_ID", "Quality", "Strand", "Gene_Name", "Gene_Type")
bed_file3 <- subset(bed_file2, Quality > 30)
merged_bed <- merge(data, bed_file3, by.y="Read_ID", by.x="readid")
final2<- vector()
for (gene in unique(merged_bed$Gene_Name)){
	subs <- subset(merged_bed ,Gene_Name==gene)
		gene_count <- nrow(subs)
		last_three_row_a <- nrow(subset(subs,  p1 =="A" &  p2 =="A" &  p3 =="A"))
		common_base_a <- nrow(subset(subs,  most_common_base=="A"))
		common_base_t <- nrow(subset(subs,  most_common_base=="T"))
		Gene_Name <- gene
		final <- as.data.frame(cbind(gene,gene_count,last_three_row_a,common_base_a,common_base_t))
		final$gene_count <- as.numeric(as.character(final$gene_count ))
		final$common_base_a <- as.numeric(as.character(final$common_base_a ))
		final$common_base_t <- as.numeric(as.character(final$common_base_t ))
		final$last_three_row_a <- as.numeric(as.character(final$last_three_row_a ))
		final$last_three_row_a_perc <- final$last_three_row_a/final$gene_count
		final$common_base_a_perc <- final$common_base_a/final$gene_count
		final$common_base_t_perc <- final$common_base_t/final$gene_count
		final2 <- rbind(final2,final)
	}
final2$Sample <- label
return(final2)
}




mergedbed.rep2.2hpf <- merge_bed(processed.rep2.2hpf,2hpf.genome11_sequin_ALLRNAs_Merged_Rep2.bed,"2HPF_Rep2")




### Rep1 vs Rep2
scatter <- function(data,label) {
      pdf(file=paste(label,"Gene_Count_Common_Base_A.pdf",sep="_"),height=4,width=4,onefile=FALSE)
		print(ggplot(data, aes(x=gene_count, y=common_base_t)) + 
		theme_bw()+ 
        ggtitle(label))+
		xlab("Gene_Count")+
        ylab("Common Base A")+
  		geom_point())
  		dev.off()
}

scatter(mergedbed.rep2.2hpf, "2HPF_Rep2")










data_a123 <- subset(data4,  p1 =="A" &  p2 =="A" &  p3 =="A")




data_U_common <- subset(data3, most_common_base=="A")
bed_file <- read.delim("2hpf.genome11_sequin_ALLRNAs_Merged_Rep2.bed",header=FALSE)
bed_file2 <- bed_file[,c("V1", "V4", "V5", "V6", "V16", "V17")]
colnames(bed_file2) <- c("Chr","Read_ID", "Quality", "Strand", "Gene_Name", "Gene_Type")
bed_file3 <- subset(bed_file2, Quality > 30)
merged_bed <- merge(data4, bed_file3, by.y="Read_ID", by.x="readid")


merged_bed_a123 <- subset(merged_bed,  p1 =="A" &  p2 =="A" )











###################################
##### ANOTHER QUESTION TO ANSWER
###################################





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







merge_bed <- function(data, bedfile, label) {
merged_bed <- merge(data, bedfile, by.y="Read_ID", by.x="readid")
final2<- vector()
for (gene in unique(merged_bed$Gene_Name)){
	subs <- subset(merged_bed ,Gene_Name==gene)
		gene_count <- nrow(subs)
		last_three_row_a <- nrow(subset(subs,  p1 =="A" &  p2 =="A" &  p3 =="A"))
		common_base_a <- nrow(subset(subs,  most_common_base=="A"))
		common_base_t <- nrow(subset(subs,  most_common_base=="T"))
		Gene_Name <- gene
		final <- as.data.frame(cbind(gene,gene_count,last_three_row_a,common_base_a,common_base_t))
		final$gene_count <- as.numeric(as.character(final$gene_count ))
		final$common_base_a <- as.numeric(as.character(final$common_base_a ))
		final$common_base_t <- as.numeric(as.character(final$common_base_t ))
		final$last_three_row_a <- as.numeric(as.character(final$last_three_row_a ))
		final$last_three_row_a_perc <- final$last_three_row_a/final$gene_count
		final$common_base_a_perc <- final$common_base_a/final$gene_count
		final$common_base_t_perc <- final$common_base_t/final$gene_count
		final2 <- rbind(final2,final)
	}
final2$Sample <- label
return(final2)
}

mergedbed.rep2.2hpf <- merge_bed(processed.rep2.2hpf,ribodep_all_merged,"2HPF_Rep2")

mergedbed.rep2.2hpf_2 <- merge(mergedbed.rep2.2hpf ,ribodep_all_merged[,c("Gene_Name", "Median_Length")], by.x="gene", by.y="Gene_Name" )
mergedbed.rep2.2hpf_2  <- subset(mergedbed.rep2.2hpf_2 , )


mergedbed.rep2.2hpf_unique <-mergedbed.rep2.2hpf_2 [!duplicated(mergedbed.rep2.2hpf_2 [c("gene")]),]

mergedbed.rep2.2hpf_unique_top <- subset(mergedbed.rep2.2hpf_unique, last_three_row_a_perc >0.5)
mergedbed.rep2.2hpf_unique_top$Group <- "Last_Three_U"
mergedbed.rep2.2hpf_unique_rest <- subset(mergedbed.rep2.2hpf_unique, last_three_row_a_perc <0.5)
mergedbed.rep2.2hpf_unique_rest$Group <- "Rest"


content_groups_both <- rbind(mergedbed.rep2.2hpf_unique_top,mergedbed.rep2.2hpf_unique_rest)



pdf("Tail_Lenght_Tail_Content_Boxplot.pdf",height=4,width=4,onefile=FALSE)
print(ggplot(content_groups_both, aes(x=Group, y=Median_Length)) + 
  geom_boxplot())
dev.off()



