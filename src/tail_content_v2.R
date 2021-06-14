### Tail Content v2 ##
library(ggplot2)
library(reshape2)
library(stringr)
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


filtered_3bases_rep1_alltimepoints <- rbind(filtered_3bases_2hpf_rep1,filtered_3bases_4hpf_rep1,filtered_3bases_6hpf_rep1)
filtered_3bases_rep2_alltimepoints <- rbind(filtered_3bases_2hpf_rep2,filtered_3bases_4hpf_rep2,filtered_3bases_6hpf_rep2)

boxplot <- function(data, label){
pdf(paste(label, "Tail_Length_Tail_Content_Boxplot_Per_Read.pdf",sep="_"),height=4,width=10,onefile=FALSE)
print(ggplot(data, aes(x=Sample, y=tail_length,colour=Group)) + 
  geom_boxplot())
dev.off()
}

boxplot(filtered_3bases_rep1_alltimepoints, "Rep1")
boxplot(filtered_3bases_rep2_alltimepoints, "Rep2")







