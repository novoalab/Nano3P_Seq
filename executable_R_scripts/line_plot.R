#How to execute
#Rscript --vanilla line_plot.R tails_file.csv 2hpf.bed 4hpf.bed 6hpf.bed gene_list
#Library needed
library(ggplot2)
library(ggbeeswarm)
library(EnvStats)
library(ggpubr)

# Reading arguments from command line
args = commandArgs(trailingOnly=TRUE)

#Arguments
tail <- read.delim(args[1],sep=",")  #1st variable


bed_2hpf <- read.delim(args[2], header=FALSE) # 2nd variable
bed_4hpf <- read.delim(args[3], header=FALSE) # 3rd variable
bed_6hpf <- read.delim(args[4], header=FALSE) # 4th variable



gene_list <- read.delim(args[5]) #5th variable



#Process the Tail
manipulate_tail<- function(data) { 
  #Convert some columns into numeric
  data$tail_length <- as.numeric(as.character(data$tail_length))
  data$tail_start <- as.numeric(as.character(data$tail_start))
  data$tail_end <- as.numeric(as.character(data$tail_end))
  data$samples_per_nt <- as.numeric(as.character(data$samples_per_nt))
  #Correction of the poly(A) tail estimation
  data$samples_per_nt <- data$samples_per_nt*1.2
  data$tail_length <- (data$tail_end-data$tail_start)/data$samples_per_nt
  #Convert NA values into 0
  data[["tail_length"]][is.na(data[["tail_length"]])] <- 0
  #keep the rows with valid tail
  data_filt <- subset(data, tail_is_valid=="TRUE")
  #Keep the rows with polyT value
  data_filt_T <- subset(data_filt, read_type=="polyT")
  return(data_filt_T)
}

tail.processed <- manipulate_tail(tail)



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
  merged3 <- merge(merged2, gene_type, by.x="Gene_Type", by.y="Gene_Type")
  merged4 <-  merged3[!duplicated(merged3[c("Gene_Name")]),]
  merged5 <- merged4[,c("Gene_Type", "Gene_Name", "Median_Length", "Gene_Count_Norm")]
  merged6 <- subset(merged5 , Gene_Type=="protein_coding")
  return(merged6)
}

processed.2hpf <- reshape(bed_2hpf,tail.processed,"2hpf")
processed.4hpf <- reshape(bed_4hpf,tail.processed,"4hpf")
processed.6hpf <- reshape(bed_6hpf,tail.processed,"6hpf")



columns <-c("Gene_Type", "Gene_Name")

#Merge 2 and 4 hpf
ribodep_all_unique_2hpf_4hpf <- merge(ribodep_all_unique_2hpf, ribodep_all_unique_4hpf, all.x=TRUE, all.y=TRUE, by.x=columns, by.y=columns, suffixes = c(".2hpf", ".4hpf"))
#Merge 2/4 and 6 hpf
ribodep_all_unique_246hpf <- merge(ribodep_all_unique_2hpf_4hpf, ribodep_all_unique_6hpf, all.x=TRUE, all.y=TRUE, by.x=columns, by.y=columns)
#Rename the column names
colnames(ribodep_all_unique_246hpf) <- c(colnames(ribodep_all_unique_2hpf_4hpf),"Median_Length.6hpf", "Gene_Count_Norm.6hpf")
#Seperate character and numeric columns
ribodep_all_unique_246hpf_characters <- ribodep_all_unique_246hpf[,c(1,2)]
ribodep_all_unique_246hpf_numbers <- ribodep_all_unique_246hpf[,-c(1,2)]
#Convert NA columns into 0
ribodep_all_unique_246hpf_numbers[is.na(ribodep_all_unique_246hpf_numbers)] <- 0
#Rebind the columns
ribodep_all_unique_246hpf_final <- cbind(ribodep_all_unique_246hpf_characters,ribodep_all_unique_246hpf_numbers )


#Merge the data with gene list containing top maternal decay genes
ribodep_all_unique_246hpf_final2 <- merge(ribodep_all_unique_246hpf_final, gene_list,all.x=TRUE, by.x="Gene_Name", by.y="Ensembl_Gene_ID")
#Assign the rows containing NA into "Rest" group
ribodep_all_unique_246hpf_final2$Group1[is.na(ribodep_all_unique_246hpf_final2$Group)] <- "Rest"
#Keep only the Gene Counts
ribodep_all_unique_246hpf_genecounts <- ribodep_all_unique_246hpf_final2[,c("Gene_Name","Gene_Type","Group", "Gene_Count_Norm.2hpf","Gene_Count_Norm.4hpf","Gene_Count_Norm.6hpf" )]


#Normalize by 2hpf
ribodep_all_unique_246hpf_genecounts$Gene_Count_Norm2.2hpf <- (ribodep_all_unique_246hpf_genecounts$Gene_Count_Norm.2hpf+1)/(ribodep_all_unique_246hpf_genecounts$Gene_Count_Norm.2hpf+1)
ribodep_all_unique_246hpf_genecounts$Gene_Count_Norm2.4hpf <- (ribodep_all_unique_246hpf_genecounts$Gene_Count_Norm.4hpf+1)/(ribodep_all_unique_246hpf_genecounts$Gene_Count_Norm.2hpf+1)
ribodep_all_unique_246hpf_genecounts$Gene_Count_Norm2.6hpf <- (ribodep_all_unique_246hpf_genecounts$Gene_Count_Norm.6hpf+1)/(ribodep_all_unique_246hpf_genecounts$Gene_Count_Norm.2hpf+1)

#KEep the normalized columns
ribodep_all_unique_246hpf_genecounts_relativeto2hpf <- ribodep_all_unique_246hpf_genecounts[,c("Gene_Name","Gene_Type", "Group1", "Gene_Count_Norm2.2hpf" ,"Gene_Count_Norm2.4hpf", "Gene_Count_Norm2.6hpf")]

#Melt the table
ribodep_all_unique_246hpf_genecounts_melted_relativeto2hpf <- melt(ribodep_all_unique_246hpf_genecounts_relativeto2hpf)


pdf(file="Zebrafish_Embryos_GeneCount_Line_Plot_Maternal_vs_Rest_Relative_to_2hpf.pdf",height=4,width=5,onefile=FALSE)
print(ggplot(ribodep_all_unique_246hpf_genecounts_melted_mrna_relativeto2hpf, aes(x=variable, y=log(value), group=Gene_Name)) +
  geom_line(data=subset(ribodep_all_unique_246hpf_genecounts_melted_mrna_relativeto2hpf, Group=="Rest"), alpha=1/5, color="grey")+
  geom_line(data=subset(ribodep_all_unique_246hpf_genecounts_melted_mrna_relativeto2hpf, Group=="Maternal"), alpha=1/2, color="red")+
  theme_bw()+
  geom_point(data=subset(ribodep_all_unique_246hpf_genecounts_melted_mrna_relativeto2hpf, Group=="Rest"), color="grey")+
    geom_point(data=subset(ribodep_all_unique_246hpf_genecounts_melted_mrna_relativeto2hpf, Group=="Maternal"), color="red"))

dev.off()


