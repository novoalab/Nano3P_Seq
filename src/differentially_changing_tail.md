########################################################
######## Differentially changing tail analysis #########
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




 

#### INVESTIGATE THE MATERNAL RNAs GENE COUNTS
ribodep_all_unique <-ribodep_all_merged[!duplicated(ribodep_all_merged[c("Gene_Name", "Sample")]),]
ribodep_all_unique <- subset(ribodep_all_unique, Gene_Count>30)


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


ribodep_all_unique_246hpf_final_mRNA <- subset(ribodep_all_unique_246hpf_final, Gene_Type=="protein_coding")

scatter_data <- ribodep_all_unique_246hpf_final_mRNA[,c("Gene_Name", "Median_Length.2hpf", "Median_Length.4hpf", "Median_Length.6hpf")]


scatter_data_melted <- melt(scatter_data)
#scatter_data_melted$value <- log(scatter_data_melted$value+1)



#Calculate THRESHOLD for the specificity
genemean<- aggregate(scatter_data_melted[, 3], list(scatter_data_melted$Gene_Name), mean) #Row mean grouped by Gene
colnames(genemean)<- c("Gene_Name", "genemean") 
scatter2<- plyr::join(scatter_data_melted, genemean, by="Gene_Name") #Add Rowmeans to the original data
scatter2$abs<- scatter2$value- scatter2$genemean #absolute distance of gene's expresion in tissue A from mean expression of this gene ins all tissues
res_vec <- vector() #res_vec file is a vector file that will contain all the residuals (tissue.vs.all) for all of the genes in all tissues
for (tissue in unique(scatter2$variable)){ #for every single tissue
  subset <- scatter2[with(scatter2, scatter2$variable %in% tissue),] #extract the data for a specific tissue
  res<- rlm(subset$value ~0 + subset$genemean) #linear model for that tissue 
  res_vec= c(res$residuals,res_vec)#this contains residuals for every gene in every tissuevsall combination
}
threshold <- 2.5*sd(res_vec) #The threshold is 2.5 times the standard deviation of all the residuals


##Seperate plots and calculations for each tissue
specific_genes<-vector() 
for (tissue in unique(scatter2$variable)){ #for each tissue
subset <- scatter2[with(scatter2, scatter2$variable %in% tissue),] #extract the data for a specific tissue
res<- rlm(subset$value ~0 + subset$genemean)#linear model for that tissue 
subset$res<- res$residuals #add residual values to the matrix
subset$diff<- abs(subset$res)-threshold #difference between gene's residual and threshold
spec<-subset(subset,  abs>3) #extract specific genes in each tissue
specific_genes<- rbind(spec,specific_genes) #add these genes to the initial data
pdf(file=paste(tissue,"tail_difference.pdf",sep="."),height=5,width=5)
print(ggplot(subset, aes(x=genemean, y=value,label=Gene_Name)) + 
    #scale_x_continuous(limits= c(0,5))+
    #scale_y_continuous(limits= c(0,5))+
    geom_point(data=subset, col="black",size=0.5)+ #All data points will be black
    geom_point(data=subset(subset, abs>50),col="red",size=2)+ #Except the specific genes
    #geom_text_repel(data=subset(subset,  abs>50),segment.size  = 0.4,segment.color = "grey50",)+ #Add text to the specific genes
    geom_smooth(method=rlm, formula = y ~0 + x, size=0.5,fullrange=TRUE)+ #abline will be from rlm function that passes through 0,0
    xlab("Mean Tail Length All Timepoints")+
    ylab(paste("Median Tail Length",tissue,sep=" "))+
    theme(panel.background = element_blank(),
         panel.border=element_rect(fill=NA),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         strip.background=element_blank(),
         axis.text.x=element_text(colour="black"),
         axis.text.y=element_text(colour="black"),
         axis.ticks=element_line(colour="black"),
        plot.margin=unit(c(1,1,1,1),"line")))
dev.off()
}
write.table(specific_genes, file="specific_genes_encode.tsv",quote=FALSE, row.names=FALSE,sep="\t")






scatter_data_melted_2 <- subset(scatter_data_melted, variable=="Median_Length.2hpf")
outlier <- quantile(scatter_data_melted_2$value, 0.95)
scatter_data_melted_2_tops <- subset(scatter_data_melted_2, value > outlier)



scatter_data_melted_4 <- subset(scatter_data_melted, variable=="Median_Length.4hpf")
outlier <- quantile(scatter_data_melted_4$value, 0.95)
scatter_data_melted_4_tops <- subset(scatter_data_melted_4, value > outlier)


scatter_data_melted_6 <- subset(scatter_data_melted, variable=="Median_Length.6hpf")
outlier <- quantile(scatter_data_melted_6$value, 0.95)
scatter_data_melted_6_tops <- subset(scatter_data_melted_6, value > outlier)



all_top_genes <- rbind(scatter_data_melted_2_tops,scatter_data_melted_4_tops, scatter_data_melted_6_tops)
all_top_genes = all_top_genes[!duplicated(all_top_genes$Gene_Name),]
all_top_genes$value <- NULL
all_top_genes$variable <- 1




ribodep_all_unique_2hpf_tail <- ribodep_all_unique_2hpf[,c( "Gene_Name", "Median_Length")]
ribodep_all_unique_4hpf_tail <- ribodep_all_unique_4hpf[,c( "Gene_Name", "Median_Length")]
ribodep_all_unique_6hpf_tail <- ribodep_all_unique_6hpf[,c("Gene_Name", "Median_Length")]

columns <-c("Gene_Name")


top_2hpf <- merge(all_top_genes, ribodep_all_unique_2hpf_tail,  by.x=columns, by.y=columns)

top_2hpf_4hpf <- merge(top_2hpf, ribodep_all_unique_4hpf_tail, by.x=columns, by.y=columns,suffixes = c(".2hpf", ".4hpf"))

top_2hpf_4hpf_6hpf <- merge(top_2hpf_4hpf, ribodep_all_unique_6hpf_tail, by.x=columns, by.y=columns)
top_2hpf_4hpf_6hpf$variable <- NULL
colnames(top_2hpf_4hpf_6hpf) <- c("Gene_Name", "Median_Length.2hpf", "Median_Length.4hpf", "Median_Length.6hpf")





