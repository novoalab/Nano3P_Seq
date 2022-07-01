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


label <- as.character(args[5])  #5th variable



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
  merged2$Time_point <- label
  merged3 <-  merged2[!duplicated(merged2[c("Gene_Name")]),]
  return(merged3)
}


processed.2hpf <- reshape(bed_2hpf,tail.processed,"2hpf")
processed.4hpf <- reshape(bed_4hpf,tail.processed,"4hpf")
processed.6hpf <- reshape(bed_6hpf,tail.processed,"6hpf")



comparison_all <- rbind(processed.2hpf,processed.4hpf,processed.6hpf)


my_comparisons <- list( c("2hpf", "4hpf"), c("2hpf", "6hpf"), c("4hpf", "6hpf") )


pdf(file= paste(label, "dotplot.pdf", sep="_"),height=10,width=8,onefile=FALSE)
      print(ggplot(comparison_all, aes(x=Time_point, y=Median_Length)) + 
        geom_quasirandom(varwidth = TRUE, aes(color=Time_point))+
        geom_boxplot(aes(alpha=0), outlier.shape=NA)+
        stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                geom = "crossbar", width = 0.7, color="#c06c84")+
        theme_bw()+
        ylim(0,180)+
        stat_compare_means(comparisons = my_comparisons, label.y = c(150, 160, 170))+
        #facet_wrap(~Time_point,nrow=1)+
        ggtitle("Zebrafish Embryo")+
        xlab("Time Points")+
        ylab("Tail length") +
        theme(axis.text=element_text(size=14),strip.text = element_text(size=13),
                axis.title=element_text(size=17,face="bold"),
                legend.title = element_text(size = 20),
                legend.text = element_text(color = "black", size=15)))
    dev.off()



