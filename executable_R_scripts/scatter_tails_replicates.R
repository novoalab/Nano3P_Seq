

# Reading arguments from command line
args = commandArgs(trailingOnly=TRUE)

#Arguments
rep1_tail <- read.delim(args[1],sep=",")  #1st variable
rep2_tail <- read.delim(args[2],sep=",")  #2nd variable


rep1_bed <- read.delim(args[3], header=FALSE) # 3rd variable
rep2_bed <- read.delim(args[4], header=FALSE) # 4th variable


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

rep1_tail.processed <- manipulate_tail(rep1_tail)
rep2_tail.processed  <- manipulate_tail(rep2_tail)




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
  return(merged3)
}


rep1.processed <- reshape(rep1_bed,rep1_tail.processed,"Rep1")
rep2.processed <- reshape(rep2_bed,rep2_tail.processed,"Rep2")





# Merge two Replicates
merge_two_rep <- function(rep1, rep2) {
	rep1_unique <-   rep1[!duplicated(rep1[c("Gene_Name", "Sample")]),]
	rep2_unique <-   rep2[!duplicated(rep2[c("Gene_Name", "Sample")]),]
	data_merged <- merge(rep1_unique,rep2_unique, by.x="Gene_Name", by.y="Gene_Name" )
	data_merged_mRNA <-subset(data_merged, Gene_Type.x =="protein_coding")
	data_merged_mRNA_min30 <- subset(data_merged_mRNA, Gene_Count.x > 20 & Gene_Count.y >20)
	return(data_merged_mRNA_min30)
}

both_reps_merged <- merge_two_rep(rep1.processed,rep2.processed)





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


plot_denscols_with_corr_pearson(paste(label,"Rep1_Rep2",sep="_"), both_reps_merged$Median_Length.x , both_reps_merged$Median_Length.y, "Rep1_Median_Length", "Rep2_Median_Length" )
