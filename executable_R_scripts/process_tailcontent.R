# Process tail content

library(stringr)
library(Biostrings)


# Reading arguments from command line
args = commandArgs(trailingOnly=TRUE)

#Arguments
tails <- read.delim(args[1],sep=",")  #1st variable
bed <- read.delim(args[2], header=FALSE) # 2nd variable
content <- read.delim(args[3], header=FALSE) # 3rd variable

label <- as.character(args[4])  #4th variable





###################################################################
##########. 1. TAIL PREDICTIONS ###################################
###################################################################
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
  data_final <- data_filt_T[,c("read_id", "tail_length")]
  return(data_final)
}

tail.processed <- manipulate_tail(tails)


###################################################################
###################################################################






###################################################################
##########     2. BED FILES     ###################################
###################################################################
process.bed <- function(data){
	#Select columns
	data2<- data[,c("V4","V16","V17")]
	#Rename columns
	colnames(data2) <- c("read_id", "Gene_Name", "Gene_Type")
	return(data2)
}

bed_file.processed <- process.bed(bed)

###################################################################
###################################################################





###################################################################
##########     3. CONTENT FILES ###################################
###################################################################

content_process <- function(data,label) {
	##### Positive Strand #####
	pos <- subset(data, strand=="+")
		#Keep the columns interested
		pos2 <- pos[,c("read", "strand", "softclip3end", "sc3ATCG")]
		#Rename columns
		colnames(pos2) <- c("read", "strand", "softclip", "ATCG")
		#Replace bases with complementary
		pos2$softclip <- gsub("T", "X", pos2$softclip)
		pos2$softclip <- gsub("A", "Y", pos2$softclip)
		pos2$softclip <- gsub("G","Z", pos2$softclip)
		pos2$softclip <- gsub("C","B", pos2$softclip)
		pos2$softclip <- gsub("X", "A", pos2$softclip)
		pos2$softclip <- gsub("Y", "U", pos2$softclip)
		pos2$softclip <- gsub("Z","C", pos2$softclip)
		pos2$softclip <- gsub("B","G", pos2$softclip)
		#Split the ATCG column into 4
		bases <- str_split_fixed(pos2$ATCG, "-", 4)
		#Rename columns
		colnames(bases) <- c("U", "A","G","C")
		#Reposition columns
		bases <- bases[,c("A", "U", "C", "G")]
		#Bind it to the original data
		pos3 <- cbind(pos2, bases)
		#Convert the columns into numeric
		pos3$A <- as.numeric(as.character(pos3$A))
		pos3$U <- as.numeric(as.character(pos3$U))
		pos3$G <- as.numeric(as.character(pos3$G))
		pos3$C <- as.numeric(as.character(pos3$C))
		#Convert NA into 0 
		pos3[["A"]][is.na(pos3[["A"]])] <- 0
		pos3[["U"]][is.na(pos3[["U"]])] <- 0
		pos3[["G"]][is.na(pos3[["G"]])] <- 0
		pos3[["C"]][is.na(pos3[["C"]])] <- 0
		#Make SUM of the Counts
		pos3$Sum <- pos3$A + pos3$U +pos3$G + pos3$C
		#calculate frequency
		pos3$A_freq <- pos3$A/pos3$Sum
		pos3$U_freq <- pos3$U/pos3$Sum
		pos3$G_freq <- pos3$G/pos3$Sum
		pos3$C_freq <- pos3$C/pos3$Sum
		#Create a column
		pos3$softclip_final<- pos3$softclip


	##### Negative Strand #####
	neg <- subset(data, strand =="-")
		#Keep the columns interested
		neg2 <- neg[,c("read", "strand", "softclip3end", "sc3ATCG")]
		#Rename columns
		colnames(neg2) <- c("read", "strand", "softclip", "ATCG")
		#Split the ATCG column into 4
		bases <- str_split_fixed(neg2$ATCG, "-", 4)
		#Rename columns
		colnames(bases) <- c("A", "U","C","G")
		#Convert T into U
		neg2$softclip <- gsub("T", "U", neg2$softclip)
		#Bind it to the original data
		neg3 <- cbind(neg2, bases)
		#Convert the columns into numeric
		neg3$A <- as.numeric(as.character(neg3$A))
		neg3$U <- as.numeric(as.character(neg3$U))
		neg3$G <- as.numeric(as.character(neg3$G))
		neg3$C <- as.numeric(as.character(neg3$C))
		#Convert NA into 0
		neg3[["A"]][is.na(neg3[["A"]])] <- 0
		neg3[["U"]][is.na(neg3[["U"]])] <- 0
		neg3[["G"]][is.na(neg3[["G"]])] <- 0
		neg3[["C"]][is.na(neg3[["C"]])] <- 0
		#Make SUM of the Counts	
		neg3$Sum <- neg3$A + neg3$U +neg3$G + neg3$C
		#calculate frequency
		neg3$A_freq <- neg3$A/neg3$Sum
		neg3$U_freq <- neg3$U/neg3$Sum
		neg3$G_freq <- neg3$G/neg3$Sum
		neg3$C_freq <- neg3$C/neg3$Sum
		#Reverse the softclip column
		neg3$softclip <- as.character(neg3$softclip)
		neg3$softclip_final<- reverse(neg3$softclip)
	### Bind both strand reads
	both_strands <- rbind(pos3, neg3)
	both_strands$softclip <- NULL
	both_strands$Timepoint <- label
	return(both_strands)
}



content.processed <- content_process(content,label)


###################################################################
###################################################################



###################################################################
##########     4. MERGE WITH TAIL LENGTH AND BED ##################
###################################################################

merge_tail_bed <- function(data) {
	content.bed <- merge(data, bed_file.processed, by.x="read", by.y="read_id")
	print(dim(content.bed))
	content.bed.tail <- merge(content.bed, tail.processed,by.x="read", by.y="read_id")
	print(dim(content.bed.tail))
	content_all_mRNA <- subset(content.bed.tail, Gene_Type=="protein_coding")
	return(content_all_mRNA)
}


tail_bed_content <- merge_tail_bed(content.processed)



###################################################################
###################################################################



###################################################################
##########     5. TAKE AU RICH TAILS     ##########################
###################################################################

au_rich <- function(data) {
	content_all_mRNA.AUrich <- subset(data, A_freq+U_freq > 0.8)
	#Split columns
	nucleotides <-  str_split_fixed(content_all_mRNA.AUrich$softclip_final, "", 11)
	#Rename columns
	colnames(nucleotides) <-c("p1", "p2","p3","p4", "p5", "p6", "p7", "p8", "p9", "p10", "rest")
	#Select four columns
	nucleotides2 <- nucleotides[,c("p1", "p2","p3","p4", "p5", "p6", "p7", "p8", "p9", "p10")]
	#Attach columns to the original data
	content_all_mRNA.AUrich2 <- cbind(content_all_mRNA.AUrich, nucleotides2)
	return(content_all_mRNA.AUrich2)
}

au_rich <- au_rich(tail_bed_content)

write.table(au_rich, file=paste(label, "AU_rich_tail_content.tsv",sep="_"), quote=FALSE, sep="\t", row.names=FALSE)

