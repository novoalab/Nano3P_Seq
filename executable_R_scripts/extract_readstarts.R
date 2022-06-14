# File needed : 
## BED converted from BAM mapped to genome and rRNAs and sorted


# Reading arguments from command line
args = commandArgs(trailingOnly=TRUE)

#Arguments
data.input <- read.delim(args[1],header=FALSE)  #1st variable
label <- as.character(args[2])  #1st label


#Create a function
manipulate <- function(data) {
	#Subset the negative strand reads
	subs_neg <- subset(data, V6 =="-")
	#Keep the columns we are interested in (Including read start position)
	subs_neg2 <- subs_neg[,c("V1", "V3", "V6", "V4")]
	#Rename the columns
	colnames(subs_neg2) <- c("Chr", "Position", "Strand", "Read_ID")
	#Subset the positive strand reads
	subs_pos <- subset(data, V6 =="+")
	#Keep the columns we are interested in (Including read start position)
	subs_pos2 <- subs_pos[,c("V1", "V2", "V6", "V4")]
	#Rename the columns
	colnames(subs_pos2) <- c("Chr", "Position", "Strand", "Read_ID")
	#Merge the tables
	merged <- rbind(subs_neg2,subs_pos2 )
	#Create two position columns +/- 10 nt from read start
	merged$Position1 <- merged$Position-10
	merged$Position2 <-  merged$Position+10
	#Choose the columns we are interested in
	final <- merged[,c("Chr","Position1", "Position2",  "Read_ID","Strand")]
	#Make sure they dont coincide with chromosome ends
	final2 <- subset(final, Position1 > 0 & Position2 > 0)
	#Return the final table
	return(final2)
}

data.processed <- manipulate(data.input)

#Export the final table
write.table(data.processed, file=paste(label, ".readstarts.bed", sep=""), quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)