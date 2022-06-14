# File needed : 
## Output of TailfindR

# Reading arguments from command line
args = commandArgs(trailingOnly=TRUE)

#Arguments
data.input <- read.delim(args[1],sep=",")  #1st variable
label <- as.character(args[2])  #1st label


#Create a function
manipulate_tail_<- function(data) { 
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

data.tails_processed <- manipulate_tail_(data.tails)


#Export the final table
write.table(data.tails_processed, file=paste(label, ".tails.processed.tsv", sep=""), quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)