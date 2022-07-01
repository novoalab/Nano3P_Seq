#Script to extracy read starts from BED

args <- commandArgs(trailingOnly = TRUE)
data <- read.delim(args[1], header=FALSE)
label <- args[2] #2nd variable
manipulate <- function(data) {
       subs_neg <- subset(data, V6 =="\"-\"")
       subs_neg2 <- subs_neg[,c('V1', 'V3', 'V6', 'V4')]
       colnames(subs_neg2) <- c('Chr', 'Position', 'Strand', 'Read_ID')
       subs_pos <- subset(data, V6 =="\"+\"")
       subs_pos2 <- subs_pos[,c('V1', 'V2', 'V6', 'V4')]
       colnames(subs_pos2) <- c('Chr', 'Position', 'Strand', 'Read_ID')
       merged <- rbind(subs_neg2,subs_pos2 )
       merged\$Position1 <- merged\$Position-10
       merged\$Position2 <-  merged\$Position+10
       final <- merged[,c('Chr','Position1', 'Position2',  'Read_ID','Strand')]
       final2 <- subset(final, Position1 > 0 & Position2 > 0)
       return(final2)
}
data.processed <- manipulate(data)
write.table(data.processed, file=paste(label, 'readstarts.bed', sep='.'), quote=FALSE, sep='\\t', col.names=FALSE, row.names=FALSE)