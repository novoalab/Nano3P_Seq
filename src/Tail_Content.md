# Tail Content
```bash
python /users/enovoa/hliu/analysis/mygit/projects_notes/oz_ribo_switch_find_adapter_in_reads/adapterFinder.py tails.csv cDNA123791_2hpf.fastq 0.7 > cDNA123791_2hpf_tailcontent.csv

python /users/enovoa/hliu/analysis/mygit/projects_notes/oz_ribo_switch_find_adapter_in_reads/adapterFinder.py tails.csv cDNA123791_4hpf.fastq 0.7 > cDNA123791_4hpf_tailcontent.csv

python /users/enovoa/hliu/analysis/mygit/projects_notes/oz_ribo_switch_find_adapter_in_reads/adapterFinder.py tails.csv cDNA123791_6hpf.fastq 0.7 > cDNA123791_6hpf_tailcontent.csv

python /users/enovoa/hliu/analysis/mygit/projects_notes/oz_ribo_switch_find_adapter_in_reads/adapterFinder.py cDNA786327_tails.csv cDNA786327_2hpf.fastq 0.7 > cDNA786327_2hpf_tailcontent.csv

python /users/enovoa/hliu/analysis/mygit/projects_notes/oz_ribo_switch_find_adapter_in_reads/adapterFinder.py cDNA786327_tails.csv cDNA786327_4hpf.fastq 0.7 > cDNA786327_4hpf_tailcontent.csv

python /users/enovoa/hliu/analysis/mygit/projects_notes/oz_ribo_switch_find_adapter_in_reads/adapterFinder.py cDNA786327_tails.csv cDNA786327_6hpf.fastq 0.7 > cDNA786327_6hpf_tailcontent.csv


qsub -cwd -q long-sl7 -l virtual_free=40G tail.sh







sed '$d'  2hpf_nonrRNA_tailcontent.csv >  2hpf_nonrRNA_tailcontent2.csv

```


```R
library(ggplot2)
library(reshape2)

library(stringr)
data.2hpf <- read.delim("cDNA123791_2hpf_tailcontent.csv")


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

processed.2hpf <- manipulate(data.2hpf)







## First making a barplot for the percentage of the bases in 20 nt

A_total <- sum(data3$A)
G_total <- sum(data3$G)
C_total <- sum(data3$C)
T_total <- sum(data3$T)
Total <- A_total+G_total+C_total+T_total
T_freq <- T_total/Total
C_freq <- C_total/Total
G_freq <- G_total/Total
A_freq <- A_total/Total








Frequencies <- as.data.frame(rbind(A_freq,G_freq,C_freq,T_freq))
Frequencies$Base <- rownames(Frequencies)
Frequencies$V1 <- as.numeric(as.character(Frequencies$V1))
Frequencies$V1 <- Frequencies$V1*100


pdf(file="Tail_Content_20nt_base_percentage.pdf",height=6,width=10,onefile=FALSE)
	print(ggplot(data=Frequencies, aes(x=Base, y=V1, fill=Base)) +
		theme_bw()+
		scale_y_continuous(trans = 'log2')+
		ylab("Log Scaled Percentage")+
  		geom_bar(stat="identity"))
	dev.off()




## Looking at each position
per_position_15 <- function(data){
data2 <- data[,c( "readid", "score", "most_common_base", "p15", "p14","p13","p12", "p11", "p10","p9","p8", "p7", "p6","p5","p4", "p3", "p2","p1")]
data3 <- melt(data2, id.vars=c("readid", "score", "most_common_base"))
data4 <- subset(data3, value =="A" | value =="C" |value =="G" |value =="T")
data4$count <- 1
position_sum <- aggregate(.~variable+value, data4[,c("variable", "value","count")], sum)
return(position_sum)
}




per_pos_15.2hpf <- per_position_15(processed.2hpf )







pdf(file="Tail_Content_15nt_Distribution_LinePlot.pdf",height=6,width=10,onefile=FALSE)
	print(ggplot(data=per_pos_15.2hpf , aes(x=variable, y=count, group=value, color=value)) +
    geom_line()+
    theme_bw()+
    scale_y_continuous(trans = 'log2')+
    ylab("Log Scaled Count"))
    dev.off()







pdf(file="Tail_Content_20nt_Distribution.pdf",height=6,width=10,onefile=FALSE)
	print(ggplot(data=data5, aes(x=variable, fill=value)) +
		theme_bw()+
		scale_y_continuous(trans = 'log2')+
		ylab("Log Scaled Percentage")+
  		geom_bar(position="dodge"))
	dev.off()


data_a123 <- subset(data4,  p1 =="A" &  p2 =="A" &  p3 =="A")




data_U_common <- subset(data3, most_common_base=="A")





bed_file <- read.delim("2hpf.genome11_sequin_ALLRNAs_Merged_Rep2.bed",header=FALSE)
bed_file2 <- bed_file[,c("V1", "V4", "V5", "V6", "V16", "V17")]
colnames(bed_file2) <- c("Chr","Read_ID", "Quality", "Strand", "Gene_Name", "Gene_Type")
bed_file3 <- subset(bed_file2, Quality > 30)
merged_bed <- merge(data4, bed_file3, by.y="Read_ID", by.x="readid")


merged_bed_a123 <- subset(merged_bed,  p1 =="A" &  p2 =="A" )



