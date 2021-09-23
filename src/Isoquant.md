```bash
#Trying out IsoQuant
isoquant=/users/enovoa/boguzhan/Software/IsoQuant/isoquant.py

gtf=/users/enovoa/boguzhan/references/danio_rerio/Danio_rerio.GRCz11.103.2.gtf

python $isoquant --genedb $gtf --complete_genedb --bam 246hpf_nonrRNA_both_aligned.bam --data_type nanopore -o OUTPUT_FOLDER





```


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
library(ggforce)

#reads
reads <- read.delim("00_246hpf_nonrRNA_both_aligned.read_assignments.tsv")
reads2 <- subset(reads, assignment_type!="ambiguous" &  assignment_type!="noninformative" )

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

ribodep_hpf2_rep1.data <- read.delim("2hpf.genome11_sequin_ALLRNAs_Merged_Rep1.bed", header=FALSE)
ribodep_hpf2_rep2.data <- read.delim("2hpf.genome11_sequin_ALLRNAs_Merged_Rep2.bed", header=FALSE)

ribodep_hpf4_rep1.data <- read.delim("4hpf.genome11_sequin_ALLRNAs_Merged_Rep1.bed", header=FALSE)
ribodep_hpf4_rep2.data <- read.delim("4hpf.genome11_sequin_ALLRNAs_Merged_Rep2.bed", header=FALSE)

ribodep_hpf6_rep1.data <- read.delim("6hpf.genome11_sequin_ALLRNAs_Merged_Rep1.bed", header=FALSE)
ribodep_hpf6_rep2.data <- read.delim("6hpf.genome11_sequin_ALLRNAs_Merged_Rep2.bed", header=FALSE)

ribodep_hpf2_merged.data<- rbind(ribodep_hpf2_rep1.data ,ribodep_hpf2_rep2.data )
ribodep_hpf4_merged.data<- rbind(ribodep_hpf4_rep1.data ,ribodep_hpf4_rep2.data )
ribodep_hpf6_merged.data<- rbind(ribodep_hpf6_rep1.data ,ribodep_hpf6_rep2.data )


# Reshape the tables and remove low quality reads
reshape<- function(data,tails,label) {
  data2 <- data[,c("V1", "V4", "V5", "V6", "V16", "V17")]
  colnames(data2) <- c("Chr","Read_ID", "Quality", "Strand", "Gene_Name", "Gene_Type")
  data3 <- subset(data2, Quality > 30)
  #Merge the data with tails
  merged <- merge(data3, tails, by.x="Read_ID",by.y=c("read_id") )
  merged <- merged[,c("Read_ID", "Chr", "Gene_Name", "Gene_Type", "tail_length")]
  return(merged)
}

ribodep_hpf2_merged.reshape <- reshape(ribodep_hpf2_merged.data,ribodepleted_merged.tails_processed,"Ribodep_2hpf_merged")
ribodep_hpf4_merged.reshape <- reshape(ribodep_hpf4_merged.data,ribodepleted_merged.tails_processed,"Ribodep_4hpf_merged")
ribodep_hpf6_merged.reshape <- reshape(ribodep_hpf6_merged.data,ribodepleted_merged.tails_processed,"Ribodep_6hpf_merged")



merge_reads <- function(data,reads,label) {
	reads_bed_tail <- merge(data, reads, by.x="Read_ID", by.y="X.read_id")
	reads_bed_tail2<- subset(reads_bed_tail, Gene_Type=="protein_coding")
	reads_bed_tail3 <- reads_bed_tail2[,c("Read_ID","isoform_id","gene_id","tail_length","assignment_type")]
    #remove transcripts less than x coverage
    df <- vector()
    for (trans in unique(reads_bed_tail3$isoform_id)){
    	subs <- subset(reads_bed_tail3, isoform_id==trans)
    	if (nrow(subs) > 10) {
    		df <- rbind(df, subs) } else {}
    	}
    #remove genes with only one isoform
    df2 <- vector()
   for (genes in unique(df$gene_id)){
    	subs <- subset(df, gene_id==genes)
    	if (length(unique(subs$isoform_id)) > 1) {
    		df2 <- rbind(df2,subs) } else {}
    	}
    #Count isoform
    df3 <- vector()
    for (genes in unique(df2$gene_id)){
    	subs <- subset(df2, gene_id==genes)
    	iso_count <- length(unique(subs$isoform_id)) 
    	subs$isoform_count <- iso_count
    	df3 <- rbind(df3, subs)
    }
    #median length per isoform
    df4 <- vector()
    for (trans in unique(df3$isoform_id)){
    	subs <- subset(df3, isoform_id==trans)
    	median_length <- median(subs$tail_length)
    	subs$median_length_isoform <- median_length 
    	df4 <- rbind(df4, subs)
    }
    write.table(df4, file=paste(label,"processed_tail_isoform.tsv", sep="_"), quote=FALSE, sep="\t")
    return(df4)
}

ribodep_hpf2_merged_READS <- merge_reads(ribodep_hpf2_merged.reshape,reads, "2HPF" )
ribodep_hpf4_merged_READS <- merge_reads(ribodep_hpf4_merged.reshape,reads, "4HPF" )
ribodep_hpf6_merged_READS <- merge_reads(ribodep_hpf6_merged.reshape,reads, "6HPF" )


ribodep_hpf2_merged_READS <- read.delim("2HPF_processed_tail_isoform.tsv")
ribodep_hpf4_merged_READS <- read.delim("4HPF_processed_tail_isoform.tsv")
ribodep_hpf6_merged_READS <- read.delim("6HPF_processed_tail_isoform.tsv")






#Kruskal-Wallis test
  compute_kruskal <- function(data,label) {
    df <- vector()
    for (genes in unique(data$gene_id)){
      subs <- subset(data, gene_id==genes)
      kruskal <- kruskal.test(tail_length ~ isoform_id, data = subs)
      p_value <- kruskal$p.value
      isoform_count <- unique(subs$isoform_count)
      table <- as.data.frame(cbind(genes,p_value,isoform_count ))
      df <- rbind(df, table)
    }
    df$p_value <- as.numeric(as.character(df$p_value))
    df$isoform_count <- as.numeric(as.character(df$isoform_count)) 
    df$p.adjust <- p.adjust(df$p_value,method="hochberg")
    df <- df[order(df$p.adjust),]
    df$TimePoint <- label
    return(df)

}

ribodep_hpf2_merged_READS_kruskal <- compute_kruskal(ribodep_hpf2_merged_READS, "2hpf")
ribodep_hpf4_merged_READS_kruskal <- compute_kruskal(ribodep_hpf4_merged_READS, "4hpf")
ribodep_hpf6_merged_READS_kruskal <- compute_kruskal(ribodep_hpf6_merged_READS, "6hpf")

ribodep_hpf246_merged_READS_kruskal <- rbind(ribodep_hpf2_merged_READS_kruskal,ribodep_hpf4_merged_READS_kruskal,ribodep_hpf6_merged_READS_kruskal)


write.table(ribodep_hpf246_merged_READS_kruskal, file="All_timepoints_kruskal.tsv", sep="\t", quote=FALSE,row.names=FALSE)



dotplot_kruskal <- function(data,label) {
pdf(file=paste(label, "kruskal_isoform_dotplot.pdf",sep="_"),height=10,width=12,onefile=FALSE)
      print(ggplot(data, aes(x=as.factor(isoform_count), y=log(1/p_value))) + 
        geom_quasirandom(varwidth = TRUE, aes())+
        geom_boxplot(aes(alpha=0), outlier.shape=NA)+
        stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                geom = "crossbar", width = 0.7, color="#c06c84")+
        geom_text_repel(data=subset(data, log(1/p_value) > 8 ), aes(y = log(1/p_value), label = genes), color = "black", size=6)+
        stat_n_text(size=5)+
        theme_bw()+
        ggtitle(label)+
        xlab("Isoform count")+
        ylab("log(1/Kruskal-Wallis p-value) ") +
        theme(axis.text=element_text(size=14),strip.text = element_text(size=13),
                axis.title=element_text(size=17,face="bold"),
                legend.title = element_text(size = 20),
                legend.text = element_text(color = "black", size=15)))
    dev.off()
}



dotplot_kruskal(ribodep_hpf2_merged_READS_kruskal,"2hpf")
dotplot_kruskal(ribodep_hpf4_merged_READS_kruskal,"4hpf")
dotplot_kruskal(ribodep_hpf6_merged_READS_kruskal,"6hpf")



## Dotplot kruskal all the timepoints
boxplot_kruskal_merged <- function(data,label) {
pdf(file=paste(label, "kruskal_isoform_boxplot.pdf",sep="_"),height=4,width=7,onefile=FALSE)
      print(ggplot(data, aes(x=as.factor(TimePoint), y=log(1/p.adjust),fill=TimePoint))+
        geom_boxplot(notch=TRUE) +
        stat_n_text(size=3)+
        theme_bw()+
        geom_hline(yintercept=log(1/0.05), linetype="dashed", size=0.5,colour="red")+
        ggtitle(label)+
        ylab("log(1/p_adjusted)")+
        facet_zoom(ylim = c(0, 10))+
        xlab("TimePoints") +
        theme(axis.text=element_text(size=14),strip.text = element_text(size=13),
                axis.title=element_text(size=17,face="bold"),
                legend.title = element_text(size = 20),
                legend.text = element_text(color = "black", size=15)))
    dev.off()
}


boxplot_kruskal_merged(ribodep_hpf246_merged_READS_kruskal, "246HPF")



## STATS OF SIGNIFICANT GENES


ribodep_hpf2_merged_READS_kruskal_sign <- subset(ribodep_hpf2_merged_READS_kruskal, p.adjust < 0.05)
ribodep_hpf4_merged_READS_kruskal_sign <- subset(ribodep_hpf4_merged_READS_kruskal, p.adjust < 0.05)
ribodep_hpf6_merged_READS_kruskal_sign <- subset(ribodep_hpf6_merged_READS_kruskal, p.adjust < 0.05)



library(VennDiagram)


list_all_timepoints_significant <- list(ribodep_hpf2_merged_READS_kruskal_sign[,1], ribodep_hpf4_merged_READS_kruskal_sign[,1],ribodep_hpf6_merged_READS_kruskal_sign[,1])
Reduce(intersect, list_all_timepoints_significant)

venn.diagram(x =list_all_timepoints_significant,
     category.names = c("2HPF" , "4HPF", "6HPF"),
     filename = 'Venn_All_Timepoints_Significant.tiff',
     lwd = 2,
     lty = 'blank',
     fill = c("#f5a25d", "#bedbbb", "#776d8a"),
     output=TRUE,
     # Output features
    imagetype="tiff" ,
    height = 6000 , 
    width = 6000 , 
    resolution = 6000,
    compression = "lzw",
   # Numbers
    cex = .5,
    fontface = "bold",
    fontfamily = "sans",
    # Set names
    cat.cex = 0.3,
    cat.fontface = "bold",
    cat.default.pos = "outer",
    cat.pos = c(-25, 0, 25 ),
    cat.dist = c(0.15, 0.15 ,-0.5),
    cat.fontfamily = "sans",
    #rotation = 1
)





list_all_timepoints <- list(ribodep_hpf2_merged_READS_kruskal[,1], ribodep_hpf4_merged_READS_kruskal[,1],ribodep_hpf6_merged_READS_kruskal[,1])
Reduce(intersect, list_all_timepoints)

venn.diagram(x =list_all_timepoints,
     category.names = c("2HPF" , "4HPF", "6HPF"),
     filename = 'Venn_All_Timepoints.tiff',
     lwd = 2,
     lty = 'blank',
     fill = c("#f5a25d", "#bedbbb", "#776d8a"),
     output=TRUE,
     # Output features
    imagetype="tiff" ,
    height = 6000 , 
    width = 6000 , 
    resolution = 6000,
    compression = "lzw",
   # Numbers
    cex = .5,
    fontface = "bold",
    fontfamily = "sans",
    # Set names
    cat.cex = 0.3,
    cat.fontface = "bold",
    cat.default.pos = "outer",
    cat.pos = c(-25, 0, 25 ),
    cat.dist = c(0.15, 0.15 ,-0.5),
    cat.fontfamily = "sans",
    #rotation = 1
)






















dotplot_individual <- function(data,label) {
for (Gene in unique(data$gene)) {
		subs <- subset( data, gene_id==Gene)
pdf(file=paste(label,Gene, "tail_length_per_isoform_dotplot.pdf",sep="_"),height=5,width=12,onefile=FALSE)
			print(ggplot(subs, aes(x=isoform_id, y=tail_length)) + 
				geom_quasirandom(varwidth = TRUE, aes())+
				geom_boxplot(aes(alpha=0), outlier.shape=NA)+
				stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
            		geom = "crossbar", width = 0.7, color="#c06c84")+
				stat_n_text(size=5)+
				theme_bw()+
				ggtitle(label)+
				xlab("Isoform")+
              	ylab("Tail length") +
				theme(axis.text=element_text(size=14),strip.text = element_text(size=13),
            		axis.title=element_text(size=17,face="bold"),
            		legend.title = element_text(size = 20),
            		legend.text = element_text(color = "black", size=15)))
		dev.off()
}
}


dotplot_individual(ribodep_hpf2_merged_READS,"2hpf")
dotplot_individual(ribodep_hpf4_merged_READS,"4hpf")
dotplot_individual(ribodep_hpf6_merged_READS,"6hpf")



ribodep_hpf2_merged_READS$Sample <- "HPF2"
ribodep_hpf4_merged_READS$Sample <- "HPF4"
ribodep_hpf6_merged_READS$Sample <- "HPF6"

ribodep_hpf246_merged_READS <- rbind(ribodep_hpf2_merged_READS,ribodep_hpf4_merged_READS,ribodep_hpf6_merged_READS)




library(ggpubr)


dotplot_merged <- function(data,label) {
for (Gene in unique(data$gene)) {
    subs <- subset( data, gene_id==Gene)
    group.colors <- c(HPF2 = "#F8766D", HPF4 = "#00BA38", HPF6 ="#619CFF")
  pdf(file=paste(label,Gene, "merged_tail_length_per_isoform_dotplot.pdf",sep="_"),height=7,width=15,onefile=FALSE)
      print(ggplot(subs, aes(x=isoform_id, y=tail_length, colour=Sample)) + 
        geom_quasirandom(varwidth = TRUE, aes())+
        geom_boxplot(aes(alpha=0), outlier.shape=NA)+
        stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                geom = "crossbar", width = 0.7, color="#c06c84")+
        stat_n_text(size=5)+
        theme_bw()+
        ggtitle(label)+
        scale_colour_manual(values=group.colors)+
        facet_wrap(~Sample, nrow=1)+
        xlab("Isoform")+
        ylab("Tail length") +
        theme(axis.text=element_text(size=14),strip.text = element_text(size=13),
                axis.title=element_text(size=17,face="bold"),
                legend.title = element_text(size = 20),
                legend.text = element_text(color = "black", size=15)))
    dev.off()
}
}

dotplot_merged(ribodep_hpf246_merged_READS,"246hpf")



dotplot_merged_specific_gene <- function(data,label,Gene,h,w) {
    subs <- subset( data, gene_id==Gene)
    group.colors <- c(HPF2 = "#F8766D", HPF4 = "#00BA38", HPF6 ="#619CFF")
  pdf(file=paste(label,Gene, "merged_tail_length_per_isoform_dotplot.pdf",sep="_"),height=h,width=w,onefile=FALSE)
      print(ggplot(subs, aes(x=isoform_id, y=tail_length, colour=Sample)) + 
        geom_quasirandom(varwidth = TRUE, aes())+
        geom_boxplot(aes(alpha=0), outlier.shape=NA)+
        stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                geom = "crossbar", width = 0.7, color="#c06c84")+
        stat_n_text(size=5)+
        theme_bw()+
        ggtitle(label)+
        scale_colour_manual(values=group.colors)+
        facet_wrap(~Sample, nrow=1)+
        xlab("Isoform")+
        ylab("Tail length") +
        theme(axis.text=element_text(size=14),strip.text = element_text(size=13),
                axis.title=element_text(size=17,face="bold"),
                legend.title = element_text(size = 20),
                legend.text = element_text(color = "black", size=15)))
    dev.off()
}

dotplot_merged_specific_gene(ribodep_hpf246_merged_READS,"246hpf","ENSDARG00000040184",4, 8.5 )
dotplot_merged_specific_gene(ribodep_hpf246_merged_READS,"246hpf","ENSDARG00000035679",4, 8.5 )
















boxplot_merged <- function(data,label) {
for (Gene in unique(data$gene)) {
    subs <- subset( data, gene_id==Gene)
    pdf(file=paste(label,Gene, "merged_tail_length_per_isoform_boxplot.pdf",sep="_"),height=5,width=15,onefile=FALSE)
    print(ggboxplot(subs, x = "isoform_id", y = "tail_length",
          color = "Sample", palette = "jco",
          add = "jitter",
          facet.by = "Sample", short.panel.labs = FALSE)+
   stat_compare_means(label = "p.format")+
    stat_compare_means(label =  "p.signif", label.x = 1.5))
   dev.off()
 }
}
boxplot_merged(ribodep_hpf246_merged_READS,"246hpf")

















#EXTRACT READS INDIVIDUALLY `PER TRANSCRIPT
extract_reads <- function(data, genes, label) {
  subs <- subset(data, gene_id==genes)
  for (trans in unique(subs$isoform_id)) {
    subs2 <- subset(subs, isoform_id == trans)
    subs3 <- subs2[,"Read_ID"]
    write.table(subs3, file=paste(label, genes, trans, "reads.tsv", sep="_"), col.names=FALSE, row.names=FALSE, quote=FALSE)
  }
}


extract_reads(ribodep_hpf2_merged_READS,"ENSDARG00000040184", "2HPF")







#Compare p values
compare_pvalues <- function(data1, data2, label1, label2) {
merged <- merge(data1, data2, by.x="genes", by.y="genes")
pdf(file=paste(label1,label2, "kruskal_scatter.pdf",sep="_"),height=6,width=6,onefile=FALSE)
print(ggplot(merged, aes(log(1/p_value.x),log(1/p_value.y)))+
	xlab(label1)+
	ylab(label2)+
	geom_text_repel(data=subset(merged, log(1/p_value.x) > 15 ), aes(x = log(1/p_value.x), y=log(1/p_value.y), label = genes), color = "black", size=6)+
	#ylim(0,30)+
	#xlim(0,30)+
	ggtitle("Kruskal-Wallis P value")+
 	geom_point())
dev.off()
}

compare_pvalues(ribodep_hpf2_merged_READS_kruskal,ribodep_hpf4_merged_READS_kruskal, "2HPF", "4HPF")
compare_pvalues(ribodep_hpf4_merged_READS_kruskal,ribodep_hpf6_merged_READS_kruskal, "4HPF", "6HPF")












#### COMPUTE KRUSKAL WALLIS BY ISOFORM
ribodep_hpf2_merged_READS$Sample <- "2HPF"
ribodep_hpf4_merged_READS$Sample <- "4HPF"
ribodep_hpf6_merged_READS$Sample <- "6HPF"

ribodep_hpf246_merged_READS <- rbind(ribodep_hpf2_merged_READS,ribodep_hpf4_merged_READS,ribodep_hpf6_merged_READS)

transcript <- "ENSDART00000158291"

#Kruskal-Wallis test
  compute_kruskal_per_transcript <- function(data) {
    df <- vector()
    for (transcript in unique(data$isoform_id)){
      subs <- subset(data, isoform_id==transcript)
      if (length(unique(subs$Sample)) >1 ) {
        kruskal <- kruskal.test(tail_length ~ Sample, data = subs)
        chi<- as.numeric(kruskal$statistic)
        p_value <- kruskal$p.value
        #isoform_count <- unique(subs$isoform_count)
        gene <- unique(subs$gene_id)
        table <- as.data.frame(cbind(transcript,p_value,chi, as.character(gene)))
        df <- rbind(df, table)
    } else {}
  }
    colnames(df) <- c("transcript", "p_value","chi", "gene")
    df$p_value <- as.numeric(as.character(df$p_value))
    #df$isoform_count <- as.numeric(as.character(df$isoform_count))
    df$chi <- as.numeric(as.character(df$chi))    
    df <- df[order(df$p_value),]
    df <- df[order(df$chi),]
    df <- df[order(df$gene),]
    return(df)
}

ribodep_hpf246_merged_kruskal_per_transcript <- compute_kruskal_per_transcript(ribodep_hpf246_merged_READS)

 ribodep_hpf246_merged_kruskal_per_transcript$p.adjust <- p.adjust(ribodep_hpf246_merged_kruskal_per_transcript$p_value,method="hochberg")


write.table(ribodep_hpf246_merged_kruskal_per_transcript, file="Comparison_across_timepoints_kruskal_padj.tsv", sep="\t", quote=FALSE,row.names=FALSE)



compute_variance_per_gene <- function(data) {
  df <- vector()
for (genes in unique(data$gene)) {
        subs <- subset(data, gene==genes)
        subs$var <- var(log(1/subs$p.adjust))
        df <- rbind(df, subs)
      }

        df <- df[order(df$var),]
       # df <- df[order(df$gene),]
return(df)
}

ribodep_hpf246_merged_kruskal_per_transcript_var <- compute_variance_per_gene(ribodep_hpf246_merged_kruskal_per_transcript)
ribodep_hpf246_merged_kruskal_per_transcript_var$Xlab <- "1"
ribodep_hpf246_merged_kruskal_per_transcript_var <- ribodep_hpf246_merged_kruskal_per_transcript_var[!duplicated(ribodep_hpf246_merged_kruskal_per_transcript_var$gene), ]


write.table(ribodep_hpf246_merged_kruskal_per_transcript_var, file="Comparison_across_timepoints_kruskal_var.tsv", sep="\t", quote=FALSE,row.names=FALSE)


pdf(file="Dotplot_Variance_Per_Gene.pdf",height=8,width=3,onefile=FALSE)
      print(ggplot(ribodep_hpf246_merged_kruskal_per_transcript_var, aes(x=Xlab, y=log(var+1))) + 
        geom_quasirandom(varwidth = TRUE, aes())+
        geom_boxplot(aes(alpha=0), outlier.shape=NA)+
        stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                geom = "crossbar", width = 0.7, color="#c06c84")+
       # geom_text_repel(data=subset(ribodep_hpf246_merged_kruskal_per_transcript_var, var > 1000 ), aes(x=Xlab, y = var, label = gene), color = "black", size=6)+
        stat_n_text(size=5)+
        theme_bw()+
        #xlab("Isoform count")+
        ylab("log(Variance of p-adj") +
        theme(axis.text=element_text(size=14),strip.text = element_text(size=13),
                axis.title=element_text(size=17,face="bold"),
                legend.title = element_text(size = 20),
                legend.text = element_text(color = "black", size=15)))
    dev.off()
















dotplot_merged_time <- function(data,label) {
for (Gene in unique(data$gene_id)) {
    subs <- subset( data, gene_id==Gene)
pdf(file=paste(label,Gene, "tail_length_per_isoform_dotplot.pdf",sep="_"),height=5,width=12,onefile=FALSE)
      print(ggplot(subs, aes(x=Sample, y=tail_length, colour=Sample)) + 
        geom_quasirandom(varwidth = TRUE, aes())+
        geom_boxplot(aes(alpha=0), outlier.shape=NA)+
        stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                geom = "crossbar", width = 0.7, color="#c06c84")+
        stat_n_text(size=5)+
        theme_bw()+
        facet_wrap(~isoform_id,nrow=1)+
        ggtitle(label)+
        xlab("TimePoint")+
        ylab("Tail length") +
        theme(axis.text=element_text(size=14),strip.text = element_text(size=13),
                axis.title=element_text(size=17,face="bold"),
                legend.title = element_text(size = 20),
                legend.text = element_text(color = "black", size=15)))
    dev.off()
}
}


dotplot_merged_time(ribodep_hpf246_merged_READS,"246hpf")


library(ggpubr)

boxplot_merged_time <- function(data,label) {
for (Gene in unique(data$gene)) {
    subs <- subset( data, gene_id==Gene)
    my_comparisons <- list( c("2HPF", "4HPF"), c("4HPF", "6HPF"), c("2HPF", "6HPF") )
    pdf(file=paste(label,Gene, "tail_length_per_isoform_boxplot.pdf",sep="_"),height=5,width=15,onefile=FALSE)
    print(ggboxplot(subs, x = "Sample", y = "tail_length",
          color = "Sample", palette = "jco",
          add = "jitter",
          facet.by = "isoform_id", short.panel.labs = FALSE)+
   stat_compare_means(comparisons = my_comparisons, label = "p.format")+
    stat_compare_means(comparisons = my_comparisons, label =  "p.signif", label.x = 1.5))
   dev.off()
 }
}
boxplot_merged_time(ribodep_hpf246_merged_READS,"246hpf")


```



```BASH
java -jar /users/enovoa/boguzhan/Software/picard/build/libs/picard.jar FilterSamReads \
       I=2hpf_nonrRNA_both_trimmed_porechop.genome11_sequin.sorted.bam \
       O=2HPF_ENSDARG00000040184_ENSDART00000058779.bam\
       READ_LIST_FILE=2HPF_ENSDARG00000040184_ENSDART00000058779_reads.tsv\
       FILTER=includeReadList

java -jar /users/enovoa/boguzhan/Software/picard/build/libs/picard.jar FilterSamReads \
       I=2hpf_nonrRNA_both_trimmed_porechop.genome11_sequin.sorted.bam \
       O=2HPF_ENSDARG00000040184_ENSDART00000152269.bam\
       READ_LIST_FILE=2HPF_ENSDARG00000040184_ENSDART00000152269_reads.tsv\
       FILTER=includeReadList



samtools index 2HPF_ENSDARG00000040184_ENSDART00000058779.bam
samtools index 2HPF_ENSDARG00000040184_ENSDART00000152269.bam




