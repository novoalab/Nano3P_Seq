
library(ggplot2)
data <- read.delim("ratio.tsv")



#### Gene Type Count Comparison ###
simple_barplot_grouped <- function(data, label){
	pdf(file=paste(label, ".pdf",sep="_"),height=4,width=7,onefile=FALSE)
	print(ggplot(data, aes( y=Value+1, x=Sample, fill=RNA)) + 
    geom_bar(position="dodge", stat="identity")+
    scale_y_continuous(trans = 'log2')+
    xlab("Condiiton")+
    ylab("Abundance")+
    theme_bw())
    dev.off()
}



simple_barplot_grouped(data, "Mt_rRNA_Ratio")
