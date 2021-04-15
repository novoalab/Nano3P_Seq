###################################################
###### ANALYSIS OF THE SEQUINS TAILFINDR #########
###################################################
########## OGUZHAN BEGIK APRIL 2020 ###############


```R
library(stringr)
library(ggplot2)
library(ggbeeswarm)
library(ggrepel)
library(scales)
library(dplyr)
library(EnvStats)
library(reshape2)



######## PROCESS THE cDNA TAIL FILE

#Import tail
cDNA_tails <- read.delim("cDNA964321_tails.csv", sep=",")

manipulate_tail <- function(data) { 
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

cDNA_tails.processed <- manipulate_tail(cDNA_tails)




cDNA.input <- read.delim("cDNA964321_SequinOnly.bed", header=FALSE)

cDNA_sequin <- cDNA.input[,c("V4", "V5", "V16", "V17")]
colnames(cDNA_sequin) <- c("Read_ID", "Quality", "Gene_Name", "Gene_Type")
cDNA_sequin$Group <- gsub("_.*","",cDNA_sequin$Gene_Name)





# First lets cleanup the data (Remove any gene that has less than 20 reads)
cleanup<- function(data) {
	data2 <- subset(data, Quality > 50)
	final<- vector()
		for (gene in unique(data2$Gene_Name)) {
			subs <- subset(data2, Gene_Name==gene)
			if (nrow(subs) > 30) {
				final <- rbind(final, subs)
			} else {
				subs <- subs
			}
		}
		return(final)
	}
cDNA_sequin_cleanup <- cleanup(cDNA_sequin)



merge_tails <- function(data, tails) { 
	merged<- merge(data, tails, by.x=c("Read_ID"),by.y=c("read_id") )
	merged$tail_is_valid <- NULL
	merged$read_type <- NULL
	return(merged)
}


cDNA_sequin_tails <- merge_tails(cDNA_sequin_cleanup, cDNA_tails.processed)


cDNA_sequin_tails_final <- cDNA_sequin_tails[,c("Read_ID", "Gene_Name", "Group", "tail_length")]
cDNA_sequin_tails_final$Method <- "EndCapture"






## dRNA Input 
dRNA.input <- read.delim("dRNA_nanopolish.tsv", header=TRUE)
dRNA.input$Group <- gsub("_.*","",dRNA.input$contig)
dRNA <- dRNA.input[,c("read_id","contig","Group", "polya_length" )]
colnames(dRNA) <- c("Read_ID", "Gene_Name", "Group", "tail_length")
dRNA$Method <- "dRNA"

cleanup<- function(data) {
	final<- vector()
		for (gene in unique(data$Gene_Name)) {
			subs <- subset(data, Gene_Name==gene)
			if (nrow(subs) > 30) {
				final <- rbind(final, subs)
			} else {
				subs <- subs
			}
		}
		return(final)
	}
dRNA_sequin_cleanup <- cleanup(dRNA)




both_merged <- rbind(cDNA_sequin_tails_final,dRNA_sequin_cleanup )
both_merged$Unique <- paste(both_merged$Group , both_merged$Method)

      pdf(file="Sequin_cDNA_vs_dRNA_Tails_Boxplot_Min30.pdf",height=4,width=8,onefile=FALSE)
            print(ggplot(both_merged, aes(x = Unique, y = tail_length, fill=Method )) + 
              geom_boxplot(aes(fill = Method),position = position_dodge(0.9), outlier.shape = NA) +
              theme_bw()+
              #coord_cartesian(ylim = c(0, 150))+
              stat_n_text() + 
              geom_hline(yintercept=30, linetype="dashed" ,size=0.3)+
              geom_hline(yintercept=60, linetype="dashed" ,size=0.3)+
              scale_fill_manual(values = c("#91c788", "#ffaaa7"))
              )
         dev.off()

      pdf(file="Sequin_cDNA_vs_dRNA_Tails_Boxplot_AxisFixed_Min30.pdf",height=4,width=8,onefile=FALSE)
            print(ggplot(both_merged, aes(x = Unique, y = tail_length, fill=Method )) + 
              geom_boxplot(aes(fill = Method),position = position_dodge(0.9), outlier.shape = NA) +
              theme_bw()+
              coord_cartesian(ylim = c(0, 150))+
              stat_n_text() + 
              geom_hline(yintercept=30, linetype="dashed" ,size=0.3)+
              geom_hline(yintercept=60, linetype="dashed" ,size=0.3)+
              scale_fill_manual(values = c("#91c788", "#ffaaa7"))
              )
         dev.off()



		both_merged <- both_merged[order(both_merged$Group),]

		both_merged$Gene_Name <- factor(both_merged$Gene_Name, levels=unique(as.character(both_merged$Gene_Name)) )

      pdf(file="Sequin_cDNA_vs_dRNA_Tails_Boxplot_AxisFixed_PerGene_Min30.pdf",height=4,width=20,onefile=FALSE)
            print(ggplot(both_merged, aes(x = Gene_Name, y = tail_length, fill=Method )) + 
              geom_boxplot(aes(fill = Method),position = position_dodge(0.9), outlier.shape = NA) +
              theme_bw()+
              coord_cartesian(ylim = c(0, 150))+
              stat_n_text() + 
              geom_hline(yintercept=30, linetype="dashed" ,size=0.3)+
              geom_hline(yintercept=60, linetype="dashed" ,size=0.3)+
              scale_fill_manual(values = c("#91c788", "#ffaaa7"))
              )
         dev.off()


         #Plot only the paired ones
         Paired_merged <- vector()
         for (gene in unique(both_merged$Gene_Name)){
         	subs <- subset(both_merged, Gene_Name==gene)
         	if (length(unique(subs$Method)) == 2 ) {
         		Paired_merged <- rbind(Paired_merged, subs)
         		} else {}
         	}


      pdf(file="Sequin_cDNA_vs_dRNA_Tails_Boxplot_AxisFixed_PerGene_Min30_PairedOnly.pdf",height=4,width=20,onefile=FALSE)
            print(ggplot(Paired_merged, aes(x = Gene_Name, y = tail_length, fill=Method )) + 
              geom_boxplot(aes(fill = Method),position = position_dodge(0.9), outlier.shape = NA) +
              theme_bw()+
              coord_cartesian(ylim = c(0, 150))+
              stat_n_text() + 
              geom_hline(yintercept=30, linetype="dashed" ,size=0.3)+
              geom_hline(yintercept=60, linetype="dashed" ,size=0.3)+
              scale_fill_manual(values = c("#91c788", "#ffaaa7"))
              )
         dev.off()




Paired_merged$Unique2 <-  paste(Paired_merged$Gene_Name , Paired_merged$Method)
variance_table <- as.data.frame(group_by(Paired_merged, Unique2) %>% 
      summarise(GroupVariance=var(tail_length)))
variance_table$Method <-gsub(".* ","" , variance_table$Unique2)
variance_table$Gene <-gsub(" .*","" , variance_table$Unique2)
variance_table2 <- dcast(variance_table, Gene ~ Method, value.var = "GroupVariance")


      pdf(file="Sequin_cDNA_vs_dRNA_Variance_Scatter.pdf",height=4,width=4,onefile=FALSE)
		print(ggplot(variance_table2, aes(x=log(dRNA), y=log(EndCapture))) + 
		geom_abline(slope=1, intercept=0,linetype="dashed", size=0.2, color= "black")+
		theme_bw()+ 
        ylim(5,8)+
        xlim(5,8)+
		xlab("dRNA log(Variance)")+
        ylab("EndCapture log(Variance)")+
  		geom_point())
  		dev.off()






      pdf(file="Sequin_cDNA_vs_dRNA_Variance_Boxplot.pdf",height=4,width=20,onefile=FALSE)
            print(ggplot(variance_table, aes(x = Method, y = GroupVariance, fill=Method )) + 
              geom_boxplot(aes(fill = Method),position = position_dodge(0.9), outlier.shape = NA) +
              theme_bw()+
              #coord_cartesian(ylim = c(0, 150))+
              stat_n_text() + 
              scale_fill_manual(values = c("#91c788", "#ffaaa7"))
              )
         dev.off()
