##################################################
### TAILFINDR COMPARISON OF TS vs ANNEALING CC ###
###################################################
########## OGUZHAN BEGIK APRIL 2020 ###############

library(stringr)
library(ggplot2)
library(ggbeeswarm)
library(EnvStats)

ann.tails <- read.delim("Annealing_tails.csv", sep=",")
cap.tails <- read.delim("Capture_tails.csv", sep=",")

ann.input <- read.delim("Annealing.read_id", header=FALSE)
cap.input <- read.delim("Capture.read_id", header=FALSE)


 manipulate_read_id <- function(data) { 
    #Rename the columns
    colnames(data) <-  c("Read_ID", "Gene_Name", "Quality")
    #Quality filter
    data_filt <- subset(data, Quality > 59)
    #Seperate the column of Gene Name into three
    column <- as.data.frame(str_split_fixed(data_filt$Gene_Name, n=3, pattern="_"))
    #Paste the info as column
    data_filt$Group <- column$V1
    return(data_filt)
 }

ann <- manipulate_read_id(ann.input)
cap <- manipulate_read_id(cap.input)

manipulate_tail <- function(data) { 
    data2 <- data[,c("read_id", "read_type", "tail_is_valid", "tail_length")]
    data2$tail_length <- as.numeric(as.character(data2$tail_length))
    data_filt <- subset(data2, tail_is_valid=="TRUE")
    data_filt_T <- subset(data_filt, read_type=="polyT")
    return(data_filt_T)
}

ann.tails_processed <- manipulate_tail(ann.tails)
cap.tails_processed <- manipulate_tail(cap.tails)


merge_tails <- function(read_data, tail_data, label) { 
    merged <- merge(read_data, tail_data, by.y=c("read_id"), by.x=c("Read_ID"))
    merged[["tail_length"]][is.na(merged[["tail_length"]])] <- 0
    merged$Sample <- label
    return(merged)
}

ann.merged <- merge_tails(ann, ann.tails_processed, "Annealing")
cap.merged <- merge_tails(cap, cap.tails_processed, "Capture")





both_merged <- rbind(ann.merged, cap.merged)
both_merged$Unique <- paste(both_merged$Gene_Name, both_merged$Sample, sep="_")




        pdf(file="AnnealingVSCapture_Tails.pdf",height=4,width=10,onefile=FALSE)
            print(ggplot(both_merged, aes(x=Unique, y=tail_length)) + 
                geom_quasirandom(varwidth = TRUE, aes(color=Group))+
                geom_boxplot(aes(alpha=0), outlier.shape=NA)+
                stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                    geom = "crossbar", width = 0.7, color="#c06c84")+
                theme_bw()+
                xlab("Group")+
                ylab("Tail length") +
                theme(axis.text=element_text(size=14),strip.text = element_text(size=13),
                    axis.title=element_text(size=17,face="bold"),
                    legend.title = element_text(size = 20),
                    legend.text = element_text(color = "black", size=15)))
        dev.off()


        pdf(file="AnnealingVSCapture_Tails_Boxplot.pdf",height=4,width=7,onefile=FALSE)
            print(ggplot(both_merged, aes(x = Unique, y = tail_length )) + 
              geom_boxplot(aes(fill = Sample),position = position_dodge(0.9), outlier.shape=NA ) +
              stat_n_text() + 
              coord_cartesian(ylim = c(0,175))+
              theme_bw()+
              scale_fill_manual(values = c("#999999", "#E69F00"))
              )
         dev.off()

            