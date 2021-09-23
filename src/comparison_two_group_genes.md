### COMPARISON OF TWO GROUP OF GENES

```R


#First Group
  maternal <- read.delim("allTranscripts_riboZero_rnaSeq.maternal.txt", header=FALSE)
  maternal$Group <- "Maternal"
  zyfir <- read.delim("allTranscripts_riboZero_rnaSeq.zyfir.txt", header=FALSE)
  zyfir$Group <- "Zygotic"
  mir430 <- read.delim("allTranscripts_riboZero_rnaSeq.mir430.txt", header=FALSE)
  mir430$Group <- "mir430"
  id_convert <- read.delim("Transcript_ID_To_Gene_ID.txt")
  groups <- rbind(maternal, zyfir, mir430)
  groups_id <- merge(groups,id_convert, by.x="V1", by.y="Transcript.stable.ID")
  groups_id$V1 <- NULL

first_maternal <- subset(groups_id, Group=="Maternal")
first_mir430 <- subset(groups_id, Group=="mir430")
first_zygotic <- subset(groups_id, Group=="Zygotic")



gr_groups <- read.delim("gr_list.txt")
gr_groups2 <- subset(gr_groups, Group != "Rest")


gr_maternal <- subset(gr_groups2, Group=="Maternal")
gr_mir430 <- subset(gr_groups2, Group=="mir430")
gr_zygotic <- subset(gr_groups2, Group=="Zygotic")





 maternal_common <- list(first_maternal[,2],gr_maternal[,1] )
 mir430_common <- list(first_mir430[,2],gr_mir430[,1] )
 zygotic_common <- list(first_zygotic[,2],gr_zygotic[,1] )

library(VennDiagram)


venn_reps <- function(data, label){
venn.diagram(x =data,
     category.names = c("First" , "Second(GR)"),
     filename = paste(label,'First_Second_Group.tiff',sep="_"),
     lwd = 2,
     lty = 'blank',
     fill = c("red", "blue"),
     output=TRUE,
       # Output features
      imagetype="tiff" ,
      height = 500 , 
      width = 500 , 
      resolution = 500,
      compression = "lzw",
     # Numbers
      cex = .3,
      fontface = "bold",
      fontfamily = "sans",
      # Set names
      cat.cex = 0.2,
      cat.fontface = "bold",
      cat.default.pos = "outer",
      #cat.pos = c(-27, 27),
      #cat.dist = c(0.055, 0.055),
      #cat.fontfamily = "sans",
      #rotation = 1
  )
}

venn_reps(maternal_common, "Maternal")
venn_reps(mir430_common, "Mir430")
venn_reps(zygotic_common, "Zygotic")




