#########
# Library
#########
library(karyoploteR)
library(rtracklayer)

setwd("~/SURFdrive/BleekerLab/11.analysis/karyoplots/20170828_srnadesc_concatenated/")

###########
# Load data
###########

# load genome information for karyotype building
custom.genome <- toGRanges("~/SURFdrive/BleekerLab/05.databases_and_refseq/Solanum_lycopersicum/ITAG3.2/karyoploteR.txt")
custom.cytobands <- toGRanges("~/SURFdrive/BleekerLab/05.databases_and_refseq/Solanum_lycopersicum/ITAG3.2/cytobands.txt") 

# load repet info
#repeats <- import.gff("~/SURFdrive/BleekerLab/05.databases_and_refseq/Solanum_lycopersicum/ITAG3.2/ITAG3.0_RepeatModeler_repeats_light.gff2")

###############
# Filter miRNAs
###############
# Get ShortStack Dicer regions
clusters <- import.gff("ShortStack_D.annotated.gff3",colnames = c("ID","MIRNA"))

# Subset de novo annotated miRNAs
mirnas = clusters[(elementMetadata(clusters)[,"MIRNA"] %in% c("Y"))]
mirnas$ID = unlist(mirnas$ID)

######
# Plot
######

######## Make the default plot: for one chromosome add 'chromosomes=SL3.ch08'
# add space for data.panels
par(oma=c(2,2,2,2))
pp <- getDefaultPlotParams(plot.type=2)
#pp$topmargin = 30
#pp$bottommargin = 30
#pp$data2height = 400
#pp$data1height = 100
pp$data2inmargin = 40

# plot
kp <- plotKaryotype(genome = custom.genome,plot.type = 2,plot.params = pp)
kpAddBaseNumbers(kp, tick.dist = 10000000, tick.len = 10, tick.col="red", cex=1,minor.tick.dist = 1000000, minor.tick.len = 5, minor.tick.col = "gray")
#kpAxis(kp,ymax = 1,side = 2)
kpDataBackground(kp, data.panel = 1,color = "white")
kpDataBackground(kp, data.panel = 2,color = "white")

# plot the data
kpPlotDensity(kp,data=clusters,data.panel = 1,col="#a6cee3",window.size = 1e6)
#kpPlotCoverage(kp,data=clusters,col = "#a6cee3")
kpPlotMarkers(kp, data=mirnas, labels=mirnas$ID,data.panel = 2,text.orientation = "vertical",cex=0.6,adjust.label.position = T)


####################################
# To save all chromosomes one by one
####################################
chroms = as.character(seqnames(custom.genome))
for (i in seq_along(chroms)){
  pdfName = paste("mirgenes",chroms[i],".pdf",sep = "")
  pdf(file = pdfName,width = 5,height = 8)
  kp <- plotKaryotype(genome = custom.genome,plot.type = 2,plot.params = pp,chromosomes = chroms[i])
  kpAddBaseNumbers(kp, tick.dist = 10000000, tick.len = 10, tick.col="red", cex=1,minor.tick.dist = 1000000, minor.tick.len = 5, minor.tick.col = "gray")
  #kpAxis(kp,ymax = 1,side = 2)
  kpDataBackground(kp, data.panel = 1,color = "white")
  kpDataBackground(kp, data.panel = 2,color = "white")
  
  # plot the data
  kpPlotDensity(kp,data=clusters,data.panel = 1,col="#a6cee3",window.size = 1e6)
  #kpPlotCoverage(kp,data=clusters,col = "#a6cee3")
  kpPlotMarkers(kp, data=mirnas, labels=mirnas$ID,data.panel = 2,text.orientation = "vertical",cex=0.6,adjust.label.position = T)
  dev.off()  
}




###############
# From tutorial
###############
# to change plot parameters
pp <- getDefaultPlotParams(plot.type=2)


# add text
#mydata <- toGRanges(data.frame(chr=rep("ch01", 3), start=rep(60e6, 3), end=rep(60e6, 3), y=c(0, 0.5, 1)))
#kpText(kp, data=mydata,  labels=c(0, 0.5, 1))


