#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

####
# This script will make stack bar graphs based on the read counts from multiple samples
# There should be one tabulated data file sample per folder and one file containing the counts in that folder.
# Table has to be like
#   sample    step          number
# ----------------------------------
#   PI127826	01.original	  4771609
#   PI127826	02.spikes	    4682739
#   PI127826	03.trimmed	  3592756
#   PI127826	04.virus	    3583722
#   PI127826	05.ncRNA	    3022151
#   PI127826	06_shortstack 1474378
#   LA1777	  01.original	  3022151
#   LA1777  	02.spikes	    2583722
#   LA1777  	03.trimmed	  1583711   
#   etc.
#

#######
# Usage
#######
# Rscript --vanilla size.R [input file] [output directory]

###########
# libraries
###########
# install if not already installed
if (is.element('optparse', installed.packages()[,1]))
{
} else
{
  install.packages('optparse')
}
if (is.element('ggplot2', installed.packages()[,1]))
{
} else
{
  install.packages('dplyr')
}

# load library
library(ggplot2)
library(dplyr)
library(stats)
library(optparse)


###########
# Arguments
###########
## Default setting when no arguments passed
if(length(args) < 1) {
  args <- c("--help")
}

option_list <- list(
  make_option(c("-i", "--inputdir"), default="./",help="directory with count files (one sample per file"),
  make_option(c("-o", "--outdir"), default="./results/",help="desired location for results")
  )

opt <- parse_args(OptionParser(option_list=option_list))

########### 
# Read data
###########
countList = list.files(opt$inputdir)
countFiles = lapply(countList,FUN = function(x) {read.delim(x,header = F,sep="\t",stringsAsFactors = F,colClasses = c("character","character","numeric"))})  
for(i in seq_along(countFiles)){
  colnames(countFiles[[i]]) = c("sample","step","nb")
}
names(countFiles) = list.files(opt$inputdir)

##################################
## Calculate percentages and merge
##################################
Percentages = sapply(countFiles,FUN = function(df) { df$nb / max(df$nb)})

for(i in seq_along(countFiles)){
  countFiles[[i]]$pct = Percentages[,i]
} 

merged.df = Reduce(function(...) merge(...,all=T), countFiles)

############################
# plot stacked bar histogram
############################
# The palette with grey:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

g <- ggplot(merged.df,mapping = aes(x = step,y=pct,fill=sample)) +
  facet_grid(. ~ sample) +
  geom_bar(stat = "identity",colour="black") +
  theme(axis.text.x = element_text(angle = 40,hjust = 1)) +
  # To use for fills, add
  scale_fill_manual(values=cbPalette) +
  labs(x = "Analysis step",y = "Percentage of sequences mapped (%)")
print(g)



