#!/usr/bin/env Rscript
#####################################################################
### Merge ShortStack miRNA results
# author: marc galland
# contact: m.galland@uva.nl
# version 1.0 created on Jan 12, 2017

## Script to merge ShortStack results output
## If n result tables, then they will be merged in one final table 
## How? By locus (coordinates on the pangenome) 

# Usage
# Rscript --vanilla merge_shortstack_results.R [path/to/shortStack_directory] [path/to/output_folder]

# Arguments: 
#   shortstack_directory: directory that contains one folder per tomato accession 
#   output_folder: a path to a folder that will contain all results

#######################################################################
# working directory
setwd("~/SURFdrive/trichome_team/06.results_from_xp/sRNA-Seq/20170113_srnadesc/")

# library
library(dplyr)

# list directories containing the results.txt file
accessions = list.dirs("shortstack/",full.names = F,recursive = F)
result.dirs = sapply(X = accessions,FUN = function(x) {file.path(getwd(),"shortstack",x,"results.txt")})

# get all results file
dfs = lapply(X = result.dirs,FUN = function(x) {read.delim(x,header = T,stringsAsFactors = F)})

# select only miRNAs
# rename column (add accession)
for (i in seq_along(dfs)){
  dfs[[i]] = filter(dfs[[i]],MIRNA == "Y")
  colNames = colnames(dfs[[i]])[2:ncol(dfs[[i]])]
  sampleName = names(dfs)[i]
  newColNames = sapply(colNames,FUN = function(x){paste(x,sampleName,sep = ".")})
  colnames(dfs[[i]])[2:ncol(dfs[[i]])] = newColNames
}

# reduce (recursively performs a full outer merge on the dataframes)
byLocus = Reduce(function(...) merge(...,by="X.Locus",all=T),dfs)
write.table(byLocus,file = "all_merged.txt",quote = F,sep="\t",row.names = F)

# same merge but with less column (only miRNA)
# rename column (add accession)
for (i in seq_along(dfs)){
  dfs[[i]] = dfs[[i]][c(1,10)]
}
byLocus.short = Reduce(function(...) merge(...,by="X.Locus",all=T),dfs)
write.table(byLocus.short,file = "all_merged.shortened.txt",quote = F,sep="\t",row.names = F)


