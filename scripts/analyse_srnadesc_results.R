#!/usr/bin/env Rscript
#####################################################################
### Analyse small RNA results from the srnadesc pipeline (various stats (counts,lengths) and ShortStack results 
# author: marc galland
# contact: m.galland@uva.nl
# Location of the srnadesc pipeline: https://github.com/BleekerLab/Solanum_sRNAs/tree/master/snakepipelines/srnadesc (version v1.0, commit #3979551)


# Usage
# Rscript --vanilla analyse_shortstack_results.R -s [path/to/shortStack_directory] -o [path/to/output_folder] -a [path/to/sample/annotation]

# Arguments: 
#   srnadesc: ShortStack result directory that contains one folder per sample (e.g. tomato accession) 
#   output_folder: a path to a folder to store all plots and other analyses

# Example of ShortStack directory structure
# shortstack/
#   sample1/
#     Results.txt
#     Counts.txt
#     Unplaced.txt
#      ...
#   sample2/
#      Results.txt
#      Counts.txt
#      Unplaced.txt
#      ...

#######################################################################
# accession to species
#accession2species = read.delim("../../../07.data/accession2species.txt",header = T,stringsAsFactors = F)


######################################
# Libraries and command-line arguments
######################################

# load library
library(RColorBrewer)
library(reshape2)
library(ggplot2)
library(optparse)
library(tidyverse,quietly = T,warn.conflicts = F)

# parse command line arguments
option_list = list(
  make_option(c("-s", "--shortstack"), type="character",action="store",help="Shortstack result directory",metavar="character"),
  make_option(c("-o", "--outdir"), type="character",default="results/",help="output directory for results"),
  make_option(c("-a","--sample2annotation"),type="character",default=NULL,help="optional sample annotation file")
) 
opt_parser = OptionParser(option_list=option_list);
args = parse_args(opt_parser)

# print arguments
print("These are the arguments you specified")
print("#####################################")
print(paste("Shortstack result directory",args$shortstack,sep = ":"))
print(paste("Where to store the results",args$outdir,sep = ":"))
print(paste("Sample annotation file",args$sample2annotation,sep = ":"))
print("#####################################")

# creates result directory
dir.create(args$outdir,showWarnings = F)

############
## Load data
############
accessions = list.dirs(args$shortstack,full.names = F,recursive = F)

#### Size distributions of small RNAs
distri.originals = list.files(path="distri",pattern = "*original*",full.names = T)
distri.trimmed = list.files(path="distri",patter = "*trimmed*",full.names = T)

distri.originals = lapply(distri.originals,FUN = function(x){read.delim(x,header=T,stringsAsFactors = F)})
distri.trimmed = lapply(distri.trimmed,FUN = function(x){read.delim(x,header=T,stringsAsFactors = F)})

##### ShortStack
# list directories containing the ShortStack results.txt file
shortstack.dirs = sapply(X = accessions,FUN = function(x) {file.path(args$shortstack,x,"results.txt")})


# get all ShortStack results file
shortstacks = map(shortstack.dirs,function(x) {read.delim(x,header = T,stringsAsFactors = F)})

#shortstacks = lapply(X = shortstack.dirs,FUN = function(x) {read.delim(x,header = T,stringsAsFactors = F)})
for (i in seq_along(shortstacks)){
  colnames(shortstacks[[i]])[1] = "locus"
}


########################################
## small RNA length distribution
########################################
all.sizes = list(distri.originals,distri.trimmed)

# proportion
for (j in seq_along(all.sizes)){
  for (i in seq_along(all.sizes[[j]])){
    df = all.sizes[[j]][[i]]
    df = as.data.frame(table(df))
    totalReads = sum(df$Freq)
    df$percentages = df$Freq / totalReads * 100
    colnames(df)=c("length","freq",accessions[i])
    all.sizes[[j]][[i]] = df[,c("length",accessions[i])]
  }
}

for (j in seq_along(all.sizes)){
  merged = Reduce(function(...) merge(...,by="length",all.x=T),all.sizes[[j]])
  m_merged = melt(merged,id.vars = "length",variable.name = "accession")
  m_merged = left_join(m_merged,accession2species,by="accession")
  m_merged = m_merged[order(m_merged$species),]
  all.sizes[[j]] = m_merged
}

########################################
## Statistics on all small RNAs clusters
########################################

#### how many sRNA clusters per accession?
clusters = shortstacks
# rename column (add accession)
for (i in seq_along(clusters)){
  accessionName = names(clusters)[i]
  clusters[[i]] = melt(clusters[[i]],id.vars = "locus",variable.name = "param",value.name = "value")
  clusters[[i]]$accession = accessionName
}

# row bind all of these dfs
all.clusters = plyr::ldply(clusters,data.frame)
all.clusters = left_join(all.clusters,accession2species,by="accession")

nClusters = all.clusters %>%
  group_by(accession,species) %>%
  summarise(nb = length(locus))

nClusters <- data.frame(lapply(nClusters,as.character),stringsAsFactors = F)
nClusters = nClusters[order(nClusters$species),]
nClusters$nb = as.numeric(nClusters$nb)
nClusters$accession = factor(nClusters$accession,levels = nClusters$accession)

#### What are the fraction of clusters per length 
dicerCalls = filter(all.clusters,param=="DicerCall")
dicerCalls.tables = list()
accessions = unique(dicerCalls$accession)
for (i in seq_along(accessions)){
  df = filter(dicerCalls,accession == accessions[i])
  df = as.data.frame(table(df$value))
  colnames(df) = c("length","freq")
  df$percent = df$freq / sum(df$freq) *100
  df$accession = accessions[i]
  dicerCalls.tables[[i]] = df
}
all.dicerCalls = plyr::ldply(dicerCalls.tables,data.frame)

dicerCalls %>%
  group_by(accession,species) %>%
  summarise(l = length(value))

#### how many phased sRNA clusters per accession?


#### total abundance (RPM) of sRNA clusters per DicerCall (I take the total number of reads for the locus)
all.dfs = plyr::ldply(shortstacks,data.frame,.id = "accession")
m_shortstacks = melt(data = all.dfs,id.vars=c("accession","DicerCall"),measure.vars = "Reads",value.name = "counts")

originalReadsFiles = list.files("counts",pattern = "*01_original*",full.names = T)
originalReadsFiles = lapply(originalReadsFiles,FUN = function(x){read.delim(x,header = F,sep = "\t")})
originalReadsFiles = plyr::ldply(originalReadsFiles,data.frame)
colnames(originalReadsFiles) = c("accession","step","nb")

rpmPerDicerCall = m_shortstacks %>%
  group_by(DicerCall,accession) %>%
  summarise(sum = sum(counts))
rpmPerDicerCall = left_join(rpmPerDicerCall,originalReadsFiles,by="accession")
rpmPerDicerCall = left_join(rpmPerDicerCall,accession2species,by="accession")
rpmPerDicerCall = mutate(rpmPerDicerCall,rpm = sum / nb * 10^5)

# for plotting (reordering levels)
rpmPerDicerCall = rpmPerDicerCall[order(rpmPerDicerCall$species),]
rpmPerDicerCall$accession = factor(rpmPerDicerCall$accession,levels = unique(rpmPerDicerCall$accession))


#########################
## miRNA analysis
########################

##### select only miRNAs and merge by locus
byLocus = shortstacks
# rename column (add accession)
for (i in seq_along(byLocus)){
  byLocus[[i]] = filter(byLocus[[i]],MIRNA == "Y")
  colNames = colnames(byLocus[[i]])[2:ncol(byLocus[[i]])]
  sampleName = names(byLocus)[i]
  newColNames = sapply(colNames,FUN = function(x){paste(x,sampleName,sep = ".")})
  colnames(byLocus[[i]])[2:ncol(byLocus[[i]])] = newColNames
}
# reduce (recursively performs a full outer merge on the dataframes)
byLocus = Reduce(function(...) merge(...,by="locus",all=T),byLocus)
write.table(byLocus,file = "all_merged.byLocus.txt",quote = F,sep="\t",row.names = F)

# same merge but with less column (only miRNA)
# rename column (add accession)
byLocus = shortstacks
for (i in seq_along(byLocus)){
  byLocus[[i]] = filter(byLocus[[i]],MIRNA == "Y")
  colNames = colnames(byLocus[[i]])[2:ncol(byLocus[[i]])]
  sampleName = names(byLocus)[i]
  newColNames = sapply(colNames,FUN = function(x){paste(x,sampleName,sep = ".")})
  colnames(byLocus[[i]])[2:ncol(byLocus[[i]])] = newColNames
  byLocus[[i]] = byLocus[[i]][c(1,10)]
}
byLocus.short = Reduce(function(...) merge(...,by="locus",all=T),byLocus)
write.table(byLocus.short,file = "all_merged.byLocus.shortened.txt",quote = F,sep="\t",row.names = F)

#### Merge by mature miRNA sequence
byMajorRNA = shortstacks
for (i in seq_along(byMajorRNA)){
  byMajorRNA[[i]] = filter(byMajorRNA[[i]],MIRNA == "Y")
  byMajorRNA[[i]] = byMajorRNA[[i]][c(10)]
}
# merge all and then extracts only unique miRNAs
byMajorRNA.all = unique(Reduce(function(...) merge(...,all=TRUE),byMajorRNA))
byMajorRNA.all$length = sapply(byMajorRNA.all$MajorRNA,nchar)
write.table(byMajorRNA.all,file = file.path(outdir,"all_merged.byMajorRNA.txt"),quote = F,sep="\t",row.names=F)


#####################################################
## Plots of all analyses
#####################################################

my_theme <- theme_grey() +
  theme(
    axis.text = element_text(size = 8,face = "plain")
  )

### Size distribution of oringinal / trimmed reads
names(all.sizes)=c("original","trimmed")

for (j in seq_along(all.sizes)){
  if (names(all.sizes[j]) == "original"){
    df = all.sizes[[j]]
    # For x axis
    maxLength = max(as.integer(df$length))
    tickForPlot = 50 # a tick every ... nucleotides
    upperLimit = maxLength %/% tickForPlot + 1 # quotient (non-float reminder of division) and adds one
    g <- ggplot(data=df,aes(x=length,y = value,fill=species)) +
      geom_bar(stat="identity") +
      facet_wrap(~ accession) +
      labs(x = "Accession",y = "% of total small RNAs") +
      scale_fill_brewer(palette = "Set3") +
      my_theme
    ggsave(filename = file.path(outdir,paste("size_",names(all.sizes[j]),".png",sep = "")),plot = g,width = 7,height = 5,dpi = 600)
    ggsave(filename = file.path(outdir,paste("size_",names(all.sizes[j]),".svg",sep = "")),plot = g,width = 7,height = 5)
  } else {
    df = all.sizes[[j]]
    g2 <- ggplot(data=df,aes(x=length,y = value,fill=species)) +
      geom_bar(stat="identity",colour="black") +
      facet_wrap(~ accession) +
      labs(x = "Accession",y = "% of trimmed small RNAs") +
      scale_fill_brewer(palette = "Set3") +
      my_theme
    ggsave(filename = file.path(outdir,paste("size_",names(all.sizes[j]),".png",sep = "")),plot = g2,width = 7,height = 5,dpi = 400)
    ggsave(filename = file.path(outdir,paste("size_",names(all.sizes[j]),".svg",sep = "")),plot = g2,width = 7,height = 5)
  }
}


#### Dicer call fraction per accession
p <- ggplot(data = all.dicerCalls,aes(x = accession,y = percent,fill=length)) +
  geom_bar(stat="identity",colour="black") +
  labs(x = "Accession",y = "Fraction of total small RNA clusters") +
  theme(axis.text.x = element_text(angle=0,hjust=0.5)) +
  my_theme
print(p)  
ggsave(filename = file.path(outdir,"fraction_dicer_call.svg"),plot = p,width = 7,height = 5)
ggsave(filename = file.path(outdir,"fraction_dicer_call.png"),plot = p,width = 7,height = 5,dpi = 400)

#### Number of small RNA clusters per accession
p1 <- ggplot(data = nClusters,aes(x = accession,y = nb,fill=species)) +
  geom_bar(stat="identity",colour="black") +
  geom_text(aes(label=nb),vjust=-0.5,size=4) +
  scale_fill_brewer(palette = "Set3") +
  labs(x = "Accession",y = "Number of small RNA clusters") +
  my_theme
print(p1)  
ggsave(filename = file.path(outdir,"number_of_small_rnas_per_accession.svg"),plot = p1,width = 7,height = 5)
ggsave(filename = file.path(outdir,"number_of_small_rnas_per_accession.png"),plot = p1,width = 7,height = 5,dpi=400)


#### Cluster abundance per DicerCall 
p2 <- ggplot(data = rpmPerDicerCall,aes(x = DicerCall,y = rpm,fill=species)) +
  geom_bar(stat="identity",colour="black") +
  facet_wrap(~ accession,nrow = 4) +
  scale_fill_brewer(palette = "Set3",guide = guide_legend(title = "Species")) +
  labs(x = "Length of small RNA cluster ('Dicer call')",y = expression("small RNA cluster abundance (RPM x 10^5)")) +
  scale_y_continuous(limits=c(0,30000)) +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.x = element_text(angle=0,hjust=0.5)) +
  my_theme
print(p2)
ggsave(filename = file.path(outdir,"cluster_abundance_per_dicercall.svg"),plot = p2,width = 7,height = 5)
ggsave(filename = file.path(outdir,"cluster_abundance_per_dicercall.png"),plot = p2,width = 7,height = 5,dpi=400)








