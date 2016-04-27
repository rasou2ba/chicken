##First set central code dir for sourcing
codedir="~/Code/isac/chicken/meth"

##Plots go here:
outdir="~/Dropbox/Data/Genetics/MethSeq/150415_chicken/"
plotdir=file.path(outdir,"plots")

##All alignment data lives here
procroot="/mithril/Data/NGS/Aligned/150415_HiSeqChick"

##Gene data
genedir="/mithril/homes/isac/Data/genes"
genefile="galGal4refGene.txt.gz"
genepath=file.path(genedir,genefile)

##Load libraries and sources
require(Biostrings)
require(plyr)
require(ggplot2)
require(bsseq)
require(reshape)
require(GenomicRanges)
library(topGO)
source("~/Code/timp_genetics/util/timp_seqtools.R")
source("~/Code/timp_genetics/util/read_tools.R")
source("~/Code/ilee/util/ilee_plot.R")

require(doMC)
registerDoMC()
options(cores=4)

# Load bsseq.R file
load(file=file.path(procroot,"bsobject.rda")) # bsobject has bismark,BS.fit.large,BS.fit.small
load(file=file.path(procroot,"dmrgene.rda"))# dmrgene has dmrtot.gr,genes.gr,distnear,proximal,combos, and bismark.samp.info

pData <- pData(BS.fit.small)
pData$lab = as.character(paste0(bismark.samp.info$type,bismark.samp.info$time))
pData$col <- as.character(bismark.samp.info$col)

pData(BS.fit.small) <- pData
pData(BS.fit.large) <- pData

combos

# starting off with dmrtot.gr, genes.gr, and distnear/proxima (hits objects)



