#!/usr/bin/Rscript
##First set central code dir for sourcing
codedir="~/Code/isac/chicken/meth"

##Plots go here:
outdir="/home/isac/Dropbox/Data/Genetics/MethSeq/170120_chicken"
plotdir=file.path(outdir,"plots")

##All alignment data lives here
datdir="/atium/Data/NGS/Aligned/170120_chicken"

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
source("~/Code/timp_genetics/util/timp_seqtools.R")
source("~/Code/timp_genetics/util/read_tools.R")
source("~/Code/ilee/util/ilee_plot.R")

require(doMC)
registerDoMC()
options(cores=4)

# Load bsseq.R file

load(file=file.path(datdir,"bsobject.rda")) # bsobject has bismark,BS.fit.large,BS.fit.small

#PCA
# get the methylation values: refer to bsseq user guide
totmeth = getMeth(BS.fit.small,type="smooth",what="perBase")
totcov = getCoverage(bismark,type="Cov",what="perBase")
idx = which(rowSums(totcov>=2)==12)
meth = totmeth[idx,]

meth.t = t(meth)

meth.pca = prcomp(meth.t,center=TRUE,retx=TRUE)

meth.x = data.frame(meth.pca$x[,1:6])
meth.x$pheno = pData$pheno
x.pca = t(meth.pca$x)
colnames(x.pca) = pData$pheno
x.melt = melt(x.pca)

require(ggplot2)
g = ggplot(meth.x,aes(x=PC1,y=PC2,colour=pheno))+geom_point()+theme_bw()
g2 = ggplot(meth.x,aes(x=PC3,y=PC4,colour=pheno))+geom_point()+theme_bw()
g3 = ggplot(meth.x,aes(x=PC5,y=PC6,colour=pheno))+geom_point()+theme_bw()

pdf(file.path(plotdir,"pcaBiplot.pdf"))
print(g)
print(g2)
print(g3)
dev.off()

pdf(file.path(plotdir,"pcaVariance.pdf"))
plot(meth.pca,type="l")
dev.off()
