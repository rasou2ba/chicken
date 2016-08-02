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
source("~/Code/timp_genetics/util/timp_seqtools.R")
source("~/Code/timp_genetics/util/read_tools.R")
source("~/Code/ilee/util/ilee_plot.R")

require(doMC)
registerDoMC()
options(cores=4)

# Load genes list
# file is bed format
#genes.comp= read.table(file=gzfile(genepath),sep="\t",header=F,stringsAsFactors=F)
#genes = genes.comp[,c(2,13,3,5,6,4)]
#colnames(genes) = c("name","ID","chrom","start","end","strand")
#genes.gr = GRanges(seqnames=genes$chrom,IRanges(start=genes$start,end=genes$end),id=genes$ID,name=genes$name)


# Load bsseq.R file
load(file=file.path(procroot,"bsobject.rda")) # bsobject has bismark,BS.fit.large,BS.fit.small
load(file=file.path(procroot,"processed.rda"))# processed has blocks,dmrs,tstat.blocks,tstat.dmrs,combos,bismark.samp.info

pData <- pData(BS.fit.small)
pData$pheno = bismark.samp.info$pheno
pData$col <- as.character(bismark.samp.info$col)

pData(BS.fit.small) <- pData
pData(BS.fit.large) <- pData

combos

# plotting tstat distrbutions

pdf(file.path(plotdir,"tstatDist.pdf"))
for (i in seq(length(tstat.dmrs))){
    print(i)
    plot(tstat.dmrs[[i]])
}
dev.off()

#PCA?
# get the methylation values: refer to bsseq user guide
totmeth = getMeth(bismark,type="raw",what="perBase")
idx = which(rowSums(!is.nan(totmeth))==12)
meth = totmeth[idx,]
loc = granges(bismark)[idx]

meth.t = t(meth)

meth.pca = prcomp(meth.t,center=TRUE,retx=TRUE)

meth.x = data.frame(meth.pca$x[,1:2])
meth.x$pheno = pData$pheno
x.pca = t(meth.pca$x)
colnames(x.pca) = pData$pheno
x.melt = melt(x.pca)

require(ggplot2)
g = ggplot(meth.x,aes(x=PC1,y=PC2,colour=pheno))+geom_point()+theme_bw()

pdf(file.path(plotdir,"pcaBiplot.pdf"))
print(g)
dev.off()

pdf(file.path(plotdir,"pcaVariance.pdf"))
plot(meth.pca,type="l")
dev.off()
