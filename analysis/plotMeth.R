##First set central code dir for sourcing
codedir="~/Code/isac/chicken/meth"

##Plots go here:
outdir="~/Dropbox/Data/Genetics/MethSeq/150415_chicken/"
plotdir=file.path(outdir,"plots")

##All alignment data lives here
procroot="/mithril/Data/NGS/Aligned/150415_HiSeqChick"

##Gene data
genedir="/mithril/homes/isac/Data/genes"
genefile="galGal4ensGene.txt.gz"
namefile="galGal4ensemblToGeneName.txt.gz"
genepath=file.path(genedir,genefile)
namepath=file.path(genedir,namefile)

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

# Load genes list
# file is bed format
genes.comp= read.table(file=gzfile(genepath),sep="\t",header=F,stringsAsFactors=F)
genes = genes.comp[,c(2,13,3,5,6,4)]
colnames(genes) = c("name","ID","chrom","start","end","strand")
names = read.table(file=gzfile(namepath),sep="\t",header=F,stringsAsFactors=F)
genes$ID = names[match(genes$name,names[,1]),2]
genes = na.omit(genes)
genes.gr = GRanges(seqnames=genes$chrom,IRanges(start=genes$start,end=genes$end),id=genes$ID,name=genes$name)

# Load bsseq.R file
load(file=file.path(procroot,"bsobject.rda")) # bsobject has bismark,BS.fit.large,BS.fit.small
load(file=file.path(procroot,"processed.rda"))# processed has blocks,dmrs,tstat.blocks,tstat.dmrs,combos,bismark.samp.info

pData <- pData(BS.fit.small)
pData$lab = as.character(paste0(bismark.samp.info$type,bismark.samp.info$time))
pData$col <- as.character(bismark.samp.info$col)

pData(BS.fit.small) <- pData
pData(BS.fit.large) <- pData

combos

## load candidates file
candidates = read.table(file=file.path(outdir,"candidateGenes.csv"),sep=",",header=T,stringsAsFactors=F)

require(reshape2)
require(ggplot2)
width = 10000
cplot = ggplot() + theme_bw()
for (i in seq(length(candidates[1,]))){
    names = candidates[,i]
    genes.cand = genes.gr[na.omit(pmatch(toupper(names),toupper(genes.gr$id)))]
    pdf(file.path(plotdir,paste0(colnames(candidates)[i],".pdf")))
    for (j in seq(length(genes.cand))){
        nwid = width+width(genes.cand[j])
        genestart = start(genes.cand[j])
        geneend = end(genes.cand[j])
        region = resize(genes.cand[j],width=nwid,fix="center")
        regy = granges(BS.fit.small) %within% region
        meth = getMeth(BS.fit.small[regy,],type="smooth",what="perBase")
        cov = as.data.frame(getCoverage(BS.fit.small[regy,],type="Cov",what="perBase"))
        colnames(meth)=colnames(cov)=sampleNames(BS.fit.small)
        cov$loc = as.data.frame(granges(BS.fit.small[regy]))[,2]
        meth.melt= melt(meth)
        cov.melt = melt(cov,id.vars="loc")
        dat = data.frame(loc=cov.melt$loc,sample=cov.melt$variable,meth=meth.melt$value,cov=cov.melt$value)
        dat$lab = pData(BS.fit.small)$lab[match(dat$sample,rownames(pData(BS.fit.small)))]
        
        qmeth = cplot + stat_smooth(data=dat,aes(x=loc,y=meth,group=sample,color=lab),se=F,method="loess")+
            geom_rug(data=dat,aes(x=loc))+scale_y_continuous(limits=c(0,1))+
            geom_rect(mapping=aes(xmin=genestart,xmax=geneend,ymin=0,ymax=1,fill="gene"),alpha=0.2)+ labs(title=genes.cand[j]$id,x="",y="methylation")
        qcov = cplot + stat_smooth(data=dat,aes(x=loc,y=cov,group=sample,color=lab),se=F,method="loess")+
            labs(x=seqnames(genes.cand[j]),y="coverage")
        multiplot(qmeth,qcov,cols=1)
    }
    dev.off()
}
