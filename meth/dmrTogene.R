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

# Load genes list
# file is bed format
genes.comp= read.table(file=gzfile(genepath),sep="\t",header=F,stringsAsFactors=F)
genes = genes.comp[,c(2,13,3,5,6,4)]
colnames(genes) = c("name","ID","chrom","start","end","strand")
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
# Create GRanges objects from the DMRs
dmrs.gr = GRanges()
blocks.gr = GRanges()
for (i in seq(length(combos[,1]))){
    size=length(dmrs[[i]]$chr)
    d.gr = GRanges(seqnames=dmrs[[i]]$chr,IRanges(start=dmrs[[i]]$start,end=dmrs[[i]]$end),stat=dmrs[[i]]$areaStat,direction=dmrs[[i]]$direction,one=rep(combos$one[i],size),two=rep(combos$two[i],size),label=rep(combos$label[i],size))#
    dmrs.gr = append(dmrs.gr,d.gr)
    size = length(blocks[[i]]$chr)
    b.gr = GRanges(seqnames=blocks[[i]]$chr,IRanges(start=blocks[[i]]$start,end=blocks[[i]]$end),stat=blocks[[i]]$areaStat,direction=blocks[[i]]$direction,one=rep(combos$one[i],size),two=rep(combos$two[i],size),label=rep(combos$label[i],size))#
    blocks.gr = append(blocks.gr,b.gr)
}

dmrs.gr$type = rep("dmr",length(dmrs.gr))
blocks.gr$type = rep("block",length(blocks.gr))
dmrtot.gr = append(dmrs.gr,blocks.gr)
dmrtot.gr = dmrtot.gr[order(-abs(dmrtot.gr$stat))]

## use nearest to do matches
distnear = as.data.frame(distanceToNearest(dmrtot.gr,genes.gr))
proximal = distnear[which(distnear$distance<10000),]

save(list=c("dmrtot.gr","genes.gr","distnear","proximal","combos","bismark.samp.info"),file=file.path(procroot,"dmrgene.rda"))


# plotting
if (FALSE){
    bins= c(seq(0,2.6e6,10000))
    distribution = hist(distnear$distance,breaks=bins,plot=FALSE)
    dist.df = data.frame(distance=distribution$breaks[-length(distribution$breaks)],freq=distribution$density,count=distribution$counts)

    dist0 = hist(distnear$distance[which(distnear$distance!=0)],breaks=bins,plot=FALSE)
    dist0.df = data.frame(distance=dist0$breaks[-length(dist0$breaks)],freq=dist0$density,count=dist0$counts)
    
    distsub = dist.df[dist.df$distance<1e6,]
    dist0sub = dist0.df[dist0.df$distance<1e6,]
    
    c = ggplot()+theme_bw()+theme(legend.position="none") 
    q = c + geom_bar(data=distsub,aes(x=distance,y=freq),stat="identity",fill="white",colour="blue")+
        labs(title="distance distribution")
    qnozero = c + geom_bar(data=dist0sub,aes(x=distance,y=freq),stat="identity",fill="white",colour="blue")+
        labs(title="Excluding distance = 0")

    pdf(file.path(plotdir,"distance_distributions.pdf"),width=8.5,height=11)
    print(q)
    print(qnozero)
    dev.off() 
    # plot methylation/coverage for top dmrs
    types = c("dmr","block")
    comps = combos$label
    dmrlist = list()
    dmrlist[["dmr"]]=BS.fit.small
    dmrlist[["block"]]=BS.fit.large
    plotnum = 10
    
    for (i in comps){
        pdf(file.path(plotdir,paste0(i,".pdf")),width=8.5,height=11)
        pri = proximal[which(dmrtot.gr[proximal$queryHits]$label==i),]
        for (j in types){
            prj = pri[which(dmrtot.gr[pri$queryHits]$type==j),]
            q = bsseqPlot(dmrtot.gr,genes.gr,prj,dmrlist[[j]],plotnum)
        }
        dev.off()
    }
    
}

# starting off with dmrtot.gr, genes.gr, and distnear/proxima (hits objects)



