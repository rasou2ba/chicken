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
source("~/Code/timp_genetics/util/timp_seqtools.R")
source("~/Code/timp_genetics/util/read_tools.R")
source("~/Code/ilee/util/ilee_plot.R")

require(doMC)
registerDoMC()
options(cores=4)

# Load genes list
# file is bed format
genelab = c("bin","name","chrom","strand","txStart","txEnd","cdsStart","cdsEnd","exonCount","exonStarts","exonEnds","score","name2","cdsStartStat","cdsEndStat","exonFrames")
Genes.comp= read.table(file=gzfile(genepath),sep="\t",header=F,stringsAsFactors=F,col.names=genelab)
genes = genes.comp[,c(2,13,3,5,6,4)]
colnames(genes) = c("name","ID","chrom","start","end","strand")
names = read.table(file=gzfile(namepath),sep="\t",header=F,stringsAsFactors=F)
genes$ID = names[match(genes$name,names[,1]),2]
genes = na.omit(genes)
genes.gr = GRanges(seqnames=genes$chrom,IRanges(start=genes$start,end=genes$end),strand=genes$strand,id=genes$ID,name=genes$name)

save(list=c("genes.comp","genes.gr"),file=file.path(procroot,"galGal4_genes.rda"))


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
    d.gr = GRanges(seqnames=dmrs[[i]]$chr,IRanges(start=dmrs[[i]]$start,end=dmrs[[i]]$end),meanDiff=dmrs[[i]]$meanDiff,areaStat=dmrs[[i]]$areaStat,direction=dmrs[[i]]$direction,one=rep(combos$one[i],size),two=rep(combos$two[i],size),label=rep(combos$label[i],size))#
    dmrs.gr = append(dmrs.gr,d.gr)
    size = length(blocks[[i]]$chr)
    b.gr = GRanges(seqnames=blocks[[i]]$chr,IRanges(start=blocks[[i]]$start,end=blocks[[i]]$end),meanDiff=blocks[[i]]$meanDiff,areaStat=blocks[[i]]$areaStat,direction=blocks[[i]]$direction,one=rep(combos$one[i],size),two=rep(combos$two[i],size),label=rep(combos$label[i],size))#
    blocks.gr = append(blocks.gr,b.gr)
}

dmrs.gr$type = rep("dmr",length(dmrs.gr))
blocks.gr$type = rep("block",length(blocks.gr))
dmrtot.gr = append(dmrs.gr,blocks.gr)
dmrtot.gr = dmrtot.gr[order(-abs(dmrtot.gr$areaStat))]

## use nearest to do matches
distnear = as.data.frame(distanceToNearest(dmrtot.gr,genes.gr))
proximal = distnear[which(distnear$distance<10000),]

# using promoters fxn to exact match promoters
prom.gr = promoters(genes.gr)
ovlmat.prom = as.data.frame(findOverlaps(dmrtot.gr,prom.gr))
ovlmat.body = as.data.frame(findOverlaps(dmrtot.gr,genes.gr))
length(ovlmat.prom[,1])
length(ovlmat.body[,1])
length(dmrtot.gr)

## save
save(list=c("dmrtot.gr","genes.gr","prom.gr","ovlmat.prom","ovlmat.body","distnear","proximal","combos","bismark.samp.info"),file=file.path(procroot,"dmrgene.rda"))

# stats analysis
load(file.path(procroot,"dmrgene.rda"))
cutoff = length(dmrtot.gr)/10
promcnt = bodycnt=tot=tot.top=bodycnt.top=promcnt.top=numeric()
for (i in unique(dmrtot.gr$label)){
    for (j in unique(dmrtot.gr$type)){
        jdx = which(dmrtot.gr$type==j)
        idx = jdx[which(dmrtot.gr$label[jdx]==i)]
        promidx = unique(match(ovlmat.prom[,1],idx))
        bodyidx = unique(match(ovlmat.body[,1],idx))
        promcnt.top = append(promcnt.top,length(which(idx[promidx]<cutoff)))
        bodycnt.top = append(bodycnt.top,length(which(idx[bodyidx]<cutoff)))
        tot.top = append(tot.top,length(which(idx<cutoff)))
        promcnt = append(promcnt,sum(!is.na(promidx)))
        bodycnt = append(bodycnt,sum(!is.na(bodyidx)))
        tot = append(tot,length(idx))
    }
}
dmrcnt.df = data.frame(comp=rep(unique(dmrtot.gr$label),each=2),type=rep(unique(dmrtot.gr$type),6),promcnt,bodycnt,none=tot-promcnt-bodycnt)
top.df = dmrcnt.df
top.df[,c(3,4,5)]=cbind(promcnt.top,bodycnt.top,tot.top)
top.melt = melt(top.df,measure.var=c(3,4,5))
dmrcnt.melt = melt(dmrcnt.df,measure.var=c(3,4,5))
# top 10%


# stats plotting (counts)
if (TRUE){
    require(ggplot2)
    g = ggplot()+theme_bw()
    g.cnt = g + geom_bar(aes(x=comp,y=value,fill=variable),data=dmrcnt.melt,position="stack",stat="identity")+
        facet_wrap(~type)
    g.top.block = g + geom_bar(aes(x=comp,y=value,fill=variable),data=top.melt[top.melt$type=="block",],position="stack",stat="identity")    
    g.top.dmr = g + geom_bar(aes(x=comp,y=value,fill=variable),data=top.melt[top.melt$type=="dmr",],position="stack",stat="identity")    
    pdf(file.path(plotdir,"dmrstats.pdf"))
    print(g.cnt)
    print(g.top.block)
    print(g.top.dmr)
    dev.off()

}


# plotting
if (TRUE){
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
    genenames = data.frame(rank=rep(seq(plotnum),2))
    for (i in comps){
        pdf(file.path(plotdir,paste0(i,".pdf")),width=8.5,height=11)
        pri = proximal[which(dmrtot.gr[proximal$queryHits]$label==i),]
        names = c()
        for (j in types){
            prj = pri[which(dmrtot.gr[pri$queryHits]$type==j),]
            names = c(names,bsseqPlot(dmrtot.gr,genes.gr,prj,dmrlist[[j]],plotnum))
        }
        genenames = cbind(genenames,names)
        dev.off()
    }
    colnames(genenames) = c("rank",comps)
    write.table(genenames,file=file.path(outdir,"methgenes.csv"),quote=FALSE,row.names=F,col.names=T,sep=",")
}

# starting off with dmrtot.gr, genes.gr, and distnear/proxima (hits objects)



