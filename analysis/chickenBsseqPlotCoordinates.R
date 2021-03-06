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

# load genes rda
load(file=file.path(procroot,"galGal4_genes.rda")) # objects: genes.comp, genes.gr

# Load bsseq.R file
load(file=file.path(procroot,"bsobject.rda")) # bsobject has bismark,BS.fit.large,BS.fit.small
load(file.path(procroot,"dmrgene.rda"))# dmrtot.gr,genes.gr,prom.gr,ovlmat.prom,ovlmat.body,distnear,proximal,combox,bismark.samp.info

pData <- pData(BS.fit.small)
pData$lab = as.character(paste0(bismark.samp.info$type,bismark.samp.info$time))
pData$col <- as.character(bismark.samp.info$col)

pData(BS.fit.small) <- pData
pData(BS.fit.large) <- pData

combos

# stats analysis
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
if (F){
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


#plotting
if (T){
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
            #test
#            source("~/Code/ilee/util/ilee_plot.R")
 #           pdf(file.path(plotdir,"test.pdf"))
  #          bsseqPlotchicken(dmrtot.gr,genes.gr,prj,dmrlist[[j]],1)
   #         dev.off()
            #end test
            names = c(names,bsseqPlotchicken(dmrtot.gr,genes.gr,prj,dmrlist[[j]],plotnum))
        }
        genenames = cbind(genenames,names)
        dev.off()
    }
    colnames(genenames) = c("rank",comps)
    write.table(genenames,file=file.path(outdir,"methgenes.csv"),quote=FALSE,row.names=F,col.names=T,sep=",")
}

# starting off with dmrtot.gr, genes.gr, and distnear/proxima (hits objects)



