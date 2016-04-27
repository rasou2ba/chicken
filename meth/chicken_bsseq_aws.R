##First set central code dir for sourcing
codedir="~/Code/timp_genetics"


##Plots go here:
plotdir="~/data/plots"


##All alignment data lives here 
procroot="~/data"

##Load libraries and sources
require(Biostrings)
require(plyr)
require(ggplot2)
require(bsseq)
require(reshape)
require(GenomicRanges)
source("~/data/timp_genetics/util/timp_seqtools.R")
source("~/data/timp_genetics/util/read_tools.R")


require(doMC)
registerDoMC()
options(cores=4)

if (TRUE) {
    bismark.samp.info=data.frame(project="chicken", sample=c("ACTTGA","AGTCAA","AGTTCC","ATGTCA","CAGATC","CCGTCC","CTTGTA","GATCAG","GGCTAC","GTCCGC","GTGAAA","TAGCTT"),
        label=c("E8retina2","E18cornea3","E18cornea4","E18cornea5","E8retina1","E18brain1","E18retina3","E8retina3","E18retina2","E18brain2","E18brain3","E18retina1"),
        quant=NA, type=c("retina", "cornea", "cornea", "cornea", "retina", "brain", "retina", "retina", "retina", "brain", "brain", "retina"), time=c("E8", "E18", "E18", "E18", "E8", "E18", "E18", "E8", "E18", "E18", "E18", "E18"), stringsAsFactors=F) 
    
    bismark.samp.info$filepath=file.path(procroot, bismark.samp.info$sample, paste0(bismark.samp.info$sample, ".cyto.txt"))

    bismark=read.bismark(files=bismark.samp.info$filepath, sampleNames=bismark.samp.info$label)
}


if (TRUE) {
    ##Find CGs to keep
    ##Keep with at least coverage of 2 in at least two of cancer/normal
    BS.cov=getCoverage(bismark, type="Cov", what="perBase")
    colnames(BS.cov)=bismark.samp.info$label
    
    colMeans(BS.cov)
    
                                        #pData(bismark)=bismark.samp.info
    
                                        #from amy
                                        #prolly based on cancer as well for Block finding
    
    BS.fit.large<-BSmooth(bismark,mc.cores=4,verbose=TRUE,ns=500,h=20000)
    
                                        #kasper optimized based on cancer data dont got smaller than this cuz takes forever to smooth to get DMRs
    
    BS.fit.small<-BSmooth(bismark,mc.cores=2,verbose=TRUE,ns=20,h=1000)
    
    save(list=c("bismark", "BS.fit.large", "BS.fit.small"), file="~/data/bsobject.rda")
}

if (TRUE) {

    bismark.samp.info$pheno=paste(bismark.samp.info$type, bismark.samp.info$time, sep=".")

    bismark.samp.info$col=factor(bismark.samp.info$pheno)
    levels(bismark.samp.info$col)=c("blue", "red", "green", "orange")
    
    upheno=unique(bismark.samp.info$pheno)
    
    combos=data.frame(one=upheno[combn(4,2)[1,]], two=upheno[combn(4,2)[2,]])

    combos$label=paste(combos$one, combos$two, sep=".v.")

}
    

if (TRUE) {

    tstat.blocks=list()
    tstat.dmrs=list()
    
    for (i in 1:6)  {

        ##i=1
        
        oney=bismark.samp.info$pheno==combos$one[i]
        twoy=bismark.samp.info$pheno==combos$two[i]
        
        keepCGs=which(rowSums(BS.cov[, oney] >=2) >= 2 & rowSums(BS.cov[, twoy] >= 2) >=2)

        keepy=oney|twoy
        
        BS.fit.large.trim=BS.fit.large[keepCGs,keepy]
        BS.fit.small.trim=BS.fit.small[keepCGs,keepy]
        
        ##T-statistics
        tstat.blocks[[i]]=BSmooth.tstat(BS.fit.large.trim, group1=which(oney[keepy]), group2=which(twoy[keepy]), estimate.var="same", local.correct=F, verbose=T, mc.cores=6)
        tstat.dmrs[[i]]=BSmooth.tstat(BS.fit.small.trim, group1=which(oney[keepy]), group2=which(twoy[keepy]), estimate.var="same", local.correct=T, verbose=T, mc.cores=6)
        
        pdf(file.path(plotdir, paste0(combos$label[i], "_tplot.pdf")))
                                        #plot(bismark.tstat)
        plot(density(getStats(tstat.blocks[[i]])[,5], na.rm=T))
        abline(v=c(-2.75, 2.75))
        plot(density(getStats(tstat.dmrs[[i]])[,5], na.rm=T))
        dev.off()
        
    }
}
    
if (TRUE) {

    blocks=list()
    dmrs=list()
    
    for (i in 1:6) {
    
        blocks[[i]] <- dmrFinder(tstat.blocks[[i]], cutoff = c(-2, 2), stat="tstat", maxGap=1000)  #Amy uses 2 for large block finding
        dmrs[[i]] <- dmrFinder(tstat.dmrs[[i]], qcutoff=c(0.025,0.975), stat="tstat", maxGap=500)
        write.csv(blocks[[i]], file.path(plotdir, paste0(combos$label[i], "_blocks.csv")))
        write.csv(dmrs[[i]], file.path(plotdir, paste0(combos$label[i], "_dmrs.csv")))
    }

    save(file="~/data/processed.rda", list=c("blocks", "dmrs", "tstat.blocks", "tstat.dmrs", "combos", "bismark.samp.info"))
    
}




if (TRUE) {

    load(file="~/data/bsobject.rda")
    load(file="~/data/processed.rda")
    
    pData <- pData(BS.fit.small)
    pData$col <- as.character(bismark.samp.info$col)

    pData(BS.fit.small) <- pData
    pData(BS.fit.large) <- pData
    
    
    for (i in 1:6) {
        
        pdf(file.path(plotdir, paste0(combos$label[i], "_dmrs.pdf")))
        for (j in 1:100) {
            plotRegion(BS.fit.small, dmrs[[i]][j,], extend = 5000, addRegions = dmrs[[i]])
        }
        dev.off()
        
        pdf(file.path(plotdir, paste0(combos$label[i], "_blocks.pdf")))
        for (j in 1:100) {
            plotRegion(BS.fit.large, blocks[[i]][j,], extend = 2e4, addRegions = blocks[[i]])
        }
        dev.off()
        
    }
}


if (TRUE) {

    source("~/data/timp_genetics/util/timp_seqtools.R")

    
    ##DMRS BED
    for (i in 1:6) {

        bsseqdmr2bed(dmrs[[i]], namey=paste0(combos$label[i], "_dmrs"), outdir=plotdir)
        bsseqdmr2bed(blocks[[i]], namey=paste0(combos$label[i], "_blocks"), outdir=plotdir)
        
        
    }

    ##Samples wig
    for (i in 1:12) {
        wig.bsseq(BS.fit.small[,i], filedir=plotdir, modif=bismark.samp.info$label[i], smooth=T)        
    }
    
}
