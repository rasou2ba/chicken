#!/cm/shared/apps/R/3.2.2/bin/Rscript
##All alignment data lives here 
args = commandArgs(trailingOnly=TRUE)
samp = args[1]
datdir="/scratch/groups/wtimp1/170119_chicken/aligned"
procroot="/scratch/groups/wtimp1/170119_chicken"
##Load libraries and sources
require(bsseq)
library(parallel)

detectCores()
##read in the data
if (T) {
    bismark.samp.info=read.csv(file=file.path(procroot,"infotable.csv"),row.names=1,colClasses="character")
    bismark.samp.info$filepath=file.path(datdir, bismark.samp.info$sample, paste0(bismark.samp.info$sample, ".cyto.txt.gz"))
    ind = which(bismark.samp.info$sample==samp)
    bismark=read.bismark(files=bismark.samp.info$filepath[ind],sampleNames=bismark.samp.info$label[ind],fileType="cytosineReport",verbose=T)
}

## smoothing
if (T) {
    ##from amy
    ##prolly based on cancer as well for Block finding
    BS.fit.large<-BSmooth(bismark,mc.cores=6,parallelBy="chromosome",verbose=TRUE,ns=500,h=20000)
    ##kasper optimized based on cancer data dont got smaller than this cuz takes forever to smooth to get DMRs
    BS.fit.small<-BSmooth(bismark,mc.cores=6,parallelBy="chromosome",verbose=TRUE,ns=20,h=1000)
    save(list=c("bismark.samp.info","bismark", "BS.fit.large", "BS.fit.small"), file=file.path(datdir,paste0("bsobject",samp,".rda")))
}
