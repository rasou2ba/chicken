#!/usr/bin/Rscript
##All alignment data lives here 
datdir="/atium/Data/NGS/Aligned/170120_chicken"

##root for processing
procroot="/home/isac/Dropbox/Data/Genetics/MethSeq/170120_chicken"
##Plots go here:
plotdir=file.path(procroot,"plots")

##Load libraries and sources
require(ggplot2)
require(reshape)

##read in the data
if (TRUE) {
    bismark.samp.info=read.csv(file=file.path(procroot,"infotable.csv"),row.names=1,colClasses="character")
    bismark.samp.info$filepath=file.path(datdir, bismark.samp.info$sample, paste0(bismark.samp.info$sample, ".full.M-bias.txt"))
    
    for (ind in 1:12){
        samp=bismark.samp.info$label[ind]
        mbias=read.delim(file=bismark.samp.info$filepath[ind],header=F,skip=3,colClasses="character")
        cpg1=mbias[1:121,]
        cpg2=mbias[373:(373+119),]
        cpg=data.frame(sapply(rbind(cpg1,cpg2),function(x) as.numeric(x)))
        cpg=cbind(cpg,c(rep(times=121,"read 1"),rep(times=120,"read2")))
        colnames(cpg)=c("pos","countM","countUM","Mperc","cov","read")
        q = ggplot(cpg,aes(x=pos,y=Mperc,group=read,color=read))+ ylim(0,100) +
            labs(title=samp,x="Sequence Position",y="Methylation Percentage")+
            geom_line()+theme_bw()
        
        if (T){
            pdf(file.path(plotdir,paste0(samp,"_mbias.pdf")))
            print(q)
            dev.off()
        }
        
    }
}
