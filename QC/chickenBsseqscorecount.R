#!/usr/bin/Rscript
datdir="/home/isac/Dropbox/Data/Genetics/MethSeq/170120_chicken/QC"

##root for processing
procroot="/home/isac/Dropbox/Data/Genetics/MethSeq/170120_chicken"
##Plots go here:
plotdir=file.path(procroot,"plots")

##Load libraries and sources
require(ggplot2)
require(reshape)

##read in the data
if (TRUE) {
    dat=read.delim(file=file.path(datdir,"qualitycount.txt"),header=T,colClasses="numeric")
    #dat[,2]=as.numeric(dat[,2])
    #dat[,3]=as.numeric(dat[,3])
    #dat$base=pos
    df = melt(dat,measure.vars=c("read1","read2"))    
    
    q = ggplot(df,aes(x=score,y=value,group=variable,color=variable))+
	geom_vline(xintercept=28,color="orange",alpha=0.7,linetype=2,show.legend=TRUE,size=0.5)+
        labs(x="Score",y="Average Count")+
        geom_line(size=1)+theme_bw()
    pdf(file.path(plotdir,"scorecount.pdf"),width=8,height=6)
    print(q)
    dev.off()
}
