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
    dat=read.delim(file=file.path(datdir,"meanqual.txt"),header=T,colClasses="character")
    pos=as.numeric(sapply(strsplit(dat[,1],"-"),"[",1))
    dat[,2]=as.numeric(dat[,2])
    dat[,3]=as.numeric(dat[,3])
    dat$base=pos
    df = melt(dat,measure.vars=c("read1","read2"))    
    
    q = ggplot(df,aes(x=base,y=value,group=variable,color=variable))+ylim(0,38)+
	geom_hline(yintercept=28,color="orange",alpha=0.7,linetype=2,show.legend=TRUE,size=0.5)+
        labs(x="Sequence Position",y="Mean Quality Score")+
        geom_line(size=1)+theme_bw()
    pdf(file.path(plotdir,"qualmean.pdf"),width=8,height=6)
    print(q)
    dev.off()
}
