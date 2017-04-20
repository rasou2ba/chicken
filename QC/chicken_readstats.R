##Plots go here:
outdir="~/Dropbox/Data/Genetics/MethSeq/170120_chicken"
plotdir=file.path(outdir,"plots")
statdir=file.path(outdir,"QC")

##Load libraries and sources
#require(Biostrings)
#require(plyr)
require(ggplot2)
#require(bsseq)
require(reshape)
#require(GenomicRanges)
#source("~/Code/timp_genetics/util/timp_seqtools.R")
#source("~/Code/timp_genetics/util/read_tools.R")
#source("~/Code/ilee/util/ilee_plot.R")

# read in the files
reads = read.table(file=file.path(statdir,"readnum.csv"),sep=",",stringsAsFactors=F)
reads=na.omit(reads)
mapped = read.table(file=file.path(statdir,"alignnum.csv"),sep=",",stringsAsFactors=F)
mapped[,3]=mapped[,2]/2
sampinfo=read.csv(file=file.path(outdir,"infotable.csv"),row.names=1,colClasses="character")

# Process the reads to combine the reads with same barcode
totread=vector()
for (lab in mapped[,1]){
    barcode = substr(lab,1,6)
    idx = grep(barcode,reads[,1])
    readnum = sum(reads[idx,2])
    totread = append(totread,readnum)
}

# create the df
df = data.frame(lab=sampinfo[pmatch(sampinfo[,2],mapped[,1]),3],raw=totread,map=mapped[,3])
df = df[order(df$lab,decreasing=F),]

# output the stats as csv
write.table(x=df,file=file.path(outdir,"readstat.csv"),sep=",",quote=F,row.names=F)

df.melt = melt(df)
df.melt = na.omit(df.melt)

g = ggplot(df.melt,aes(x=lab,y=value))+geom_bar(stat="identity",aes(fill=variable),position="dodge")+theme_bw()

# plot bar graph
pdf(file.path(plotdir,"readstat.pdf"))
print(g)
dev.off()
