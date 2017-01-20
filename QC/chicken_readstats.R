##Plots go here:
outdir="~/Dropbox/Data/Genetics/MethSeq/150415_chicken"
plotdir=file.path(outdir,"plots")
statdir=file.path(outdir,"QC")

##All alignment data lives here
procroot="/mithril/Data/NGS/Aligned/150415_HiSeqChick"

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

# Process the reads to combine the reads with same barcode
totread=vector()
for (lab in mapped[,1]){
    barcode = substr(lab,1,6)
    idx = grep(barcode,reads[,1])
    readnum = sum(reads[idx,2])
    totread = append(totread,readnum)
}

# sample information manually inputted
sampinfo= data.frame(sample=c("ACTTGA","AGTCAA","AGTTCC","ATGTCA","CAGATC","CCGTCC","CTTGTA","GATCAG","GGCTAC","GTCCGC","GTGAAA","TAGCTT"),label=c("E8retina2","E18cornea3","E18cornea4","E18cornea5","E8retina1","E18brain1","E18retina3","E8retina3","E18retina2","E18brain2","E18brain3","E18retina1"),stringsAsFactors=F)
# create the df
df = data.frame(lab=sampinfo[pmatch(sampinfo[,1],mapped[,1]),2],raw=totread,map=mapped[,3])
# output the stats as csv
write.table(x=df,file=file.path(outdir,"readstat.csv"),sep=",",quote=F,row.names=F)

df.melt = melt(df)
df.melt = na.omit(df.melt)

g = ggplot(df.melt,aes(x=lab,y=value))+geom_bar(stat="identity",aes(fill=variable),position="dodge")+theme_bw()

# plot bar graph
pdf(file.path(plotdir,"readstat.pdf"))
print(g)
dev.off()
