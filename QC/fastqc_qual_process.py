#!/usr/bin/python

import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def getlines(fname):
    with open(fname) as fh:
        lines = fh.readlines()
    return lines

def getquals(lines,which="mean"):
    titles=lines[12]
    qlines = lines[13:80]
    fields=[s.strip().split("\t") for s in qlines]
    idx=titles.strip().split("\t").index(which)
    data=[line[idx] for line in fields]
    return data

def getcount(lines,which="Count"):
    titles = lines[152]
    endline = [i for i in range(len(lines)) if lines[i][0]==">" and i>152][0]
    clines = lines[153:endline]
    fields = [s.strip().split("\t") for s in clines]
    idx = titles.strip().split("\t").index(which)
    data=[line[idx] for line in fields]
    return data

def listMean(numlist):
    numlist=np.array(numlist,dtype="float")
    return np.mean(numlist,axis=0)

def plotData(x,ylist):
    xnum = [float(c.split("-")[0]) for c in x]
    fig=plt.figure()
    plt.plot(xnum,ylist[0],label="Read 1")
    plt.plot(xnum,ylist[1],label="Read 2")
    plt.ylim([0,38])
    plt.xlim([0,max(xnum)])
    plt.legend(loc=3)
    fig.savefig('means.pdf',format='pdf')

def countMean(qscore,qcount):
    scorelist=[int(item) for sublist in qscore for item in sublist]
    countlist = [float(item) for sublist in qcount for item in sublist]
    minq = min(scorelist)
    maxq = max(scorelist)
    qlist = list(range(minq,maxq+1))
    clist = []
    for i,q in enumerate(qlist):
        clist.append([])
        [clist[i].append(countlist[j]) for j in range(len(countlist)) if scorelist[j]==q]
    cmean = [np.mean(l) for l in clist]
    return qlist,cmean


def main():
    # get files from cwd
    cwd = "/home/isac/Dropbox/Data/Genetics/MethSeq/150415_chicken/QC/zips/extracted/data"
    outd = "/home/isac/Dropbox/Data/Genetics/MethSeq/170120_chicken/QC"
    os.chdir(cwd)
    files = [f for f in os.listdir(cwd) if f.split(".")[-1]=="txt"] 
    read1 = [f for f in files if f.split("_")[2]=="1"]
    read2 = [f for f in files if f.split("_")[2]=="2"]
    meanlist=[] 
    countlist=[]
    for reads in [read1,read2]:
        lineslist= [getlines(s) for s in reads]
        base=getquals(lineslist[0],which="#Base")
        means=[getquals(l,"Mean") for l in lineslist]
        meanlist.append(listMean(means))
        qscore=[getcount(l,"#Quality") for l in lineslist]
        qcount=[getcount(l,"Count") for l in lineslist]
        scores,count=countMean(qscore,qcount)
        countlist.append(count)
    #plotData(base,meanlist)
    with open(os.path.join(outd,"meanqualvcycle.txt"),"w") as fh:
        fh.write("base\tread1\tread2\n")
        [fh.write("{}\t{}\t{}\n".format(base[i],meanlist[0][i],meanlist[1][i])) for i in range(len(base))]
    with open(os.path.join(outd,"qualitycount.txt"),"w") as fh:
        fh.write("score\tread1\tread2\n")
        [fh.write("{}\t{}\t{}\n".format(scores[i],countlist[0][i],countlist[1][i])) for i in range(len(scores))]

if __name__ == "__main__":
    main()
