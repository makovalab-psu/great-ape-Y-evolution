#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
mapping_Q_file=args[1]
copy_N_file=args[2]
folder=args[3]

library(tidyr)
library(dplyr)
library(zoo)
#par(mfrow=c(3,3))


format_window <- function(window) {
  ncolumns<-length(rev(unlist(strsplit(as.character(window),"_"))))
  wnames<-rev(unlist(strsplit(as.character(window),"_")))[3:ncolumns]
  formatted_wname<-paste(rev(wnames),collapse="_")
  wstart<-as.numeric(rev(unlist(strsplit(as.character(window),"_")))[2]) + 1 #CHANGE THE NUMBERING SYSTEM
  wend<-as.numeric(rev(unlist(strsplit(as.character(window),"_")))[1])
  return(c(formatted_wname,wstart,wend))
}

mapping_Q<-read.table(mapping_Q_file, sep="\t",header=TRUE) #e.g. "/Users/polly/Desktop/projects/classification/orangutan/sorang.1991-0051.assembly.fa.MAPQ.txt"
mapping_Q<-as.data.frame(mapping_Q)
mapping_Q<-mapping_Q[mapping_Q$flag==0 | mapping_Q$flag==16,] #only keep primary mappings, remove unmapped reads (typically super short windows)

#mapping_Q<-separate(data=mapping_Q,col=window,into=c("name","start","end"),sep="_")

mapping_Q$name<-sapply(mapping_Q$window, function(x) {format_window(x)[1]})
mapping_Q$start<-sapply(mapping_Q$window, function(x) {format_window(x)[2]})
mapping_Q$end<-sapply(mapping_Q$window, function(x) {format_window(x)[3]})
  
mapping_Q$window<-paste0(mapping_Q$name,"_",as.numeric(mapping_Q$start),"_",as.numeric(mapping_Q$end))

copy_N<-read.table(copy_N_file, sep="\t",header=TRUE) #e.g. /Users/polly/Desktop/projects/classification/orangutan/sorang_subset_1_on_sorang_100mil_SE.bamwindow_CN_summary.txt
colnames(copy_N)<-c("window","CN_MAP","CN_XDEG")
copy_N$window<-gsub('\\:','_',gsub('\\-','_',copy_N$window))

mdf<-merge(mapping_Q,copy_N,all.x=TRUE)
mdf<-mdf[order(mdf$contig,as.numeric(as.character(mdf$start))),]

#boxplot(mdf$CN_MAP,outline=FALSE,col="gold",main="Boxplot of coverage across the whole dataset")
median_coverage<-median(mdf$CN_MAP,na.rm=TRUE)
#boxplot(mdf[mdf$MAPQ==0,]$CN_MAP,outline=FALSE,col="gold",main="Copy number of windows with MAPQ 0")

mdf$class<-rep(NA,nrow(mdf))
print(dim(mdf))
nonrepetitive<-mdf[mdf$CN_MAP>0,] #CN_MAP must be greater than 0, otherwise the window is at risk of being repetitive
print(dim(nonrepetitive))
#plot(table(nonrepetitive$MAPQ),main="MAPQ of presumabely non-repetitive regions")

#in this subset we can finally flag those with MAPQ=0 as ampliconic
copy_number_threshold<-0.3
intersimilar_windows<-mdf[mdf$CN_MAP>0 & mdf$MAPQ==0,]$window
collapsed_windows<-mdf[mdf$MAPQ==60 & (mdf$CN_MAP>(2*copy_number_threshold)),]$window

unique_windows<-mdf[mdf$MAPQ==60 & (mdf$CN_MAP<copy_number_threshold),]$window

tryCatch(
        {
        mdf[mdf$window %in% intersimilar_windows,]$class<-2 #set as ampliconic
        },
        error=function(error_message) {
            message("NO intersimilar windows exist.")
            message("And below is the error message from R:")
            message(error_message)
        }
    )

tryCatch(
        {
        mdf[mdf$window %in% collapsed_windows,]$class<-2 #set as ampliconic
        },
        error=function(error_message) {
            message("NO collapsed windows exist.")
            message("And below is the error message from R:")
            message(error_message)
        }
    )

mdf[mdf$window %in% na.omit(unique_windows),]$class<-1 #set as other
table(mdf$class)

#create a plot that will be used to pick the threshold
h<-hist(mdf[mdf$CN_MAP<10 & mdf$CN_MAP>0,]$CN_MAP,breaks=200)
plot(h, main="Histogram of the copy number")