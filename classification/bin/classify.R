#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
mapping_Q_file=args[1]
copy_N_file=args[2]
folder=args[3]
copy_number_threshold=as.numeric(as.character(args[4]))

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
intersimilar_windows<-mdf[mdf$CN_MAP>0 & mdf$MAPQ==0,]$window
collapsed_windows<-mdf[mdf$MAPQ==60 & (mdf$CN_MAP>=(copy_number_threshold)),]$window

unique_windows<-mdf[mdf$MAPQ==60 & (mdf$CN_MAP<copy_number_threshold),]$window

#keep going if either group does not exist
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

palette(c("gold","cornflowerblue"))
h<-hist(mdf[mdf$CN_MAP<10 & mdf$CN_MAP>0,]$CN_MAP,breaks=200)
cuts <- cut(h$breaks, c(-Inf,copy_number_threshold,Inf))
plot(h, col=cuts, main=folder,xlab="Copy number")
abline(v=copy_number_threshold,lw=1,col="red")

##interpolate missing values
interpolated<-mdf %>%
  group_by(contig) %>%
  mutate(zoo_interpolate = na.approx(class, na.rm=FALSE))  

#filled 
filled<-interpolated
filled$filledVal<-filled$zoo_interpolate

#final decision regarding the contig (median and rounding)
final<-filled %>%
  group_by(contig) %>%
  mutate(finalGroup = round(median(filledVal,na.rm=TRUE)))

contig_statistics<-aggregate(finalGroup ~ contig, data=final, function(x) {round(median(x,na.rm=TRUE))})

print("final classification")
table(contig_statistics$finalGroup)
length(unique(sort(final$contig)))
#View(final)

ampliconic_names<-contig_statistics[contig_statistics$finalGroup==2,]$contig
XDEG_names<-contig_statistics[contig_statistics$finalGroup==1,]$contig

XDEG_subset<-mdf[mdf$contig %in% XDEG_names,]
ampl_subset<-mdf[mdf$contig %in% ampliconic_names,]

XDEG_subset<-XDEG_subset[XDEG_subset$CN_MAP>0,]
ampl_subset<-ampl_subset[ampl_subset$CN_MAP>0,]

boxplot(XDEG_subset$CN_MAP,col="gold",outline=FALSE,ylim=c(0,5),main="Copy number of X-degenerate scaffolds")
boxplot(ampl_subset$CN_MAP,col="cyan",outline=FALSE,ylim=c(0,5),main="Copy number of ampliconic scaffolds")

write.table(sort(XDEG_names), file = paste0(folder,"/other.txt"), sep = "\t", quote=FALSE, col.names = FALSE, row.names = FALSE)
write.table(sort(ampliconic_names), file = paste0(folder,"/ampliconic.txt"), sep = "\t", quote=FALSE, col.names = FALSE, row.names = FALSE)

gff<-paste(final$contig,"classification","class",(as.numeric(as.character(final$start))+1),final$end,final$filledVal,".",".",final$contig,sep='\t')
write.table(gff, file = paste0(folder,"/classification.gff"), sep = "\t", quote=FALSE, col.names = FALSE, row.names = FALSE)