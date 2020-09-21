#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
assembly_name=args[1]
mappability=args[2]
repeatmasked=args[3]
gcconentcalc=args[4]
total_length=args[5]

print(paste("assembly_name: ", assembly_name, "mappability: ",mappability, " repeatmasked: ", repeatmasked, "gcconentcalc:", gcconentcalc, " total_length: ", total_length))

chrY_M<-as.data.frame(read.table(mappability, sep="\t", header=TRUE, stringsAsFactors=FALSE))
chrY_RM<-as.data.frame(read.table(repeatmasked, sep="\t", header=TRUE, stringsAsFactors=FALSE))
chrY_GC<-as.data.frame(read.table(gcconentcalc, sep="\t", header=TRUE, stringsAsFactors=FALSE))
chrY_Map<-as.data.frame(rep(1,total_length)) #mock array
chrY_TRF<-as.data.frame(rep(1,total_length)) #mock array
colnames(chrY_Map)="MapCount"
colnames(chrY_TRF)="TRF"

print(paste("mappability:",nrow(chrY_M)))
print(paste("repeatmasked:",nrow(chrY_RM)))
print(paste("gcconentcalc:",nrow(chrY_GC)))

stopifnot(identical(nrow(chrY_M),nrow(chrY_RM),nrow(chrY_GC),nrow(chrY_Map))==TRUE)

Position<-seq(1,total_length)
SummaryY<-cbind(Position,chrY_M,chrY_RM,chrY_TRF,chrY_GC,chrY_Map)

print(head(SummaryY))
print(tail(SummaryY))

id1=which(SummaryY$RepeatMasked==1)
id2=which(SummaryY$TRF==1)
id=intersect(id1,id2)

x=SummaryY[id,]
x_M=x[which(x$Mappability>0),]
x_MG=x_M[which(x_M$GC>0),]

output_file=paste0(assembly_name,"/","Summary_",assembly_name,"_chrY_Map_RM_TRF_GC_MapC.tab")
write.table(SummaryY,file=output_file,sep = "\t",quote = FALSE,row.names = FALSE,col.names = TRUE)