cd /nfs/brubeck.bx.psu.edu/scratch6/rahul/Palindrome/analysis/plot_CN
source activate /galaxy/home/rxv923/anaconda2/envs/renv/
R

#################################PALINDROME COPY NUMBER HUMAN

palindromeCN_Gdir="/nfs/brubeck.bx.psu.edu/scratch6/rahul/Gorilla_Y/analysis/palindrome/human/gcCor_results/"
palindromeCN_Bdir="/nfs/brubeck.bx.psu.edu/scratch6/rahul/Bonobo_Y/analysis/palindrome/human/gcCor_results/"
palindromeCN_Odir="/nfs/brubeck.bx.psu.edu/scratch6/rahul/Orangutan_Y/analysis/palindrome/human/gcCor_results/"


get_coverage<-function(data_align){
if (length(data_align)> 1){
	temp=data_align[which(data_align$V4<200 & data_align$V13 > 0.8),]
	return (temp$V10/temp$V12)
} else{ return=0}
}

get_CN<-function(coverage,XDG){
control=median(XDG)
if (length(coverage) > 1){
	temp=coverage/control
}else{temp=0}
return (temp)
}



######################################GORILLA
species='gorilla'
path=palindromeCN_Gdir


#If file is empty does not print out error.
#Trying to overcome  missing data for few palindromes
palin=palin1=palin2=palin3=palin4=palin5=palin6=palin7=palin8=XDG=0

palin=read.table(paste(path,"Coverage_",species,"_msY_WGS_palindromeAll_80per.bed",sep=""), sep="\t", stringsAsFactors=F, header=F)
palin1=read.table(paste(path,"Coverage_",species,"_msY_WGS_palindrome1_80per.bed",sep=""), sep="\t", stringsAsFactors=F, header=F)
palin2=read.table(paste(path,"Coverage_",species,"_msY_WGS_palindrome2_80per.bed",sep=""), sep="\t", stringsAsFactors=F, header=F)
palin3=read.table(paste(path,"Coverage_",species,"_msY_WGS_palindrome3_80per.bed",sep=""), sep="\t", stringsAsFactors=F, header=F)
palin4=read.table(paste(path,"Coverage_",species,"_msY_WGS_palindrome4_80per.bed",sep=""), sep="\t", stringsAsFactors=F, header=F)
palin5=read.table(paste(path,"Coverage_",species,"_msY_WGS_palindrome5_80per.bed",sep=""), sep="\t", stringsAsFactors=F, header=F)
palin6=read.table(paste(path,"Coverage_",species,"_msY_WGS_palindrome6_80per.bed",sep=""), sep="\t", stringsAsFactors=F, header=F)
palin7=read.table(paste(path,"Coverage_",species,"_msY_WGS_palindrome7_80per.bed",sep=""), sep="\t", stringsAsFactors=F, header=F)
palin8=read.table(paste(path,"Coverage_",species,"_msY_WGS_palindrome8_80per.bed",sep=""), sep="\t", stringsAsFactors=F, header=F)
XDG=read.table(paste(path,"Coverage_",species,"_msY_WGS_XDG_80per.bed",sep=""), sep="\t", stringsAsFactors=F, header=F)



#Parse the coverage from bedtools output
pall=get_coverage(palin)
p1=get_coverage(palin1)
p2=get_coverage(palin2)
p3=get_coverage(palin3)
p4=get_coverage(palin4)
p5=get_coverage(palin5)
p6=get_coverage(palin6)
p7=get_coverage(palin7)
p8=get_coverage(palin8)
x=get_coverage(XDG)

#Calculate the copy number
CNpall=get_CN(pall,x)
CNp1=get_CN(p1,x)
CNp2=get_CN(p2,x)
CNp3=get_CN(p3,x)
CNp4=get_CN(p4,x)
CNp5=get_CN(p5,x)
CNp6=get_CN(p6,x)
CNp7=get_CN(p7,x)
CNp8=get_CN(p8,x)
CNx=get_CN(x,x)


#Generate ggplot friendly matrix
out<-rbind(
cbind(CNx,rep("XDG",length(CNx)),species),
cbind(CNp1,rep("P1",length(CNp1)),species),
cbind(CNp2,rep("P2",length(CNp2)),species),
cbind(CNp3,rep("P3",length(CNp3)),species),
cbind(CNp4,rep("P4",length(CNp4)),species),
cbind(CNp5,rep("P5",length(CNp5)),species),
cbind(CNp6,rep("P6",length(CNp6)),species),
cbind(CNp7,rep("P7",length(CNp7)),species),
cbind(CNp8,rep("P8",length(CNp8)),species)

)

#Median number for Supplemtary table
round(median(CNp1),digits=2)
round(median(CNp2),digits=2)
round(median(CNp3),digits=2)
round(median(CNp4),digits=2)
round(median(CNp5),digits=2)
round(median(CNp6),digits=2)
round(median(CNp7),digits=2)
round(median(CNp8),digits=2)

gor_cov<-out


######################################BONOBO
species='bonobo'
path=palindromeCN_Bdir

#If file is empty does not print out error.
#Trying to overcome  missing data for few palindromes
palin=palin1=palin2=palin3=palin4=palin5=palin6=palin7=palin8=XDG=0

palin=read.table(paste(path,"Coverage_",species,"_msY_WGS_palindromeAll_80per.bed",sep=""), sep="\t", stringsAsFactors=F, header=F)
palin1=read.table(paste(path,"Coverage_",species,"_msY_WGS_palindrome1_80per.bed",sep=""), sep="\t", stringsAsFactors=F, header=F)
palin2=read.table(paste(path,"Coverage_",species,"_msY_WGS_palindrome2_80per.bed",sep=""), sep="\t", stringsAsFactors=F, header=F)
palin3=read.table(paste(path,"Coverage_",species,"_msY_WGS_palindrome3_80per.bed",sep=""), sep="\t", stringsAsFactors=F, header=F)
palin4=read.table(paste(path,"Coverage_",species,"_msY_WGS_palindrome4_80per.bed",sep=""), sep="\t", stringsAsFactors=F, header=F)
palin5=read.table(paste(path,"Coverage_",species,"_msY_WGS_palindrome5_80per.bed",sep=""), sep="\t", stringsAsFactors=F, header=F)
palin6=read.table(paste(path,"Coverage_",species,"_msY_WGS_palindrome6_80per.bed",sep=""), sep="\t", stringsAsFactors=F, header=F)
palin7=read.table(paste(path,"Coverage_",species,"_msY_WGS_palindrome7_80per.bed",sep=""), sep="\t", stringsAsFactors=F, header=F)
palin8=read.table(paste(path,"Coverage_",species,"_msY_WGS_palindrome8_80per.bed",sep=""), sep="\t", stringsAsFactors=F, header=F)
XDG=read.table(paste(path,"Coverage_",species,"_msY_WGS_XDG_80per.bed",sep=""), sep="\t", stringsAsFactors=F, header=F)



#Parse the coverage from bedtools output
pall=get_coverage(palin)
p1=get_coverage(palin1)
p2=get_coverage(palin2)
p3=get_coverage(palin3)
p4=get_coverage(palin4)
p5=get_coverage(palin5)
p6=get_coverage(palin6)
p7=get_coverage(palin7)
p8=get_coverage(palin8)
x=get_coverage(XDG)

#Calculate the copy number
CNpall=get_CN(pall,x)
CNp1=get_CN(p1,x)
CNp2=get_CN(p2,x)
CNp3=get_CN(p3,x)
CNp4=get_CN(p4,x)
CNp5=get_CN(p5,x)
CNp6=get_CN(p6,x)
CNp7=get_CN(p7,x)
CNp8=get_CN(p8,x)
CNx=get_CN(x,x)


#Generate ggplot friendly matrix
out<-rbind(
cbind(CNx,rep("XDG",length(CNx)),species),
cbind(CNp1,rep("P1",length(CNp1)),species),
cbind(CNp2,rep("P2",length(CNp2)),species),
cbind(CNp3,rep("P3",length(CNp3)),species),
cbind(CNp4,rep("P4",length(CNp4)),species),
cbind(CNp5,rep("P5",length(CNp5)),species),
cbind(CNp6,rep("P6",length(CNp6)),species),
cbind(CNp7,rep("P7",length(CNp7)),species),
cbind(CNp8,rep("P8",length(CNp8)),species)

)

#Median number for Supplemtary table
round(median(CNp1),digits=2)
round(median(CNp2),digits=2)
round(median(CNp3),digits=2)
round(median(CNp4),digits=2)
round(median(CNp5),digits=2)
round(median(CNp6),digits=2)
round(median(CNp7),digits=2)
round(median(CNp8),digits=2)

bon_cov<-out

######################################ORANGUTAN

species='orangutan'
path=palindromeCN_Odir

#If file is empty does not print out error.
#Trying to overcome  missing data for few palindromes
palin=palin1=palin2=palin3=palin4=palin5=palin6=palin7=palin8=XDG=0

palin=read.table(paste(path,"Coverage_",species,"_msY_WGS_palindromeAll_80per.bed",sep=""), sep="\t", stringsAsFactors=F, header=F)
palin1=read.table(paste(path,"Coverage_",species,"_msY_WGS_palindrome1_80per.bed",sep=""), sep="\t", stringsAsFactors=F, header=F)
palin2=read.table(paste(path,"Coverage_",species,"_msY_WGS_palindrome2_80per.bed",sep=""), sep="\t", stringsAsFactors=F, header=F)
palin3=read.table(paste(path,"Coverage_",species,"_msY_WGS_palindrome3_80per.bed",sep=""), sep="\t", stringsAsFactors=F, header=F)
palin4=read.table(paste(path,"Coverage_",species,"_msY_WGS_palindrome4_80per.bed",sep=""), sep="\t", stringsAsFactors=F, header=F)
palin5=read.table(paste(path,"Coverage_",species,"_msY_WGS_palindrome5_80per.bed",sep=""), sep="\t", stringsAsFactors=F, header=F)
palin6=read.table(paste(path,"Coverage_",species,"_msY_WGS_palindrome6_80per.bed",sep=""), sep="\t", stringsAsFactors=F, header=F)
palin7=read.table(paste(path,"Coverage_",species,"_msY_WGS_palindrome7_80per.bed",sep=""), sep="\t", stringsAsFactors=F, header=F)
palin8=read.table(paste(path,"Coverage_",species,"_msY_WGS_palindrome8_80per.bed",sep=""), sep="\t", stringsAsFactors=F, header=F)
XDG=read.table(paste(path,"Coverage_",species,"_msY_WGS_XDG_80per.bed",sep=""), sep="\t", stringsAsFactors=F, header=F)



#Parse the coverage from bedtools output
pall=get_coverage(palin)
p1=get_coverage(palin1)
p2=get_coverage(palin2)
p3=get_coverage(palin3)
p4=get_coverage(palin4)
p5=get_coverage(palin5)
p6=get_coverage(palin6)
p7=get_coverage(palin7)
p8=get_coverage(palin8)
x=get_coverage(XDG)

#Calculate the copy number
CNpall=get_CN(pall,x)
CNp1=get_CN(p1,x)
CNp2=get_CN(p2,x)
CNp3=get_CN(p3,x)
CNp4=get_CN(p4,x)
CNp5=get_CN(p5,x)
CNp6=get_CN(p6,x)
CNp7=get_CN(p7,x)
CNp8=get_CN(p8,x)
CNx=get_CN(x,x)


#Generate ggplot friendly matrix
out<-rbind(
cbind(CNx,rep("XDG",length(CNx)),species),
cbind(CNp1,rep("P1",length(CNp1)),species),
cbind(CNp2,rep("P2",length(CNp2)),species),
cbind(CNp3,rep("P3",length(CNp3)),species),
cbind(CNp4,rep("P4",length(CNp4)),species),
cbind(CNp5,rep("P5",length(CNp5)),species),
cbind(CNp6,rep("P6",length(CNp6)),species),
cbind(CNp7,rep("P7",length(CNp7)),species),
cbind(CNp8,rep("P8",length(CNp8)),species)

)

#Median number for Supplemtary table
round(median(CNp1),digits=2)
round(median(CNp2),digits=2)
round(median(CNp3),digits=2)
round(median(CNp4),digits=2)
round(median(CNp5),digits=2)
round(median(CNp6),digits=2)
round(median(CNp7),digits=2)
round(median(CNp8),digits=2)

ora_cov<-out

##################################

#Colors.name <-  c(bo="#A7C855", ch="#74A8D3" , hu="6A6A6A", go="#BD3231", or="#E19F5A")
Colors.name <-  c(bo="#ABD9E9", ch="#2C7BB6" , hu="D7191C", go="#FFFFBF", or="#FDAE61")
darkColors.name <-  c(bo="#799630", ch="#3778ad" , hu="#4a4a4a", go="#842222", or="#ba6f21")


CN_human_palindrome<-rbind(gor_cov,bon_cov,ora_cov)
colnames(CN_human_palindrome)<-c("CopyNumber","Palindrome","Species")
CN_human_palindrome<-as.data.frame(CN_human_palindrome,stringsAsFactors = FALSE)
CN_human_palindrome$CopyNumber<-as.numeric(CN_human_palindrome$CopyNumber)
CN_human_palindrome$Species<-substr(CN_human_palindrome$Species,0,2)
require(ggplot2)


CN_log_human_palindrome<-CN_human_palindrome
CN_log_human_palindrome$CopyNumber<-log(as.numeric(CN_log_human_palindrome$CopyNumber))


human1_obj<-ggplot(CN_human_palindrome, aes(x=Species, y=CopyNumber))+
	geom_boxplot(aes(fill=Species),size=0,outlier.shape = NA, width=0.5)+scale_colour_manual(values=Colors.name)+ scale_fill_manual(values=Colors.name)+
	#geom_point(size=0.5,aes(color=Species),shape = 18)+scale_colour_manual(values=darkColors.name)+
	geom_hline(yintercept = 1, linetype = 2,size=0.1)+
	geom_hline(yintercept = 2, linetype = 2,color='red',size=0.1)+
	ylim(c(0,25))+
	theme_classic() + 
	#theme_bw() + 
	#theme(axis.title.x = element_blank(),axis.text.x = element_text(),plot.title = element_blank())+
	#labs(title = "CN")+
	#labs(y = "Copy number", color = "Species\n")+
	facet_wrap( ~ Palindrome, nrow = 1)
	
pdf("Fig_Human_copynumber.pdf")
human1_obj
dev.off()

humanLog_obj<-ggplot(CN_log_human_palindrome, aes(x=Species, y=CopyNumber))+
	geom_boxplot(aes(fill=Species),size=0.1,outlier.shape = NA, width=0.8)+scale_colour_manual(values=Colors.name)+ scale_fill_manual(values=Colors.name)+
	geom_hline(yintercept = 0, linetype = 2,size=0.1)+
	geom_hline(yintercept = 0.6931472, linetype = 2,color='red',size=0.1)+
	scale_y_continuous(limits=c(-0.75,3.218876),breaks=c(-0.6931472,0,0.6931472,1.609438,2.302585,2.70805,2.995732),labels=c("-0.6931472"="0.5","0" = "1","0.6931472" = "2", "1.609438"="5","2.302585"="10", "2.70805" = "15", "2.995732" = "20"))+
	theme_classic(base_size = 10) + 
	theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())+
	facet_wrap( ~ Palindrome, nrow = 1)

pdf("Fig_Human_copynumberLOG.pdf")
humanLog_obj
dev.off()

