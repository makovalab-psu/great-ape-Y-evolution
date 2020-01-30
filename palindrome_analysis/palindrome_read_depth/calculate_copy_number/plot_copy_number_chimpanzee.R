cd /nfs/brubeck.bx.psu.edu/scratch6/rahul/Palindrome/analysis/plot_CN
source activate /galaxy/home/rxv923/anaconda2/envs/renv/
R

#################################PALINDROME COPY NUMBER CHIMPANZEE

palindromeCN_Gdir="/nfs/brubeck.bx.psu.edu/scratch6/rahul/Gorilla_Y/analysis/palindrome/chimpanzee/gcCor_results/"
palindromeCN_Bdir="/nfs/brubeck.bx.psu.edu/scratch6/rahul/Bonobo_Y/analysis/palindrome/chimpanzee/gcCor_results/"
palindromeCN_Odir="/nfs/brubeck.bx.psu.edu/scratch6/rahul/Orangutan_Y/analysis/palindrome/chimpanzee/gcCor_results/"


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
palin=palin1=palin2=palin3=palin4=palin5=palin6=palin7=palin8=palin9=palin10=palin11=palin12=palin13=palin14=palin15=palin16=palin17=palin18=palin19=XDG=0

palin=read.table(paste(path,"Coverage_",species,"_msY_WGS_palindromeAll_80per.bed",sep=""), sep="\t", stringsAsFactors=F, header=F)
palin1=read.table(paste(path,"Coverage_",species,"_msY_WGS_palindrome1_80per.bed",sep=""), sep="\t", stringsAsFactors=F, header=F)
palin2=read.table(paste(path,"Coverage_",species,"_msY_WGS_palindrome2_80per.bed",sep=""), sep="\t", stringsAsFactors=F, header=F)
palin3=read.table(paste(path,"Coverage_",species,"_msY_WGS_palindrome3_80per.bed",sep=""), sep="\t", stringsAsFactors=F, header=F)
palin4=read.table(paste(path,"Coverage_",species,"_msY_WGS_palindrome4_80per.bed",sep=""), sep="\t", stringsAsFactors=F, header=F)
palin5=read.table(paste(path,"Coverage_",species,"_msY_WGS_palindrome5_80per.bed",sep=""), sep="\t", stringsAsFactors=F, header=F)
palin6=read.table(paste(path,"Coverage_",species,"_msY_WGS_palindrome6_80per.bed",sep=""), sep="\t", stringsAsFactors=F, header=F)
palin7=read.table(paste(path,"Coverage_",species,"_msY_WGS_palindrome7_80per.bed",sep=""), sep="\t", stringsAsFactors=F, header=F)
palin8=read.table(paste(path,"Coverage_",species,"_msY_WGS_palindrome8_80per.bed",sep=""), sep="\t", stringsAsFactors=F, header=F)
palin9=read.table(paste(path,"Coverage_",species,"_msY_WGS_palindrome9_80per.bed",sep=""), sep="\t", stringsAsFactors=F, header=F)
palin10=read.table(paste(path,"Coverage_",species,"_msY_WGS_palindrome10_80per.bed",sep=""), sep="\t", stringsAsFactors=F, header=F)
palin11=read.table(paste(path,"Coverage_",species,"_msY_WGS_palindrome11_80per.bed",sep=""), sep="\t", stringsAsFactors=F, header=F)
palin12=read.table(paste(path,"Coverage_",species,"_msY_WGS_palindrome12_80per.bed",sep=""), sep="\t", stringsAsFactors=F, header=F)
palin13=read.table(paste(path,"Coverage_",species,"_msY_WGS_palindrome13_80per.bed",sep=""), sep="\t", stringsAsFactors=F, header=F)
palin14=read.table(paste(path,"Coverage_",species,"_msY_WGS_palindrome14_80per.bed",sep=""), sep="\t", stringsAsFactors=F, header=F)
palin15=read.table(paste(path,"Coverage_",species,"_msY_WGS_palindrome15_80per.bed",sep=""), sep="\t", stringsAsFactors=F, header=F)
palin16=read.table(paste(path,"Coverage_",species,"_msY_WGS_palindrome16_80per.bed",sep=""), sep="\t", stringsAsFactors=F, header=F)
palin17=read.table(paste(path,"Coverage_",species,"_msY_WGS_palindrome17_80per.bed",sep=""), sep="\t", stringsAsFactors=F, header=F)
palin18=read.table(paste(path,"Coverage_",species,"_msY_WGS_palindrome18_80per.bed",sep=""), sep="\t", stringsAsFactors=F, header=F)
palin19=read.table(paste(path,"Coverage_",species,"_msY_WGS_palindrome19_80per.bed",sep=""), sep="\t", stringsAsFactors=F, header=F)
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
p9=get_coverage(palin9)
p10=get_coverage(palin10)
p11=get_coverage(palin11)
p12=get_coverage(palin12)
p13=get_coverage(palin13)
p14=get_coverage(palin14)
p15=get_coverage(palin15)
p16=get_coverage(palin16)
p17=get_coverage(palin17)
p18=get_coverage(palin18)
p19=get_coverage(palin19)

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
CNp9=get_CN(p9,x)
CNp10=get_CN(p10,x)
CNp11=get_CN(p11,x)
CNp12=get_CN(p12,x)
CNp13=get_CN(p13,x)
CNp14=get_CN(p14,x)
CNp15=get_CN(p15,x)
CNp16=get_CN(p16,x)
CNp17=get_CN(p17,x)
CNp18=get_CN(p18,x)
CNp19=get_CN(p19,x)

CNx=get_CN(x,x)

#Cluster similar palindromes
C1<-c(CNp1,CNp6,CNp8,CNp10,CNp14,CNp16)
C2<-c(CNp2,CNp11,CNp15)
C3<-c(CNp3,CNp12)
C4<-c(CNp4,CNp13)
C5<-c(CNp5,CNp7,CNp9)

#Generate ggplot friendly matrix
out<-rbind(
cbind(CNx,rep("XDG",length(CNx)),species),
cbind(C1,rep("C1",length(C1)),species),
cbind(C2,rep("C2",length(C2)),species),
cbind(C3,rep("C3",length(C3)),species),
cbind(C4,rep("C4",length(C4)),species),
cbind(C5,rep("C5",length(C5)),species),
cbind(CNp17,rep("C17",length(CNp17)),species),
cbind(CNp18,rep("C18",length(CNp18)),species),
cbind(CNp19,rep("C19",length(CNp19)),species)
)

#Median number for Supplemtary table
round(median(C1[which(C1>0)]),digits=2)
round(median(C2[which(C2>0)]),digits=2)
round(median(C3[which(C3>0)]),digits=2)
round(median(C4[which(C4>0)]),digits=2)
round(median(C5[which(C5>0)]),digits=2)
round(median(CNp17[which(CNp17>0)]),digits=2)
round(median(CNp18[which(CNp18>0)]),digits=2)
round(median(CNp19[which(CNp19>0)]),digits=2)

gor_cov<-out


######################################BONOBO
species='bonobo'
path=palindromeCN_Bdir

#If file is empty does not print out error.
#Trying to overcome  missing data for few palindromes
palin=palin1=palin2=palin3=palin4=palin5=palin6=palin7=palin8=palin9=palin10=palin11=palin12=palin13=palin14=palin15=palin16=palin17=palin18=palin19=XDG=0

palin=read.table(paste(path,"Coverage_",species,"_msY_WGS_palindromeAll_80per.bed",sep=""), sep="\t", stringsAsFactors=F, header=F)
palin1=read.table(paste(path,"Coverage_",species,"_msY_WGS_palindrome1_80per.bed",sep=""), sep="\t", stringsAsFactors=F, header=F)
palin2=read.table(paste(path,"Coverage_",species,"_msY_WGS_palindrome2_80per.bed",sep=""), sep="\t", stringsAsFactors=F, header=F)
palin3=read.table(paste(path,"Coverage_",species,"_msY_WGS_palindrome3_80per.bed",sep=""), sep="\t", stringsAsFactors=F, header=F)
palin4=read.table(paste(path,"Coverage_",species,"_msY_WGS_palindrome4_80per.bed",sep=""), sep="\t", stringsAsFactors=F, header=F)
palin5=read.table(paste(path,"Coverage_",species,"_msY_WGS_palindrome5_80per.bed",sep=""), sep="\t", stringsAsFactors=F, header=F)
palin6=read.table(paste(path,"Coverage_",species,"_msY_WGS_palindrome6_80per.bed",sep=""), sep="\t", stringsAsFactors=F, header=F)
palin7=read.table(paste(path,"Coverage_",species,"_msY_WGS_palindrome7_80per.bed",sep=""), sep="\t", stringsAsFactors=F, header=F)
palin8=read.table(paste(path,"Coverage_",species,"_msY_WGS_palindrome8_80per.bed",sep=""), sep="\t", stringsAsFactors=F, header=F)
palin9=read.table(paste(path,"Coverage_",species,"_msY_WGS_palindrome9_80per.bed",sep=""), sep="\t", stringsAsFactors=F, header=F)
palin10=read.table(paste(path,"Coverage_",species,"_msY_WGS_palindrome10_80per.bed",sep=""), sep="\t", stringsAsFactors=F, header=F)
palin11=read.table(paste(path,"Coverage_",species,"_msY_WGS_palindrome11_80per.bed",sep=""), sep="\t", stringsAsFactors=F, header=F)
palin12=read.table(paste(path,"Coverage_",species,"_msY_WGS_palindrome12_80per.bed",sep=""), sep="\t", stringsAsFactors=F, header=F)
palin13=read.table(paste(path,"Coverage_",species,"_msY_WGS_palindrome13_80per.bed",sep=""), sep="\t", stringsAsFactors=F, header=F)
palin14=read.table(paste(path,"Coverage_",species,"_msY_WGS_palindrome14_80per.bed",sep=""), sep="\t", stringsAsFactors=F, header=F)
palin15=read.table(paste(path,"Coverage_",species,"_msY_WGS_palindrome15_80per.bed",sep=""), sep="\t", stringsAsFactors=F, header=F)
palin16=read.table(paste(path,"Coverage_",species,"_msY_WGS_palindrome16_80per.bed",sep=""), sep="\t", stringsAsFactors=F, header=F)
palin17=read.table(paste(path,"Coverage_",species,"_msY_WGS_palindrome17_80per.bed",sep=""), sep="\t", stringsAsFactors=F, header=F)
palin18=read.table(paste(path,"Coverage_",species,"_msY_WGS_palindrome18_80per.bed",sep=""), sep="\t", stringsAsFactors=F, header=F)
palin19=read.table(paste(path,"Coverage_",species,"_msY_WGS_palindrome19_80per.bed",sep=""), sep="\t", stringsAsFactors=F, header=F)
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
p9=get_coverage(palin9)
p10=get_coverage(palin10)
p11=get_coverage(palin11)
p12=get_coverage(palin12)
p13=get_coverage(palin13)
p14=get_coverage(palin14)
p15=get_coverage(palin15)
p16=get_coverage(palin16)
p17=get_coverage(palin17)
p18=get_coverage(palin18)
p19=get_coverage(palin19)

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
CNp9=get_CN(p9,x)
CNp10=get_CN(p10,x)
CNp11=get_CN(p11,x)
CNp12=get_CN(p12,x)
CNp13=get_CN(p13,x)
CNp14=get_CN(p14,x)
CNp15=get_CN(p15,x)
CNp16=get_CN(p16,x)
CNp17=get_CN(p17,x)
CNp18=get_CN(p18,x)
CNp19=get_CN(p19,x)

CNx=get_CN(x,x)

#Cluster similar palindromes
C1<-c(CNp1,CNp6,CNp8,CNp10,CNp14,CNp16)
C2<-c(CNp2,CNp11,CNp15)
C3<-c(CNp3,CNp12)
C4<-c(CNp4,CNp13)
C5<-c(CNp5,CNp7,CNp9)

#Generate ggplot friendly matrix
out<-rbind(
cbind(CNx,rep("XDG",length(CNx)),species),
cbind(C1,rep("C1",length(C1)),species),
cbind(C2,rep("C2",length(C2)),species),
cbind(C3,rep("C3",length(C3)),species),
cbind(C4,rep("C4",length(C4)),species),
cbind(C5,rep("C5",length(C5)),species),
cbind(CNp17,rep("C17",length(CNp17)),species),
cbind(CNp18,rep("C18",length(CNp18)),species),
cbind(CNp19,rep("C19",length(CNp19)),species)
)

#Median number for Supplemtary table
round(median(C1[which(C1>0)]),digits=2)
round(median(C2[which(C2>0)]),digits=2)
round(median(C3[which(C3>0)]),digits=2)
round(median(C4[which(C4>0)]),digits=2)
round(median(C5[which(C5>0)]),digits=2)
round(median(CNp17[which(CNp17>0)]),digits=2)
round(median(CNp18[which(CNp18>0)]),digits=2)
round(median(CNp19[which(CNp19>0)]),digits=2)

bon_cov<-out

######################################ORANGUTAN

species='orangutan'
path=palindromeCN_Odir

#If file is empty does not print out error.
#Trying to overcome  missing data for few palindromes
palin=palin1=palin2=palin3=palin4=palin5=palin6=palin7=palin8=palin9=palin10=palin11=palin12=palin13=palin14=palin15=palin16=palin17=palin18=palin19=XDG=0

palin=read.table(paste(path,"Coverage_",species,"_msY_WGS_palindromeAll_80per.bed",sep=""), sep="\t", stringsAsFactors=F, header=F)
palin1=read.table(paste(path,"Coverage_",species,"_msY_WGS_palindrome1_80per.bed",sep=""), sep="\t", stringsAsFactors=F, header=F)
palin2=read.table(paste(path,"Coverage_",species,"_msY_WGS_palindrome2_80per.bed",sep=""), sep="\t", stringsAsFactors=F, header=F)
palin3=read.table(paste(path,"Coverage_",species,"_msY_WGS_palindrome3_80per.bed",sep=""), sep="\t", stringsAsFactors=F, header=F)
palin4=read.table(paste(path,"Coverage_",species,"_msY_WGS_palindrome4_80per.bed",sep=""), sep="\t", stringsAsFactors=F, header=F)
palin5=read.table(paste(path,"Coverage_",species,"_msY_WGS_palindrome5_80per.bed",sep=""), sep="\t", stringsAsFactors=F, header=F)
palin6=read.table(paste(path,"Coverage_",species,"_msY_WGS_palindrome6_80per.bed",sep=""), sep="\t", stringsAsFactors=F, header=F)
palin7=read.table(paste(path,"Coverage_",species,"_msY_WGS_palindrome7_80per.bed",sep=""), sep="\t", stringsAsFactors=F, header=F)
palin8=read.table(paste(path,"Coverage_",species,"_msY_WGS_palindrome8_80per.bed",sep=""), sep="\t", stringsAsFactors=F, header=F)
palin9=read.table(paste(path,"Coverage_",species,"_msY_WGS_palindrome9_80per.bed",sep=""), sep="\t", stringsAsFactors=F, header=F)
palin10=read.table(paste(path,"Coverage_",species,"_msY_WGS_palindrome10_80per.bed",sep=""), sep="\t", stringsAsFactors=F, header=F)
palin11=read.table(paste(path,"Coverage_",species,"_msY_WGS_palindrome11_80per.bed",sep=""), sep="\t", stringsAsFactors=F, header=F)
palin12=read.table(paste(path,"Coverage_",species,"_msY_WGS_palindrome12_80per.bed",sep=""), sep="\t", stringsAsFactors=F, header=F)
palin13=read.table(paste(path,"Coverage_",species,"_msY_WGS_palindrome13_80per.bed",sep=""), sep="\t", stringsAsFactors=F, header=F)
palin14=read.table(paste(path,"Coverage_",species,"_msY_WGS_palindrome14_80per.bed",sep=""), sep="\t", stringsAsFactors=F, header=F)
palin15=read.table(paste(path,"Coverage_",species,"_msY_WGS_palindrome15_80per.bed",sep=""), sep="\t", stringsAsFactors=F, header=F)
palin16=read.table(paste(path,"Coverage_",species,"_msY_WGS_palindrome16_80per.bed",sep=""), sep="\t", stringsAsFactors=F, header=F)
palin17=read.table(paste(path,"Coverage_",species,"_msY_WGS_palindrome17_80per.bed",sep=""), sep="\t", stringsAsFactors=F, header=F)
palin18=read.table(paste(path,"Coverage_",species,"_msY_WGS_palindrome18_80per.bed",sep=""), sep="\t", stringsAsFactors=F, header=F)
palin19=read.table(paste(path,"Coverage_",species,"_msY_WGS_palindrome19_80per.bed",sep=""), sep="\t", stringsAsFactors=F, header=F)
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
p9=get_coverage(palin9)
p10=get_coverage(palin10)
p11=get_coverage(palin11)
p12=get_coverage(palin12)
p13=get_coverage(palin13)
p14=get_coverage(palin14)
p15=get_coverage(palin15)
p16=get_coverage(palin16)
p17=get_coverage(palin17)
p18=get_coverage(palin18)
p19=get_coverage(palin19)

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
CNp9=get_CN(p9,x)
CNp10=get_CN(p10,x)
CNp11=get_CN(p11,x)
CNp12=get_CN(p12,x)
CNp13=get_CN(p13,x)
CNp14=get_CN(p14,x)
CNp15=get_CN(p15,x)
CNp16=get_CN(p16,x)
CNp17=get_CN(p17,x)
CNp18=get_CN(p18,x)
CNp19=get_CN(p19,x)

CNx=get_CN(x,x)

#Cluster similar palindromes
C1<-c(CNp1,CNp6,CNp8,CNp10,CNp14,CNp16)
C2<-c(CNp2,CNp11,CNp15)
C3<-c(CNp3,CNp12)
C4<-c(CNp4,CNp13)
C5<-c(CNp5,CNp7,CNp9)

#Generate ggplot friendly matrix
out<-rbind(
cbind(CNx,rep("XDG",length(CNx)),species),
cbind(C1,rep("C1",length(C1)),species),
cbind(C2,rep("C2",length(C2)),species),
cbind(C3,rep("C3",length(C3)),species),
cbind(C4,rep("C4",length(C4)),species),
cbind(C5,rep("C5",length(C5)),species),
cbind(CNp17,rep("C17",length(CNp17)),species),
cbind(CNp18,rep("C18",length(CNp18)),species),
cbind(CNp19,rep("C19",length(CNp19)),species)
)

#Median number for Supplemtary table
round(median(C1[which(C1>0)]),digits=2)
round(median(C2[which(C2>0)]),digits=2)
round(median(C3[which(C3>0)]),digits=2)
round(median(C4[which(C4>0)]),digits=2)
round(median(C5[which(C5>0)]),digits=2)
round(median(CNp17[which(CNp17>0)]),digits=2)
round(median(CNp18[which(CNp18>0)]),digits=2)
round(median(CNp19[which(CNp19>0)]),digits=2)

ora_cov<-out

##################################

#Colors.name <-  c(bo="#A7C855", ch="#74A8D3" , hu="6A6A6A", go="#BD3231", or="#E19F5A")
Colors.name <-  c(bo="#ABD9E9", ch="#2C7BB6" , hu="D7191C", go="#FFFFBF", or="#FDAE61")
darkColors.name <-  c(bo="#799630", ch="#3778ad" , hu="#4a4a4a", go="#842222", or="#ba6f21")


CN_chimp_palindrome<-rbind(gor_cov,bon_cov,ora_cov)
colnames(CN_chimp_palindrome)<-c("CopyNumber","Palindrome","Species")
CN_chimp_palindrome<-as.data.frame(CN_chimp_palindrome,stringsAsFactors = FALSE)
CN_chimp_palindrome$CopyNumber<-as.numeric(CN_chimp_palindrome$CopyNumber)
CN_chimp_palindrome$CopyNumber<-as.numeric(CN_chimp_palindrome$CopyNumber)
CN_chimp_palindrome$Species<-substr(CN_chimp_palindrome$Species,0,2)
CN_chimp_palindrome$Palindrome<- factor(CN_chimp_palindrome$Palindrome, levels = c("C1","C2","C3","C4","C5","C17","C18","C19","XDG"))
require(ggplot2)

CN_log_chimp_palindrome<-CN_chimp_palindrome
CN_log_chimp_palindrome$CopyNumber<-log(as.numeric(CN_log_chimp_palindrome$CopyNumber))

chimpLog_obj<-ggplot(CN_log_chimp_palindrome, aes(x=Species, y=CopyNumber))+
	geom_boxplot(aes(fill=Species),size=0.1,outlier.shape = NA, width=0.8)+scale_colour_manual(values=Colors.name)+ scale_fill_manual(values=Colors.name)+
	geom_hline(yintercept = 0, linetype = 2,size=0.1)+
	geom_hline(yintercept = 0.6931472, linetype = 2,color='red',size=0.1)+
	scale_y_continuous(limits=c(-0.75,3.218876),breaks=c(-0.6931472,0,0.6931472,1.609438,2.302585,2.70805,2.995732),labels=c("-0.6931472"="0.5","0" = "1","0.6931472" = "2", "1.609438"="5","2.302585"="10", "2.70805" = "15", "2.995732" = "20"))+
	theme_classic(base_size = 10) + 
	theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())+
	facet_wrap( ~ Palindrome, nrow = 1)

pdf("Fig_Chimp_copynumberLOG.pdf")
chimpLog_obj
dev.off()
	