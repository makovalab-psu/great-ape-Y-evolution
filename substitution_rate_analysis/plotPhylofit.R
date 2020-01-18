setwd("/Users/polly/Desktop/projects/phylofit")
files<-list.files(path = ".", pattern = "*.mod")

for (filename in files) {
  Ytree<-read.tree(file=filename)
  tm<-read.tm(filename)
  
  rates<-tm[["rate.matrix"]]
  diag(rates)<-NA
  
  print(paste(basename(filename),sum(as.vector(rates),na.rm=TRUE)))
  species_col<-c("black","chartreuse3","blue","red","orange")
  
  tiff(file=paste0(basename(filename),".phylofit",".tiff"),units="cm", width=8.9, height=8.9, res=300,compression="rle")
  plot(Ytree,main=basename(filename),cex=1.5, lwd = 5, tip.color=species_col, show.tip.label = FALSE)
  #edgelabels(round(Ytree$edge.length,digits=4),frame="none",adj = c(0.5, -0.25))
  #tiplabels(Ytree$tip.label,frame="none",pch=1.5,col=species_col,cex=1.2, adj=-0.04)
  add.scale.bar(0.001,0.95,length=0.01,cex=1)
  dev.off()

}
