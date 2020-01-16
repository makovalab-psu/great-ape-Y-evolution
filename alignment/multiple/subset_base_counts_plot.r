plotFilename   = "apes_Y.subset_base_counts.pdf"
doOneClassBars = T
includeTitle   = F

speciesNames = c("human","chimpanzee","bonobo","gorilla","orangutan")
speciesNicks = c("human","chimp",     "bonobo","gorilla","orang")
speciesToColor        = c("#D7191C","#2C7BB6","#ABD9E9","#FFFFBF","#FDAE61")
names(speciesToColor) = speciesNames
names(speciesNames)   = speciesNicks

hatchDensity = 30

dd = read.table("fiveY.Anc0_centric.subset_base_counts",header=T,comment.ch="")
colnames(dd)[1] = "class"

ddOneToOne  = dd[dd$class=="onetoone",]
ddUnique    = dd[dd$class=="unique",]
ddDuplicate = dd[dd$class=="duplicate",]
if (nrow(ddOneToOne) != nrow(ddUnique))    print("HEY, SOMETHING'S WRONG!!!!!")
if (nrow(ddOneToOne) != nrow(ddDuplicate)) print("HEY, SOMETHING'S WRONG!!!!!")

speciesColumns = 2:ncol(dd)
species = speciesNames[colnames(dd)[speciesColumns]]
numSpecies = length(species)
colors = speciesToColor[species]

numBars = 0
for (subsetIx in 1:nrow(ddOneToOne))
	{
	subsetSize = sum(!is.na(ddOneToOne[subsetIx,speciesColumns]))
	numBars = numBars + 1 + subsetSize
	}

bbNames  = rep(NA,numBars)
bbColors = rep(NA,numBars)
ooCounts = rep(NA,numBars)
uuCounts = rep(NA,numBars)
ddCounts = rep(NA,numBars)
barIx = 1
for (subsetIx in 1:nrow(ddOneToOne))
	{
	ooSubsetCounts = ddOneToOne [subsetIx,speciesColumns]
	uuSubsetCounts = ddUnique   [subsetIx,speciesColumns]
	ddSubsetCounts = ddDuplicate[subsetIx,speciesColumns]
	subsetActive   = (!is.na(ooSubsetCounts))
	subsetSize     = sum(subsetActive)
	names          = colnames(ooSubsetCounts)[subsetActive]
	countsForOO    = ooSubsetCounts[subsetActive]
	countsForUU    = uuSubsetCounts[subsetActive]
	countsForDD    = ddSubsetCounts[subsetActive]
	nums           = (speciesColumns-1)[subsetActive]
	barIx = barIx+1
	for (ix in 1:subsetSize)
		{
		bbNames [barIx] = names[ix]
		bbColors[barIx] = colors[nums[ix]]
		ooCounts[barIx] = countsForOO[ix]
		uuCounts[barIx] = countsForUU[ix]
		ddCounts[barIx] = countsForDD[ix]
		barIx = barIx+1
		}
	}
maxCount = max(uuCounts+ddCounts,na.rm=T)
minCount = min(uuCounts+ddCounts,na.rm=T)
bbLabels = bbNames
bbLabels[is.na(bbNames)] = "-----"
bbLabels[1] = ""

turnDeviceOff = F
if (is.null(plotFilename)) {
	quartz(width=12.365,height=8.66)
} else {
	print(paste("drawing to",plotFilename))
	pdf(file=plotFilename,width=12.365,height=8.66,pointsize=9)
	turnDeviceOff = T
	}

extraBars = 70   # (width of the legend)
par(mar=c(3,6.5,1,0)+0.1)     # BLTR
options(scipen=10)
xLim = c(0,numBars+extraBars)
yLim = c(0,maxCount)

plotTitle = ""
if (includeTitle)
	plotTitle = "species subsets in ape Y progressiveCactus alignments"

if (doOneClassBars)
	{
	barPos = barplot(uuCounts+ddCounts,
	           xlim=xLim,ylim=yLim,las=1,
	           col=bbColors,border="black",
	           main=plotTitle,
	           xaxt="n")
} else {
	barPos = barplot(uuCounts,
	           xlim=xLim,ylim=yLim,las=1,
	           col=bbColors,border="black",
	           main=plotTitle,
	           xaxt="n")
	barplot(uuCounts,add=T,xlim=xLim,ylim=yLim,col=bbColors,density=hatchDensity,border="black",yaxt='n')
	barplot(ooCounts,add=T,xlim=xLim,ylim=yLim,col=bbColors,border="black",yaxt='n')
	barplot(uuCounts+ddCounts,add=T,xlim=xLim,ylim=yLim,col=NA,border="black",yaxt='n')
	}
title(ylab="aligned bases",mgp=c(4.5,1,0))
axis(1,at=barPos,labels=bbLabels,tick=F,las=2,cex.axis=0.5)
dividers = barPos[is.na(bbNames)]
dividers = c(dividers,2*barPos[length(barPos)]-barPos[length(barPos)-1])
for (x in dividers)
	{
	lines(c(x,x),c(0,maxCount*0.99),lty=3)
	}

if (doOneClassBars)
	{
	legend("topright",bg="white",legend=species,
	       pch=22,col="black",pt.bg=colors)
} else {
	legNames   = rep(NA,3*numSpecies)
	legColors  = rep(NA,3*numSpecies)
	legDensity = rep(NA,3*numSpecies)
	for (ix in 1:numSpecies)
		{
		legNames  [             ix] = paste(species[ix],", one-to-one",sep="")
		legNames  [  numSpecies+ix] = paste(species[ix],", unique minus one-to-one",sep="")
		legNames  [2*numSpecies+ix] = paste(species[ix],", duplicate",sep="")
		legColors [             ix] = colors[ix]
		legColors [  numSpecies+ix] = colors[ix]
		legColors [2*numSpecies+ix] = "white"
		legDensity[  numSpecies+ix] = hatchDensity
		}
	# note that this legend has to be replaced, manually
	legend("topright",bg="white",legend=legNames,
	       pch=22,col="black",pt.bg=legColors,density=legDensity)
	}

if (turnDeviceOff) dev.off()
