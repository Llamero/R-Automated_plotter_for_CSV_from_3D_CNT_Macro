#Set working directory to file directory with data
setwd("E:/Grady Lab/Hessam/140623/All sample statistics for automated analysis")

#Create a subdirectory within the working directory to store the final plots
dir.create("Standard Plots/")

#Load 3D aggregate information into matrix
read.csv("CNT 3D Quantitation.csv")->PercentCNTMatrix

#This function is to generate the box/line plot of CNT distribution in X, y, and Z axes
DataPlotter<-function(CNTvector) {

#----------------Retrieve the default info for the sample and get the volume and position statistics generated in Imaris-------------------------------------------------------------
	#Get the sample name, number of bins, and axis limits output from the macro
	Sample<-CNTvector[2]
	nBin<-as.numeric(CNTvector[3])
	Xmin<-as.numeric(CNTvector[4])
	Xmax<-as.numeric(CNTvector[5])
	Ymin<-as.numeric(CNTvector[6])
	Ymax<-as.numeric(CNTvector[7])
	Zmin<-as.numeric(CNTvector[8])
	Zmax<-as.numeric(CNTvector[9])
	
	#Open the bundles position and volume statistics csv file for this sample
	#Add on file name specific to statistics csv for each sample
	#NOTE: sub removes .tif from end of sample name: REGEX code:
		#\\. - "." is a metacharacter that needs to be escaped.  R needs a double escape, rather than the conventional single escape "\" in regex
		# Therefore the regex code is looking for ".tif" in the name and then replacing it with nothing "", making a new string iwthout the ".tif"
	filePatternVolume<-paste(sub("\\.tif", "", Sample), "_Volume.csv", sep = "", collapse = "") 
	volumeFile<-list.files(path = ".", pattern = filePatternVolume, recursive = TRUE, ignore.case = TRUE, include.dirs = FALSE, no.. = FALSE)
	
	#NOTE:Function will crash here if CSV is empty, check to make sure *.csv file contains data!
	#Skip to line 3 - this removes excess rows and starts the read at the headers
	read.csv(volumeFile, skip = 3)->volumeMatrix

	filePatternPosition<-paste(sub("\\.tif", "", Sample), "_Position.csv", sep = "", collapse = "") 
	positionFile<-list.files(path = ".", pattern = filePatternPosition, recursive = TRUE, ignore.case = TRUE, include.dirs = FALSE, no.. = FALSE)

	#Skip to line 3 - this removes excess rows and starts the read at the headers
	read.csv(positionFile, skip = 3)->positionMatrix
	
	
	#---------------------Calculate bundle volume distribution per bin along each axis----------------------------------------------------------------------------------------------------------------------

	#Calculate the X bin size based on the same parameters used for measuring bulk CNT content
	XbinSize<-((Xmax-Xmin)/nBin)
	#Assemble a vector of bin increments
	Xbinseq<-seq(Xmin, Xmax, XbinSize)

	#Calculate the Y bin size based on the same parameters used for measuring bulk CNT content
	YbinSize<-((Ymax-Ymin)/nBin)
	#Assemble a vector of bin increments
	Ybinseq<-seq(Ymin, Ymax, YbinSize)
	
	#Calculate the Z bin size based on the same parameters used for measuring bulk CNT content
	ZbinSize<-((Zmax-Zmin)/nBin)
	#Assemble a vector of bin increments
	Zbinseq<-seq(Zmin, Zmax, ZbinSize)
		
	#Assign bin # to X, Y, and Z dimesnions
	#NOTE: If "check name" is not used when table is loaded, then spaces are replaced with periods, so "Poaition X" becomes "Position.X"
	XDataBin<-findInterval(positionMatrix$Position.X, Xbinseq)
	YDataBin<-findInterval(positionMatrix$Position.Y, Ybinseq)
	ZDataBin<-findInterval(positionMatrix$Position.Z, Zbinseq)

	#Build matrix of empty values to fill with sorted volume data by each axis - rows is equal to max number of aggregates within dataset
	XSortedData<-matrix(data=NA,nrow=length(volumeMatrix[,1]),ncol=(nBin + 1))
	#X axis - Sort aggregate volume data into column bins
	for (i in 1:length(XDataBin)){
		XSortedData[i,XDataBin[i]]<-volumeMatrix$Value[i]
	}
	#Remove the excess column, which contains last sliver of image volume, in order to keep all bin volumes identical
	XSortedData<-XSortedData[,-(nBin+1)]
	
	#Build matrix of empty values to fill with sorted volume data by each axis - rows is equal to max number of aggregates within dataset
	YSortedData<-matrix(data=NA,nrow=length(volumeMatrix[,1]),ncol=(nBin + 1))
	#Y axis - Sort aggregate volume data into column bins
	for (i in 1:length(YDataBin)){
		YSortedData[i,YDataBin[i]]<-volumeMatrix$Value[i]
	}
	#Remove the excess column, which contains last sliver of image volume, in order to keep all bin volumes identical
	YSortedData<-YSortedData[,-(nBin+1)]
	
	#Build matrix of empty values to fill with sorted volume data by each axis - rows is equal to max number of aggregates within dataset
	ZSortedData<-matrix(data=NA,nrow=length(volumeMatrix[,1]),ncol=(nBin + 1))
	#Z axis - Sort aggregate volume data into column bins
	for (i in 1:length(ZDataBin)){
		ZSortedData[i,ZDataBin[i]]<-volumeMatrix$Value[i]
	}
	#Remove the excess column, which contains last sliver of image volume, in order to keep all bin volumes identical
	ZSortedData<-ZSortedData[,-(nBin+1)]
	
	
	
	#Load measured total aggregate volume for each bin in each dimension from csv file
	#NOTE: Right index is to +9 rather than +10, because start point is effectively at 10+0  therefore end would be 10+nBin-1 which equates to nBin + 9
	Xpercenttotal<-as.numeric(CNTvector[10:(9+nBin)])
	Ypercenttotal<-as.numeric(CNTvector[(10+nBin):(2*nBin+9)])
	Zpercenttotal<-as.numeric(CNTvector[(2*nBin+10):(3*nBin+9)])
	
#--------------------------Generate final plots of the results---------------------------------------------------------------------------------------------------
#Offset and round binseq vectors (used as X axis label) so that labels start at 0 and are integers
Xbinseq<-round(Xbinseq-Xmin)
Ybinseq<-round(Ybinseq-Ymin)
Zbinseq<-round(Zbinseq-Zmin)

#Extract only sample name from full sample ID returned from ImageJ macro
#NOTE: This is specific for a given set of sample IDs.  Therefore, this will have to be modified per dataset naming pattern
#This REGEX epxression removes the leading iteration count and date in the name, such as "10_140624 "
SampleName<-sub("^[0-9]*_[0-9]* ", "", Sample)

#This REGEX expression then removes everything after the first " - " which records additional image info, not sample info, such as " - Sample 2 redo - Substack 8-214 - 8 bit.tif"
SampleName<-sub(" - .*", "", SampleName)

#Create a subdirectory within the "plots" subdirectory with just the dispersion method name (not sample name)
SubDirName<-paste("Standard Plots/", sub("~[0-9]*", "", SampleName))
dir.create(SubDirName, showWarnings = TRUE)

#NOTE: Ylim for log axis on ALL plots needs to be specified, because on some the notch is outside the box
#X axis plot
#Prompt for plot to be saved as tiff with sample name in the appropriate dispersion method subdirectory
ImageName<-paste(SubDirName, "/", SampleName, " - X Distribution.tiff", sep = "")
tiff(filename = ImageName)

par(mar = c(5,6,3,5))
boxplot(XSortedData, log="y", xaxt = "n",  xlab = "", notch = FALSE, ylim = c(min(XSortedData, na.rm = TRUE), max(XSortedData, na.rm = TRUE)), cex.lab=2, cex.axis=1.5, cex.main=2, cex.sub=2, bty="n")
axis(1, at=seq(0.5, (nBin+0.5), 1), labels=Xbinseq, cex.axis = 1.5)
par(new=T)
plot(Xpercenttotal, type="o", axes = F, xlab = NA, ylab = NA, pch = 2, cex=1.5, xlim = c((0.5),(nBin+0.5)), ylim = c(0, round((max(Xpercenttotal)+0.005), 2)), lty=5, lwd=2, bty="n")
axis(4, cex.axis=1.5)
mtext(side=1, line = 3, "X Axis Position (μm)", cex = 2)
mtext(expression(paste( plain("Bundle Volume Distribution (μm") ^ plain("3"), plain(")") )), side=2, line = 3, cex = 2)
mtext(side=4, line = 3, "CNT Volume : Sample Volume", cex = 2)
title(main = "X Axis CNT Bundle Distribution", cex.main = 1.8, line = 1.5)

#dev.off both prompts for the plot to be closed, as well as prompts for tiff to save an image of everything between the tiff prompt and dev.off
dev.off() 

#Y axis plot
#Prompt for plot to be saved as tiff with sample name in the appropriate dispersion method subdirectory
ImageName<-paste(SubDirName, "/", SampleName, " - Y Distribution.tiff", sep = "")
tiff(filename = ImageName)

par(mar = c(5,6,3,5))
boxplot(YSortedData, log="y", xaxt = "n",  xlab = "", notch = FALSE, ylim = c(min(XSortedData, na.rm = TRUE), max(XSortedData, na.rm = TRUE)), cex.lab=2, cex.axis=1.5, cex.main=2, cex.sub=2, bty="n")
axis(1, at=seq(0.5, (nBin+0.5), 1), labels=Ybinseq, cex.axis = 1.5)
par(new=T)
plot(Ypercenttotal, type="o", axes = F, xlab = NA, ylab = NA, pch = 2, cex = 1.5, xlim = c((0.5),(nBin+0.5)), ylim = c(0, round((max(Ypercenttotal)+0.005), 2)), lty=5, lwd=2, bty="n")
axis(4, cex.axis=1.5)
mtext(side=1, line = 3, "Y Axis Position (μm)", cex = 2)
mtext(expression(paste( plain("Bundle Volume Distribution (μm") ^ plain("3"), plain(")") )), side=2, line = 3, cex = 2)
mtext(side=4, line = 3, "CNT Volume : Sample Volume", cex = 2)
title(main = "Y Axis CNT Bundle Distribution", cex.main = 1.8, line = 1.5)

#dev.off both prompts for the plot to be closed, as well as prompts for tiff to save an image of everything between the tiff prompt and dev.off
dev.off()

#Z axis plot
#Prompt for plot to be saved as tiff with sample name in the appropriate dispersion method subdirectory
ImageName<-paste(SubDirName, "/", SampleName, " - Z Distribution.tiff", sep = "")
tiff(filename = ImageName)

par(mar = c(5,6,3,5))
boxplot(ZSortedData, log="y", xaxt = "n",  xlab = "", notch = FALSE, ylim = c(min(XSortedData, na.rm = TRUE), max(XSortedData, na.rm = TRUE)), cex.lab=2, cex.axis=1.5, cex.main=2, cex.sub=2, bty="n")
axis(1, at=seq(0.5, (nBin+0.5), 1), labels=Zbinseq, cex.axis = 1.5)
par(new=T)
plot(Zpercenttotal, type="o", axes = F, xlab = NA, ylab = NA, pch = 2, cex = 1.5, xlim = c((0.5),(nBin+0.5)), ylim = c(0, round((max(Zpercenttotal)+0.005), 2)), lty=5, lwd=2, bty="n")
axis(4, cex.axis=1.5)
mtext(side=1, line = 3, "Z Axis Position (μm)", cex = 2)
mtext(expression(paste( plain("Bundle Volume Distribution (μm") ^ plain("3"), plain(")") )), side=2, line = 3, cex = 2)
mtext(side=4, line = 3, "CNT Volume : Sample Volume", cex = 2)
title(main = paste(SampleName, "Z Axis CNT Bundle Distribution"), cex.main = 1.8, line = 1.5)	

#dev.off both prompts for the plot to be closed, as well as prompts for tiff to save an image of everything between the tiff prompt and dev.off
dev.off()

}

apply(PercentCNTMatrix, 1, DataPlotter)

