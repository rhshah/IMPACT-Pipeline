#!/ifs/opt/bin/Rscript

###############################
## Create metrics documents for sequencing runs
## Includes Cluster Density and alignment rate, Base quality score, Insert Size Distribution,
##		Fingerprint matches, Contamination, Capture specificity, Duplication rate, Library Size,
##		Mean target coverage, and Mean exon coverage
## KRIS KANIA, August 2012
## Modified by A. Rose Brannon
## Modified by Ronak Shah, Jan 2013
##
## submit as
##	/home/brannona/software/Rbin/R-2.15.0/bin/Rscript AllMetrics_20120827.r "prefix='/path/to/experiment_prefix'" 
##  I.e. /home/brannona/software/Rbin/R-2.15.0/bin/Rscript AllMetrics_20120827.r "prefix='/home/user/analysis/ColonTrio11'" 
###############################

#Get prefix from command line call
args=(commandArgs(TRUE))
for(i in 1:length(args)){
	eval(parse(text=args[i]))
}

source("/home/brannona/software/textplot.R") #Fix textplot issues
library("gplots")
library("RColorBrewer")
#prefix="/Users/KrisKania/Workspace/fileprefix"
out_file= paste(prefix,"_allmetrics.pdf", sep="")
pdf(out_file, height=8.5, width=11)


##DATA TABLE
titlefile<-"_title.txt"
title<- paste(prefix,titlefile,sep="")
title <- read.delim(title)
title<-title[complete.cases(title),] #Removes cases with null values
barcode <- title$Sample_ID #create a list of the Sample ID's to replace BC's later 

#Set up to have both barcode and sample_ID
combinedname<-vector()
for(i in 1:nrow(title)){
combinedname<-append(combinedname,paste(title$Barcode[i]," :",barcode[i],sep=""))}


##Title Page
#2 objects on one page

layout(matrix(c(1,1,2,2),2,2,byrow=TRUE), 1, c(1,2.3))
plot(0:10, type = "n", xaxt="n", yaxt="n", bty="n", xlab = "", ylab = "")
text(11, 9, "Memorial Sloan Kettering",adj=1)
text(11, 6, "Berger Lab",adj=1)
text(11, 3, "Metric Script Version:1.2",adj=1)
text(11, 1, date(),adj=1)
text(1, 9, "Metric Analysis",adj=0,cex=1.5,font=2)
text(1, 6, paste("Project pool: ", title$Pool[1],sep=""),adj=0,cex=1.5,font=2) #Add project name from title file data
text(1, 3, paste("Baits version: ", title$Bait_version[1],sep=""),adj=0,cex=1.5,font=2) #Add baits version from title file data
#x<- tail(unlist(strsplit(prefix, "/")),1); text(4, 2, paste("Project: ",x,sep=""), adj=0,cex=1.5) #Add project name from file name.

#Set up title page data:
size=1.4
if(is.null(title$Sample_type)==FALSE){
	if(identical(title$Sample_ID,title$Collab_ID)==TRUE){
		trimtitle<-data.frame(as.character(title$Barcode), as.character(title$Sample_ID), as.character(title$Class), as.character(title$Sample_type),title$Input_ng, title$Library_yield,title$Pool_input)
		colnames(trimtitle)<-c("Barcode","Sample_ID","Class","Sample_type","Input_ng","Library_yield","Pool_input")
	} else{
		trimtitle<-data.frame(as.character(title$Barcode), as.character(title$Sample_ID), as.character(title$Collab_ID),as.character(title$Class), as.character(title$Sample_type), title$Input_ng, title$Library_yield,title$Pool_input)
		colnames(trimtitle)<-c("Barcode","Sample_ID","Collab_ID","Class","Sample_type","Input_ng","Library_yield","Pool_input")
		size=1.2
}}else {
	if(identical(title$Sample_ID,title$Collab_ID)==TRUE){
		trimtitle<-data.frame(as.character(title$Barcode), as.character(title$Sample_ID), as.character(title$Class), title$Input_ng, title$Library_yield,title$Pool_input)
		colnames(trimtitle)<-c("Barcode","Sample_ID","Class","Input_ng","Library_yield","Pool_input")
	} else{
		trimtitle<-data.frame(data.frame(title$Barcode), as.character(title$Sample_ID), as.character(title$Collab_ID),as.character(title$Class), title$Input_ng, title$Library_yield,title$Pool_input)
		colnames(trimtitle)<-c("Barcode","Sample_ID","Collab_ID","Class","Input_ng","Library_yield","Pool_input")
}}
textplot(trimtitle, show.rownames=FALSE,hadj=0,valign="top",cex=size,mar=c(0,0,0,0)) # print title information
par(mfrow=c(1,1))
###Plot Configration File
PCF<-"RunConfigration.txt"
base <- dirname(prefix)
x<-paste(base,PCF,sep="/")
xdata<-readLines(x)
par(mai=c(3,0.82,0.82,0.42))
textplot(xdata,show.rownames = FALSE, halign="left",valign="top")
title("Version of programs & files used")



#####CLUSTER DENSITY, ALIGNMENT RATE & Trimming Stats
hsmetricsfile<-"_ALL_HSmetrics.txt"
x<- paste(prefix,hsmetricsfile,sep="")
base <- read.delim(x)

Period<-combinedname #x-axis values

both_reads_align<-base$both.reads.align
one_read_aligns<-base$one.read.aligns
neither_read_aligns<-base$neither.read.aligns
Read1<-base$PerRead1Trimmed
Read2<-base$PerRead2Trimmed
#setting up data sets
col.topo<-colors()[c(556,142)] #color scheme


#Clipping Stats
plotdata <- rbind(Read1,Read2)
par(mai=c(3,0.82,0.82,0.42))
barplot(plotdata,beside=TRUE, main="Percentage of Reads Trimmed", ylab="(Percentage)",names=Period,las=3, cex.names=1, col=col.topo,legend=TRUE, args.legend = list(x="topright", bg= "white"))

col.topo<-colors()[c(556,142,17)] #color scheme
#Alignment Stats
par(mai=c(2,0.82,0.82,0.42))
v <- rbind(both_reads_align,one_read_aligns,neither_read_aligns)
EasyNumbers<-v/1000000
maximum<- max(colSums(EasyNumbers, na.rm=TRUE))
barplot(EasyNumbers,beside=FALSE, main="Cluster Density & Alignment Rate", ylab="(millions)", ylim=c(0,1.5*max(maximum)),names=Period,las=3, cex.names=1, col=col.topo,legend=TRUE, args.legend = list(x="topright", bg= "white"))

#Cluster density
Total_Clusters_Per_Lane<-(sum(v))/1000000
text(7,1.3*max(maximum),labels=paste("Total Clusters Per Lane=",round(Total_Clusters_Per_Lane),"Million",sep=" "))
Percent_With_Both_Reads_Aligning<-(sum(both_reads_align)/sum(v))*100
height<-c(1.2*max(maximum))
text(7,height,labels=paste("Percent With Both Reads Aligning=",round(Percent_With_Both_Reads_Aligning, digits = 2),"%",sep=""))



#Percentage graph of read alignment
data.perc <- apply(v, 2, function(x){x/sum(as.numeric(x))})
x<-min(data.perc[1,], na.rm=TRUE)

par(mai=c(3,0.82,0.82,0.42))
barplot(data.perc,beside=FALSE,ylim=c(x,1), main="Cluster Density & Alignment Rate", names=Period,col=col.topo,las=3, cex.names=1, lwd=2,legend=TRUE, args.legend = list(x="bottomright", bg="white"))


#### ORG BASE QUALITY SCORES
basequalityfile<-"_ALL_orgbasequalities.txt"
x<- paste(prefix,basequalityfile,sep="")
base <- read.delim(x)
qualityscore<-base[-nrow(base),-1] #qualityscore (y axis) will not refer to cycle or N/A row

cycle=seq(1:nrow(qualityscore)) #cycle (the x axis) will be rows 1-150
col.topo<-colors()[c(501,503,508,514,523,528,546,547,551,552,566,641,639,657,654,601,450,448,491,393,51,68,100,153)] #color catalog rainbow!

par(mai=c(.45,0.82,0.82,0.42))
plot(cycle, qualityscore[,1], main="Original Base Quality Score", xlab= "Cycle Number", type="l",ylim=c(0,max(qualityscore,na.rm=TRUE)),col=col.topo,ylab="Mean Quality Score",lwd=2) #plotting the data for just barcode 1
for(i in 2:ncol(qualityscore)) {lines(cycle,qualityscore[,i],col=col.topo[i],lwd=2)} #for looping the data for every barcode
legend(20,0.98*min(qualityscore)+3,legend=barcode,lty=rep(1,length(barcode)),lwd=rep(5,length(barcode)),col=col.topo[1:ncol(qualityscore)],cex=0.75) 

####BASE QUALITY SCORES
basequalityfile<-"_ALL_basequalities.txt"
x<- paste(prefix,basequalityfile,sep="")
base <- read.delim(x)
qualityscore<-base[-nrow(base),-1] #qualityscore (y axis) will not refer to cycle or N/A row

cycle=seq(1:nrow(qualityscore)) #cycle (the x axis) will be rows 1-150
col.topo<-colors()[c(501,503,508,514,523,528,546,547,551,552,566,641,639,657,654,601,450,448,491,393,51,68,100,153)] #color catalog rainbow!

par(mai=c(.45,0.82,0.82,0.42))
plot(cycle, qualityscore[,1], main="Recalibrated Base Quality Score", xlab= "Cycle Number", type="l",ylim=c(0,max(qualityscore,na.rm=TRUE)),col=col.topo,ylab="Mean Quality Score",lwd=2) #plotting the data for just barcode 1
for(i in 2:ncol(qualityscore)) {lines(cycle,qualityscore[,i],col=col.topo[i],lwd=2)} #for looping the data for every barcode
legend(20,0.98*min(qualityscore)+5,legend=barcode,lty=rep(1,length(barcode)),lwd=rep(5,length(barcode)),col=col.topo[1:ncol(qualityscore)],cex=0.75) 

###INSERT SIZE DISTRIBUTIONS
INSERTS<-"_ALL_insertsizemetrics.txt"
x<- paste(prefix,INSERTS,sep="")
insert <- read.delim(x)

barcodesize<-insert[-501,-1] #barcodesize (y axis) will not refer to the peak values or insert size column
size=seq(1:500)  #size (the x axis) will be rows 1-500
number<-ncol(barcodesize) #number refers to the count of how many barcodes (up to 24)


plot(size, barcodesize[,1], main="Insert Size Distribution", xlab= "Insert Size",type="l",ylim=c(0,max(barcodesize)),col=col.topo,ylab="",lwd=5)
#plotting the data for just barcode 1
for(i in 2:number) {lines(size,barcodesize[,i],col=col.topo[i],lwd=5)}
#for looping the data for every barcode

peak<-insert[501,-1] #peak refers to the peak row and does not include the insert size column

legendvalues<-paste(barcode[1],peak[1],sep=", peak=")
for(i in 2:ncol(barcodesize)){legendvalues<-append(legendvalues,paste(barcode[i],peak[i],sep=", peak="))}
#appending legend with correct color and peak value

legend(315,0.995*max(barcodesize),legend=legendvalues,lty=rep(1,length(barcode)),lwd=rep(2,length(barcode)),col=col.topo[1:ncol(barcodesize)],cex=0.75)#adding legend to graph

rm(legendvalues)
#removing the legend


###FP MATCHING
if(FALSE)
{
fingerprint<-"_ALL_FPsummary.txt"
x<- paste(prefix,fingerprint,sep="") #creates master fingerprint file

##Retrieve only the genotypes
dat<-read.delim(x)
dat1=dat[1:1034,grep("Genotypes",colnames(dat))]

names(dat1)[1:ncol(dat1)]<-paste(barcode) #how to replace dat1 column headings (barcode) with sample name

sink(paste(unlist(strsplit(x,".",fixed=TRUE))[1],".matchnames.txt",sep=""))
cat("Fingerprint matches require at least 900 of the 1033 fingerprint sites to match.\n")
for (i in 1:(ncol(dat1)-1)){
  for (j in i+1:(ncol(dat1)-i)){
    num.matches=length(which(as.vector(dat1[,i])==as.vector(dat1[,j])))
    ###Can change 38 below to 0.9(nrow(dat1))
    if(num.matches>=900){ 
      cat(paste(colnames(dat1)[i]," matches ",colnames(dat1)[j],": ",num.matches,"/1033.\n",sep=""))
    }
  }
}
sink()

match<- "_ALL_FPsummary.matchnames.txt"
x<- paste(prefix,match,sep="")
dat<-read.delim(x,check.names=FALSE)
textplot(dat, show.rownames = FALSE, hadj=0)
}

###Mixups

FPX<-"_ALL_FPCsummary.txt"
x<- paste(prefix,FPX,sep="")
xt<-read.table(x, as.is = TRUE, header = TRUE, sep = "\t", row.names = 1)
xm<-as.matrix(xt)
par(mai=c(3,0.82,0.82,0.42))
Period<-combinedname
main_title="Sample Mix-Ups"
#pro = diff(range(xm))
val = isSymmetric(xm)
if(val == FALSE)
{
	heatmap.2(xm,Rowv="NA",Colv="NA",dendrogram="none",main=main_title,distfun="dist",hclustfun="hclust",labRow=Period,labCol=Period,key="TRUE",keysize=1,trace="none",density.info=c("none"),margins=c(20,20),col=brewer.pal(5,"Spectral"))
	mtext("Note: The Value below the key refers to the fraction of discordant homozygous alleles",side=1,line=10)
}
	
###Mix-ups Results

FPR<-"_ALL_FPCResults.txt"
x<- paste(prefix,FPR,sep="")
par(mai=c(3,0.82,0.82,0.42))
xdata<-readLines(x)
textplot(xdata, show.rownames = FALSE, halign="left",valign="top")
title("Fraction of positions called hetrozygous")

###Major Contamination

FPH<-"_ALL_FPhet.txt"
x<- paste(prefix,FPH,sep="")

skip <- read.delim(x)
allele<-skip$PerHetrozygosPos

par(mai=c(3,0.82,0.82,0.42))
maximum<-1.5*(as.numeric(max(allele, na.rm=TRUE)))
Period<-combinedname

barplot(t(allele), main="Major Contamination Check", ylab="Fraction of position that are hetrozygous",names=Period, ylim=c(0,maximum), cex.names=1, lwd=2,las=3, col="black")
abline(h=c(.5),col="Red")

m<-mean(skip$PerHetrozygosPos, na.rm=TRUE)
maximum<-1.25*(as.numeric(max(allele, na.rm=TRUE)))
text(5,maximum,labels=paste("Fraction of position that are hetrozygous=",round(m, digits=5),sep=""))


###Minor CONTAMINATION
FP<-"_ALL_FPavgHom.txt"
x<- paste(prefix,FP,sep="")

skip <- read.delim(x)
allele<-skip$AvgMinorHomFreq

par(mai=c(3,0.82,0.82,0.42))
maximum<-1.5*(as.numeric(max(allele, na.rm=TRUE)))
Period<-combinedname

barplot(t(allele), main="Minor Contamination Check", ylab="Avg. Minor Allele Frequency at Homozygous Position",names=Period, ylim=c(0,maximum), cex.names=1, lwd=2,las=3, col="black")
abline(h=c(.01),col="Red")

m<-mean(skip$AvgMinorHomFreq, na.rm=TRUE)
maximum<-1.25*(as.numeric(max(allele, na.rm=TRUE)))
text(5,maximum,labels=paste("Average Minor Allele Frequency=",round(m, digits=5),sep=""))

####GC Plot
GCP<-"_ALL_gcbias.txt"
x<- paste(prefix,GCP,sep="")
xt<-read.delim(x)

barcodesize<-xt[,-1] #barcodesize (y axis) will not refer to the peak values or insert size column
gc=xt[,1]  #size (the x axis) will be rows 1-500
number<-ncol(barcodesize) #number refers to the count of how many barcodes (up to 24)


plot(gc,barcodesize[,1], main="Normaized Coverage vs GC-Content", xlab="GC-Content",ylab="Normalized Coverage",type="l",xlim=c(min(gc),max(gc)+0.35),ylim=c(0,max(barcodesize)),col=col.topo,lwd=5) 
#plotting the data for just barcode 1
#title(ylab="Median of Normalized Coverage")
for(i in 2:number) {lines(gc,barcodesize[,i],col=col.topo[i],lwd=5)} #for looping the data for every barcode

legendvalues<-paste(barcode[1],sep='')
for(i in 2:ncol(barcodesize)){legendvalues<-append(legendvalues,paste(barcode[i]))}
#appending legend with correct color and peak value

legend(max(gc)+0.05,0.995*max(barcodesize)+0.1,legend=legendvalues,lty=rep(1,length(barcode)),lwd=rep(2,length(barcode)),col=col.topo[1:ncol(barcodesize)],cex=0.75)#adding legend to graph

rm(legendvalues)
#removing the legend

###CAPTURE SPECIFICITY
hsmetricsfile<-"_ALL_HSmetrics.txt"
x<- paste(prefix,hsmetricsfile,sep="")
base <- read.delim(x)


ON_BAIT_BASES<-base$ON_BAIT_BASES
NEAR_BAIT_BASES<-base$NEAR_BAIT_BASES
OFF_BAIT_BASES<-base$OFF_BAIT_BASES
#setting up data sets

col.topo<-colors()[c(556,142,17)] #color scheme

par(mai=c(3,0.82,0.82,0.42))
v <- rbind(ON_BAIT_BASES,NEAR_BAIT_BASES,OFF_BAIT_BASES)
EasyNumbers1<-v/1000000

maximum<- max(colSums(EasyNumbers1,na.rm=TRUE))
barplot(EasyNumbers1,beside=FALSE,main="Capture Specificity", ylim=c(0,(1.5*max(colSums(EasyNumbers1)))), names=Period, cex.names=1, col=col.topo,ylab="Total Bases (millions)",las=3, lwd=2,legend=TRUE, args.legend = list(x="topright", bg="white"))
#graphing barplot and legend.. how to adjust legend?


u<-sum(as.numeric(ON_BAIT_BASES)) #summation of on bait
t<-sum(as.numeric(NEAR_BAIT_BASES)) #summation of near beait
x<-(u+t) #summation of on and near
l<-sum(as.numeric(v)) #summation of on, off, near
k<-(x/l)*100 #percent on/near

text(7,1.4*maximum, labels=paste("Average % Selected On/Near Bait=",round(k, digits = 1),"%",sep=""))

f<-(u/l)*100 #percent on
text(7,1.3*maximum, labels=paste("Average % On Bait=",round(f, digits = 1),"%",sep=""))

ON_TARGET_BASES<-base$ON_TARGET_BASES #calling on target bases
s<-sum(as.numeric(ON_TARGET_BASES)) #summation of on target bases

q<-(s/l)*100 #percent on target/(on/near/off)
text(7,1.2*maximum,labels=paste("Average % On Target=",round(q, digits = 1),"%",sep=""))


#Changing the specificity to be based on percents instead.
data.perc <- apply(v, 2, function(x){x/sum(as.numeric(x))})
x<-min(data.perc[1,], na.rm=TRUE)
par(mai=c(3,0.82,0.82,0.42))
barplot(data.perc,beside=FALSE,ylim=c((0.9*x),1), main="Capture Specificity", names=Period,cex.names=1, col=col.topo,las=3, lwd=2,legend=TRUE, args.legend = list(x="bottomright", bg="white"))


###LIBRARY COMPLEXITY
hsmetricsfile<-"_ALL_HSmetrics.txt"
x<- paste(prefix,hsmetricsfile,sep="")
base <- read.delim(x)

duplication<-base$PERCENT_DUPLICATION #duplication will refer to the values in percent duplication rate Column
par(mai=c(3,0.82,0.82,0.42))
barplot(duplication, main="Estimated Duplication Rate", ylim=c(0,1.1), names.arg=Period, las=3, cex.names=1, col="black")
text(5,1,labels=paste("Average Percent Duplication Rate=",round(mean(duplication, na.rm=TRUE)*100),"%",sep=""))

complexity<-base$ESTIMATED_LIBRARY_SIZE #complexity will refer to the values in Estimated Library Size Column
par(mai=c(3,0.82,0.82,0.42))
EasyNumbers2<-complexity/1000000
barplot(EasyNumbers2, main="Estimated Library Size", ylim=c(0,1.5*max(EasyNumbers2, na.rm=TRUE)), ylab="(millions)", names.arg=Period, las=3, cex.names=1, col="darkred")
text(5,1.1*(max(EasyNumbers2,na.rm=TRUE)),labels=paste("Average Library Size=",round((mean(complexity, na.rm=TRUE)/1000000), digit=1), "million", sep=" "))
total<-(sum(complexity,na.rm=TRUE))/1000000
text(5,1.2*(max(EasyNumbers2,na.rm=TRUE)),labels=paste("Total Library Size=",round(total),"million", sep=" "))


###MEAN TARGET COVERAGE
coverage<-base$MEAN_TARGET_COVERAGE #coverage will refer to the values in Mean Target Coverage Column
par(mai=c(3,0.82,0.82,0.42))
barplot(coverage, main="Mean Target Coverage", ylim=c(0,1.8*max(coverage,na.rm=TRUE)), names.arg=Period, las=3, cex.names=1, col="black", ylab="Mean Target Coverage")
abline(h=100,col="Red")
abline(h=400,col="forestgreen")
text(4,1.6*max(coverage),labels=paste("Average Coverage=",round(mean(coverage,na.rm=TRUE),0),"X",sep=" "))

###average over normals
i<-which(title$Class=="Normal")
normalcoverage<-coverage[i]
normaldenominator<-length(normalcoverage)
text(4,1.4*max(coverage),labels=paste("Average Normal Coverage=",round(mean(normalcoverage,na.rm=TRUE),0),"X",sep=" "))

#average over tumors
i<-which(title$Class=="Tumor")
tumorcoverage<-coverage[i]
tumordenominator<-length(tumorcoverage)
text(4,1.2*max(coverage),labels=paste("Average Tumor Coverage=",round(sum(tumorcoverage,na.rm=TRUE)/tumordenominator,0),"X",sep=" "))


##MEAN EXON COVERAGE
#exonfile<-"_ALL_exoncoverage.txt"
#x<- paste(prefix,exonfile,sep="")
#base <- read.delim(x) #reading the file
#base<-base[,-(1:2)]
#exon<-colMeans(base)
#par(mai=c(3,0.82,0.82,0.42))
#barplot(exon, main="Mean Exon Coverage", ylim=c(0,1.8*max(exon,na.rm=TRUE)), names.arg=Period, las=3, cex.names=1, col="darkred", ylab="Mean Exon Coverage")
#abline(h=100,col="Red")
#abline(h=400,col="forestgreen")
#denominator <- length(exon)
#text(6,1.6*max(exon),labels=paste("Average Exon Coverage=",round(sum(exon,na.rm=TRUE)/denominator,0),"X",sep=" ")) #average over exons
#
##average over normal exons
#i<-which(title$Class=="Normal")
#normalexoncoverage<-exon[i]
#text(6,1.4*max(exon),labels=paste("Average Exon Coverage (Normals) =",round(mean(normalexoncoverage,na.rm=TRUE),0),"X",sep=" "))
#
##average over tumor exons
#i<-which(title$Class=="Tumor")
#tumorexoncoverage<-exon[i]
#text(6,1.2*max(exon),labels=paste("Average Exon Coverage (Tumors) =",round(mean(tumorexoncoverage,na.rm=TRUE),0),"X",sep=" "))
#

dev.off() 
#print(paste("Plot was saved in:", getwd())) 




