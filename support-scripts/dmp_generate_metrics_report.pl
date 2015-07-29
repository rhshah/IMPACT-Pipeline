use strict;
use warnings;

my $sample = $ARGV[0];
my $R = $ARGV[1];

my $output = $sample."_All_Metrics.rhtml";
open RSCRIPT, ">", $output || $!;
print RSCRIPT <<ENDSCRIPT;
<html>

<head>
<title>All Metrics</title>
<link href="/resources/jquery-plugins/style/bootstrap.css" rel="stylesheet" media="screen">
<script type="text/javascript" src="http://ajax.googleapis.com/ajax/libs/jquery/1.9.1/jquery.min.js"></script>  
<script class="jsbin" src="http://datatables.net/download/build/jquery.dataTables.nightly.js"></script>
<style>
.sidebar-nav-fixed {
     position:fixed;
     top:40px;
     width:15%;
 }
 </style>
 
<script type="text/javascript" charset="utf-8">
    \$(document).ready(function() {
        \$('#titleFile').dataTable({
        "bPaginate": false,
        "bLengthChange": false,
        "bFilter": true,
        "bSort": true,
        "bInfo": false,
        "bAutoWidth": true
        });
    } );
</script>
</script>
</head>

<body>
<!--begin.rcode include=FALSE
library(knitr)
library(plyr)
library(scales)
library(RColorBrewer)
library(gplots)
library(googleVis)
print(installed.packages())
my_theme <- function(base_size = 16, base_family="Helvetica", dkpanel=FALSE, legendPosition = "bottom", xAngle = 45, hJust = 1) {
  thme <- theme(
      line = element_line(colour = "#838281", size = 0.5,
                          linetype = 1, lineend = "round"),
      rect = element_rect(fill = "white",
                          colour = NA, size = 0.5, linetype = 1),
      text = element_text(family = base_family, face = "bold",
                          colour = "black", size = base_size, hjust = 0.5, vjust = 0.5,
                          angle = 0, lineheight = 1),
      axis.text = element_text(size = rel(1), face = "bold.italic"),
      axis.text.x = element_text(vjust = 1, hjust= hJust,angle = xAngle),
      axis.text.y = element_text(hjust = 1),
      axis.line = element_line(colour = "#d8d6d5"),
      axis.line.y = element_line(),
      axis.line.x = element_line(),
      axis.ticks = element_line(),
      axis.title = element_text(size = rel(1)),
      axis.title.x = element_text(),
      axis.title.y = element_text(angle = 90),
 
      axis.ticks.length = unit(base_size * 0.2 , "points"),
      axis.ticks.margin = unit(base_size * 0.5, "points"),
      
      legend.background = element_rect(color = "grey", linetype=1, size = rel(0.5)),
      legend.margin = unit(base_size * 2.5, "points"),
      legend.key = element_rect(),
      legend.key.size = unit(2, "lines"),
      legend.key.height = NULL,
      legend.key.width = NULL,
      legend.text = element_text(size = rel(1.1)),
      legend.text.align = NULL,
      legend.title = element_text(size = rel(1),  hjust = 1),
      legend.title.align = NULL,
      legend.position = legendPosition,
      legend.direction = NULL,
      legend.justification = "left",
      legend.box = "horizontal",
      
      panel.background = element_rect(linetype=0),
      panel.border = element_rect(fill = NA, color = "grey"),
      panel.grid.major = element_line(colour = "#d8d6d5", size=rel(1.25), linetype = 3),
      panel.grid.minor = element_line(colour = "#edecec", size = rel(1), linetype = 3),
      panel.margin = unit(c(1, 1, 1, 1), "lines"),
      
      strip.background = element_rect(fill = "white", colour = "#838281", linetype=1, size=1),
      strip.text = element_text(size = rel(0.8)),
      strip.text.x = element_text(),
      strip.text.y = element_text(),
      
      plot.background = element_rect(fill = "white", colour="white"),
      plot.title = element_text(face = "bold", size = rel(1.5), hjust=0.5),
      plot.margin = unit(c(6, 20, 9, 50) * 2, "points"),
      complete = TRUE)

  
  thme
}

c48 <- c("#1d915c","#5395b4",
                "#964a48",
                "#2e3b42",
                "#b14e72",
                "#402630","#f1592a",
                "#81aa90","#f79a70", # lt pink
                "#b5ddc2",
                "#8fcc8b", # lt purple
                "#9f1f63", # lt orange
                "#865444", "#a7a9ac",
                "#d0e088","#7c885c","#d22628","#343822","#231f20",
                "#f5ee31","#a99fce","#54525e","#b0accc",
                "#5e5b73","#efcd9f", "#68705d", "#f8f391", "#faf7b6", "#c4be5d", "#764c29", "#c7ac74", "#8fa7aa", "#c8e7dd", "#766a4d", "#e3a291", "#5d777a", "#299c39", "#4055a5", "#b96bac", "#d97646", "#cebb2d", "#bf1e2e", "#d89028", "#85c440", "#36c1ce", "#574a9e")
c3 <- c("#5395b4", "#f1592a", "#85c440")
end.rcode-->


<!--begin.rcode include=FALSE
prefix <- \'$sample\'
titlefile<-"_title.txt"
title <- paste(prefix,titlefile,sep="")
title <- read.delim(title)
title<-title[complete.cases(title),] #Removes cases with null values
barcode <- title\$Sample_ID #create a list of the Sample ID's to replace BC's later 
#Set up to have both barcode and sample_ID
combinedname<-vector()
basequality <- "_ALL_orgbasequalities.txt"
x <- paste(prefix, basequality, sep="")
base2 <- read.delim(x)
readLength <- (nrow(base2)-1)/2
for(i in 1:nrow(title)){
combinedname<-append(combinedname,paste(title\$Barcode[i]," :",barcode[i],sep=""))}
end.rcode-->

<div class="container-fluid">
  <div class="row-fluid">
    <div class="span2">
      <div class="sidebar-nav-fixed">
        <h3>Navigation</h3>
        <table class="table table-hover">
         <tr><td><a href="#titleFile">Title File</a></tr></td>
         <tr><td><a href="#coverage">Plot: Mean Target Coverage</a></tr></td>
         <tr><td><a href="#versions">Software & Script Versions</a></tr></td>
         <tr><td><a href="#readTrimming">Plot: Read trimming rates</a></tr></td>
         <tr><td><a href="#density">Plot: Cluster Density & Alignment Rate</a></tr></td>
         <tr><td><a href="#basequalities">Plot: Base Qualities</a></tr></td>
         <tr><td><a href="#insertsize">Plot: Insert Size Distribution</a></tr></td>
         <tr><td><a href="#hotspots">Table: Hotspot Mutations in Normals</a></tr></td>
         <tr><td><a href="#mixupsHeatmap">Plot: Sample Mix-up Heatmap</a></tr></td>
         <tr><td><a href="#mixupsTables">Table: Sample Mix-up Tables</a></tr></td>
         <tr><td><a href="#majorcontamination">Plot: Major Contamination Check</a></tr></td>
         <tr><td><a href="#minorcontamination">Plot: Minor Contamination Check</a></tr></td>
         <tr><td><a href="#poolNormal">Plot: Pool Normal VF Correlation</a></tr></td>
         <tr><td><a href="#GCcontent">Plot: GC-content</a></tr></td>
         <tr><td><a href="#capturespec">Plot: Capture Specificity</a></tr></td>
         <tr><td><a href="#duplication">Plot: Duplication Rate</a></tr></td>
         <tr><td><a href="#librarysize">Plot: Library Size</a></tr></td>
         
        </table>
      </div><!--sidebar-nav-fixed -->
    </div> <!-- span2 -->
    <div class="span8">
      <h3> Run: <!--rinline I(prefix) --> QC results</h3>
      <p><i>Memorial Sloan Kettering Cancer Center</i><br>
      Diagnostics Molecular Pathology Lab<br>
      <b>IMPACT Test with bait version:</b> <!--rinline I(title\$Bait_version[1])--><br>
      <b>Analysis Date: </b><!--rinline I(as.character(Sys.time()))--><br>
      Metrics Script Version: 1.3</p>
      <p>This sequencing run produced 2X<!--rinline I(as.character(readLength))-->pb reads.</p>
      <br><br>
      <h4><a name="titleFile">Title File</a> </h4>
      <p>Title file shows the samples that were sequenced in this run. You can click on a tumor sample to bring up the variants discovered in the IMPACT panel genes.</p>
      <!--begin.rcode echo=FALSE, results='asis', message=F
      library(xtable)
      hsmetricsfile<-"_ALL_HSmetrics.txt"
      x<- paste(prefix,hsmetricsfile,sep="")
      capture <- read.delim(x)
      title<- title[, c(1, 3, 14, 5, 6, 7, 15, 8:10)]
      title\$MeanCoverage <- capture\$MEAN_TARGET_COVERAGE
      canonical <- "_ALL_Canonical_exoncoverage.txt"
      x<- paste(prefix,canonical,sep="")
      canonical <- read.delim(x)
      a <- colwise(median)(canonical[, c(3:ncol(canonical))])
      canonical <- melt(a)
      canonical\$variable <- paste(capture\$Sample, barcode, sep = " :")
      title\$MedianCoverage <- canonical\$value
      kable(xtable(title),  "html", table.attr="class=\\"table table-striped table-condensed\\", id=\\"titleFile\\"")
      end.rcode-->

      <br><br>
      <h4><a name="coverage">Mean Target Coverage</a></h4>
      <!--begin.rcode echo=FALSE, message=FALSE
      canonical <- "_ALL_Canonical_exoncoverage.txt"
      x<- paste(prefix,canonical,sep="")
      canonical <- read.delim(x)
      a <- colwise(median)(canonical[, c(3:ncol(canonical))])
      canonical <- melt(a)
      canonical\$variable <- paste(capture\$Sample, barcode, sep = " :")
      colnames(canonical) <- c("Samples", "Cov")
      a<-canonical[title\$Class == "Normal", ]["Cov"]
      normal_cov <- mean(a\$Cov)
      a <- canonical[title\$Class == "Tumor", ]["Cov"]
      tumor_cov <- mean(a\$Cov)
      avrg_cov <- mean(canonical\$Cov)
      end.rcode-->
      <p>This chart shows the median canonical exon coverage for each sample. Median coverage across all samples is <b><!--rinline I(round(avrg_cov))-->X.</b> Median coverage across normal samples is <b><!--rinline I(round(normal_cov))-->X</b> and tumor samples is <b><!--rinline I(round(tumor_cov))-->X</b></p>
      <!--begin.rcode coverage, fig.width=30, fig.height=20, out.extra='style="display:block;width:95%"', dev=c('svg', 'pdf'), echo=FALSE
    ggplot(canonical, aes(x = Samples, y = Cov)) + geom_bar(stat = "identity") + my_theme(base_size=24, xAngle=45) + labs(title="Median Canonical Exon Coverage") + xlab("") + ylab("Median X Coverage") + geom_hline(aes(yintercept=200), color="orange", size=2) + geom_hline(aes(yintercept=50), color="red", size=2)
    end.rcode-->
    <br><br>
    <h4><a name="versions">Software & Script Versions</a></h4>
      <p>The following table shows the software and script versions:</p>
      <!--begin.rcode echo=FALSE, results='asis'
    library(xtable)
    PCF<-"RunConfigration.txt"
    base <- dirname(prefix)
    x<-paste(base,PCF,sep="/")
    xdata<-read.delim(PCF)
    kable(xdata, "html", table.attr="class=\\"table table-striped table-condensed\\", id=\\"titleFile\\"")
    end.rcode-->
      <br><br>
      <h4><a name="readTrimming">Read Trimming Rates</a></h4>
      <a href="figure/trimmed_reads.pdf", target="_blank">
      <!--begin.rcode trimmed_reads, fig.width=30, fig.height=20, out.extra='style="display:block;width:95%"', dev=c('svg', 'pdf'), echo=FALSE, message=FALSE, dpi=300
    hsmetricsfile<-"_ALL_HSmetrics.txt"
    x<- paste(prefix,hsmetricsfile,sep="")
    base <- read.delim(x)
    Period<-combinedname #x-axis values
    reads <- data.frame(Samples = Period, Read1 = base\$PerRead1Trimmed, Read2 = base\$PerRead2Trimmed)
    reads.m <- melt(reads)
    reads.m\$value <- reads.m\$value/100
    a <- ggplot(reads.m, aes(x = Samples, y = value, fill = variable))
    a + geom_bar(stat = "identity", position = position_dodge(width=0.7), width = 0.7, color = "black") + my_theme(base_size=24) + scale_fill_manual(name="Type", values = c3) + labs(title="Percentage of Reads Trimmed") + xlab("") + ylab("") + scale_y_continuous(labels=percent)
    end.rcode-->
      </a>
      <br><br>
      <h4><a name="density">Cluster Density</a></h4>
      <!--begin.rcode echo=FALSE, warnings = FALSE, message=FALSE
    density <- data.frame(Samples = Period, BothAlign = base\$both.reads.align, OneAlign = base\$one.read.aligns, NeitherAlign = base\$neither.read.aligns)
    density.m <- melt(density)
    density.m\$value <- density.m\$value/1000000
    end.rcode-->
      <p> This plot shows the numbers of alignment of the read pairs. There were a total of <b><!--rinline I(round(sum(density.m\$value)))--></b> million clusters on the lane. And <b><!--rinline I(round(sum(density.m[density.m\$variable == "BothAlign", ]["value"])/sum(density.m\$value)*100, digits=2))-->%</b> of these had both reads align to the reference genome</p>
      <!--begin.rcode cluster_density, fig.width=30, fig.height=20, out.extra='style="display:block;width:95%"', dev=c('svg', 'pdf'), echo=FALSE,
    ggplot(density.m, aes(x = Samples, y = value, fill = variable)) + geom_bar(stat="identity", width=0.7, color="black")+ my_theme(base_size=24, legendPosition="bottom") + scale_fill_manual(name="Type", values = c3, labels = c("Both Reads Aligned", "One Read Aligned", "Neither Reads Aligned")) + labs(title="Cluster Density & Alignment Rate") + xlab("") + ylab("xMillion")
    density.m2 <- ddply(density.m, .(Samples), mutate, perc = value/sum(value))
    y <- min(density.m2[density.m2\$variable == "BothAlign", ]["perc"])
    y <- y*0.90
    ggplot(density.m2, aes(x = Samples, y = perc, fill = variable)) + geom_bar(stat="identity", position="fill", width=0.7, color="black")+ my_theme(base_size=30, legendPosition="bottom") + scale_fill_manual(name="Type", values = c3, labels = c("Both Reads Aligned", "One Read Aligned", "Neither Reads Aligned")) + labs(title="Cluster Density & Alignment Rate") + xlab("") + ylab("") + scale_y_continuous(label=percent, limits = c(y,1), oob = rescale_none) 
    end.rcode-->
      <br><br>
      <h4><a name="basequalities">Base Quality</a></h4>
      <p> This plot shows the per base quality scores before and after recalibration by <code>GATK</code>.</p>
      <!--begin.rcode basequalities, fig.width=30, fig.height=20, out.extra='style="display:block;width:95%"', dev=c('svg', 'pdf'), echo=FALSE, results='asis'
    basequalityfile<-"_ALL_orgbasequalities.txt"
    x<- paste(prefix,basequalityfile,sep="")
    base2 <- read.delim(x)
    quality <- base2[1:200, ]
    colnames(quality) <- c("cycle", paste(title\$Barcode, title\$Sample_ID, sep = ": "))
    quality.m <- melt(quality, id.vars="cycle")
    quality.m\$type <- "Pre-Recalibration"
    basequalityfile2<-"_ALL_basequalities.txt"
    x<- paste(prefix,basequalityfile2,sep="")
    base3<-read.delim(x)
    colnames(base3) <- c("cycle", paste(title\$Barcode, title\$Sample_ID, sep = ": "))
    quality2 <- base3[1:200, ]
    quality2.m <- melt(quality2, id.vars="cycle")
    quality2.m\$type <- "Post-Recalibration"
    quality3.m <- rbind(quality.m, quality2.m)
    quality3.m\$type <- factor(quality3.m\$type, levels = c("Pre-Recalibration", "Post-Recalibration"))
    preCalib <- gvisLineChart(quality, options = list(width = 1500, fontName='Helvetica', chartArea="{left:50,top:30,width:\\"75%\\",height:\\"75%\\"}"))
    postCalib <- gvisLineChart(quality2, options = list(width = 1500, fontName='Helvetica', chartArea="{left:50,top:30,width:\\"75%\\",height:\\"75%\\"}"))
    print(preCalib, 'chart')
    print(postCalib, 'chart')
    #ggplot(quality3.m, aes(x = cycle, y = value, color = variable)) + geom_line(size=1.2) + my_theme(base_size=26, xAngle=0, legendPosition="right")+ facet_wrap(~type, ncol = 1,) + scale_color_manual(name="Samples", values = rep(c48, 2)) + labs(title="Pre- & Post-Recalibration Quality Scores") + xlab("Quality Score") + ylab("") + guides(colour = guide_legend(override.aes = list(size=5)))
    end.rcode-->
      <br><br>
      <h4><a name="insertsize">Insert Size Distribution</a></h4>
      <p>The following chart shows the distribution of insert sizes for each of the samples.</p>
      <!--begin.rcode insert_size, fig.width=30, fig.height=20, out.extra='style="display:block;width:95%"', dev=c('svg', 'pdf'), echo=FALSE, message=FALSE, results='asis'
    INSERTS<-"_ALL_insertsizemetrics.txt"
    x<- paste(prefix,INSERTS,sep="")
    insert <- read.delim(x)
    insert_label <- paste(paste(title\$Barcode, barcode, sep = ": "), insert[501, 2:ncol(insert)], sep = ", peak=")
    peakHeights <- insert[501, 2:ncol(insert)]
    peakHeights.m <- melt(peakHeights)
    peakHeights_label <- paste(title\$Barcode, barcode, sep = ": ")
    levels(peakHeights.m\$variable) <- peakHeights_label
    insert <- insert[1:500, ]
    insert.m <- melt(insert, id.vars="insert_size")
    insert.m\$insert_size <- as.integer(as.character(insert.m\$insert_size))
    levels(insert.m\$variable) <- insert_label
    insert.m <- insert.m[grep("NTC", insert.m\$variable, invert = T), ]  
    insert_label <- c("insert_size", insert_label)
    colnames(insert) = insert_label
    insert <- gvisLineChart(insert, options = list(width = 1500, height=1000, fontName='Helvetica', chartArea="{left:50,top:30,width:\\"75%\\",height:\\"75%\\"}"))
    print(insert, 'chart')  
    #ggplot(insert.m, aes(x = insert_size, y = value, color = variable)) +  geom_line(size=1.2) + my_theme(base_size=26, xAngle=0, legendPosition="right") + labs(title="Insert Size Distribution") + xlab("Insert size") + ylab("") + scale_color_manual(name="Samples", values = rep(c48, 2)) + guides(colour = guide_legend(override.aes = list(size=5)))
    end.rcode-->
    <!--begin.rcode insert_size_bar_chart, fig.width=30, fig.height=20, out.extra='style="display:block;width:95%"', dev=c('svg', 'pdf'), echo=FALSE, message=FALSE,
    ggplot(peakHeights.m, aes(x = variable, y = value)) + geom_bar(stat="identity") + my_theme(base_size=26, xAngle=45) + xlab("") + ylab("Insert Size") + labs(title="Peak Insert Size Values")
    end.rcode-->
      <br><br>
      <h4><a name="hotspots">Hotspot Mutations in Normals</a></h4>
      <p>The following table shows the presence of hotspot mutations with VF > 1%. </p>
      <!--begin.rcode hotspots, echo=FALSE, results='asis'
    FPR<-"_ALL_genotypehotspotnormals.txt"
    x<- paste(prefix,FPR,sep="")
    hotspots <- read.delim(x, header = F)
    if(length(hotspots\$V1) > 1){  
    colnames(hotspots) <- c("SampleType", "Sample", "Hotspot Cosmic Information", "Chr", "Position", "Ref", "Alt", "DP", "AD", "VF")
    } 
   	kable(xtable(hotspots),  "html", table.attr="class=\\"table table-striped table-condensed\\", id=\\"titleFile\\"")
    
    end.rcode-->
      <br><br>
      <h4><a name="mixupsHeatmap">Sample Mix-ups heatmap</a></h4>
      <p>The Value below the key refers to the fraction of discordant homozygous alleles. Finding a low score between unrelated samples is a red flag. </P>
      <!--begin.rcode samplemixups, fig.width=20, fig.height=20, out.extra='style="display:block;width:75%"', comment=NA, dev=c('svg', 'pdf'), echo=FALSE, error=FALSE
    FPX<-"_ALL_FPCsummary.txt"
    x<- paste(prefix,FPX,sep="")
    xt<-read.table(x, as.is = TRUE, header = TRUE, sep = "\t", row.names = 1)
    colnames(xt) <- combinedname
    rownames(xt) <- combinedname
    pt <- title[title\$Class == "PoolTumor",]["Sample_ID"]
    pn <- title[title\$Class == "PoolNormal",]["Sample_ID"]
    ntc <- title[title\$Class == "NTC",]["Sample_ID"]
    if(length(pn\$Sample_ID) > 0){
    	xt <- xt[, grep(pn\$Sample_ID, colnames(xt), invert=T)]
    	xt <- xt[grep(pn\$Sample_ID, rownames(xt), invert=T), ]
    	Period <- Period[grep(pn\$Sample_ID, Period, invert=T)]
    }
    if(length(ntc\$Sample_ID) > 0){
    	xt <- xt[, grep(ntc\$Sample_ID, colnames(xt), invert=T)]
    	xt <- xt[grep(ntc\$Sample_ID, rownames(xt), invert=T), ]
    	Period <- Period[grep(ntc\$Sample_ID, Period, invert=T)]
    }
    if(length(pt\$Sample_ID) > 0){
    	xt <- xt[, grep(pt\$Sample_ID, colnames(xt), invert=T)]
    	xt <- xt[grep(pt\$Sample_ID, rownames(xt), invert=T), ]
    	Period <- Period[grep(pt\$Sample_ID, Period, invert=T)]
    }
    xm<-as.matrix(xt)
    val = isSymmetric(xm)
    par(mar=c(50, 5, 5, 50),cex.main=2)
    nr <- nrow(xm)
    nc<-ncol(xm)
    if(val == FALSE)
    {
      heatmap.2(xm,Rowv="NA",Colv="NA",dendrogram="none",distfun="dist",hclustfun="hclust",labRow=Period,labCol=Period,cexRow = 1.3 + 1/log10(nr), cexCol = 1.3 + 1/log10(nc),key="TRUE",keysize=0.75,trace="none",density.info=c("none"),margins=c(35,35),col=brewer.pal(5,"Spectral"), main="Sample Mix-Ups", sepcolor="black", colsep=c(1:nrow(xm)), rowsep=c(1:nrow(xm)), sepwidth=c(0.01, 0.01))
    }
    end.rcode-->
    <br><br>
    <h4><a name="mixupsTables">Sample Mix-ups Tables</a></h4>
    <!--begin.rcode samplemixupsTables, fig.width=20, fig.height=20, out.extra='style="display:block;width:75%"', comment=NA, echo=FALSE, error=FALSE
    FPR<-"_ALL_FPCResultsUnMatch.txt"
    x<- paste(prefix,FPR,sep="")
    unmatch<-read.delim(x, header=T, comment.char="#")
    if(nrow(unmatch) > 0){
    	if(length(pn\$Sample_ID) > 0){
    		unmatch <- unmatch[grep(pn\$Sample_ID, unmatch\$Sample1, invert=T), ]
    		unmatch <- unmatch[grep(pn\$Sample_ID, unmatch\$Sample2, invert=T), ]
    	}
    	if(length(ntc\$Sample_ID) > 0){
    		unmatch <- unmatch[grep(ntc\$Sample_ID, unmatch\$Sample1, invert=T), ]
    		unmatch <- unmatch[grep(ntc\$Sample_ID, unmatch\$Sample2, invert=T), ]
   		 }
   	 	if(length(pt\$Sample_ID) > 0){
    		unmatch <- unmatch[grep(pt\$Sample_ID, unmatch\$Sample1, invert=T), ]
    		unmatch <- unmatch[grep(pt\$Sample_ID, unmatch\$Sample2, invert=T), ]
   	 	}
    }
    if(nrow(unmatch) == 0){
    	cat("There are no unexpected matches in this pool.")
    }else{
    	cat("The following table shows the unexpected matches.")
    	kable(unmatch, "html", table.attr="class=\\"table table-striped table-condensed\\", id=\\"titleFile\\"")
    }

    FPR<-"_ALL_FPCResultsUnMismatch.txt"
    x<- paste(prefix,FPR,sep="")
    mismatch<-read.delim(x, header=T, comment.char="#")
    

    if(nrow(mismatch) > 0) {
    	if(length(pn\$Sample_ID) > 0){
    		mismatch <- mismatch[grep(pn\$Sample_ID, mismatch\$Sample1, invert=T), ]
    		mismatch <- mismatch[grep(pn\$Sample_ID, mismatch\$Sample2, invert=T), ]
    	}
    	if(length(ntc\$Sample_ID) > 0){
    		mismatch <- mismatch[grep(ntc\$Sample_ID, mismatch\$Sample1, invert=T), ]
    		mismatch <- mismatch[grep(ntc\$Sample_ID, mismatch\$Sample2, invert=T), ]
   		 }
    	if(length(pt\$Sample_ID) > 0){
    		mismatch <- mismatch[grep(pt\$Sample_ID, mismatch\$Sample2, invert=T), ]
    		mismatch <- mismatch[grep(pt\$Sample_ID, mismatch\$Sample1, invert=T), ]
    	}
    }
    if(nrow(mismatch) == 0){
    	cat("There are no unexpected mismatches in this pool.")
    }else{
    	cat("The following table shows the unexpected mismatches.")
    	kable(mismatch, "html", table.attr="class=\\"table table-striped table-condensed\\", id=\\"titleFile\\"")
    }

    end.rcode-->
      <br><br>
      <h4><a name="majorcontamination">Major Contamination Check</a></h4>
      <!--begin.rcode echo=FALSE
    FPH<-"_ALL_FPhet.txt"
    x<- paste(prefix,FPH,sep="")
    skip <- read.delim(x)
    end.rcode-->
      <p>This chart shows the major contamination check. The mean fraction of positions that are hetrozygous is <b><!--rinline I(round(mean(skip\$PerHetrozygosPos, na.rm=TRUE), digits=5))--></b>.</p>
      <!--begin.rcode major_conatmination, fig.width=30, fig.height=20, out.extra='style="display:block;width:95%"', dev=c('svg', 'pdf'), echo=FALSE, error=FALSE
    maximum<-1.5*(as.numeric(max(skip\$PerHetrozygosPos, na.rm=TRUE)))
    skip\$Sample <- paste(skip\$Sample, barcode, sep=" :")
    if(length(pn\$Sample_ID) > 0){
    		skip <- skip[grep(pn\$Sample_ID, skip\$Sample, invert=T), ]
    }
    if(length(ntc\$Sample_ID) > 0){
    	skip <- skip[grep(ntc\$Sample_ID, skip\$Sample, invert=T), ]
   	}
    if(length(pt\$Sample_ID) > 0){
    	skip <- skip[grep(pt\$Sample_ID, skip\$Sample, invert = T), ]
    }
    ggplot(skip, aes(x = Sample, y = PerHetrozygosPos)) + geom_bar(stat="identity") + ylim(0, 1) + my_theme(base_size=26, xAngle=45) + xlab("") + ylab("Fraction of position that are hetrozygous") + labs(title="Major Contamination Check") + geom_hline(aes(yintercept=0.55), color = "red", size=2)
    end.rcode-->
      <br><br>
      <h4><a name="minorcontamination">Minor Contamination Check</a></h4>
      <!--begin.rcode echo=FALSE
    FPH<-"_ALL_FPavgHom.txt"
    x<- paste(prefix,FPH,sep="")
    skip2 <- read.delim(x)
    end.rcode-->
      <p>This chart shows the minor contamination check. This is based on... Average minor allele frequency is <b><!--rinline I(round(mean(skip2\$AvgMinorHomFreq, na.rm=TRUE), digits=5))--></b>.</p>
      <!--begin.rcode minor_conatmination, fig.width=30, fig.height=20, out.extra='style="display:block;width:95%"', dev=c('svg', 'pdf'), echo=FALSE, errro=FALSE
    maximum<-1.5*(as.numeric(max(skip2\$AvgMinorHomFreq, na.rm=TRUE)))
    skip2\$Sample <- paste(skip2\$Sample, barcode, sep=" :")
    if(length(pn\$Sample_ID) > 0){
    	skip2 <- skip2[grep(pn\$Sample_ID, skip2\$Sample, invert=T), ]
    }
    if(length(ntc\$Sample_ID) > 0){
    	skip2 <- skip2[grep(ntc\$Sample_ID, skip2\$Sample, invert=T), ]
   	}
    skip2 <- skip2[grep(pt\$Sample_ID, skip2\$Sample, invert = T), ]
    ggplot(skip2, aes(x = Sample, y = AvgMinorHomFreq)) + geom_bar(stat="identity") + ylim(0, maximum) + my_theme(base_size=26, xAngle=45) + xlab("") + ylab("Avg. Minor Allele Frequency at Homozygous Position") + labs(title="Minor Contamination Check") + geom_hline(aes(yintercept=0.02), color = "red", size=2) + geom_hline(aes(yintercept=0.01), color = "orange", size=2) + scale_color_hue(name="Samples")
    end.rcode-->
     <br><br>
     
     <h4><a name="FPsummary">Distribution of variant frequencies for samples with average MAF > 1% in minor contamination</a></h4>
 	<!--begin.rcode minor_conatmination_fpsummary, fig.width=30, fig.height=20, out.extra='style="display:block;width:95%"', dev=c('svg', 'pdf'), echo=FALSE,
	FPH<-"_ALL_FPavgHom.txt"
    x<- paste(prefix,FPH,sep="")
    skip2 <- read.delim(x)
	
    poolNormal <- title[title\$Class == "PoolNormal", ]["Barcode"]
   # skip2 <- skip2[grep(poolNormal\$Barcode, skip\$Sample, invert = T), ]
    
    poolTumor <- title[title\$Class == "PoolTumor", ]["Barcode"]
   # skip2 <- skip2[grep(poolTumor\$Barcode, skip\$Sample, invert = T), ]
    
    ntc <- title[title\$Class == "NTC", ]["Barcode"]
    skip2 <- skip2[grep(ntc\$Barcode, skip2\$Sample, invert = T), ]
   
	skip3 <- skip2[skip2\$AvgMinorHomFreq > 0.01, ]
	skip3 <- na.omit(skip3)
	
	if(nrow(skip3) > 0){
		FPS <- "_ALL_FPsummary.txt"
		x <- paste(prefix, FPS, sep="")
		fpsummary <- read.delim(x)
		fpsum <- data.frame()
		for (i in 1:nrow(skip3)){
		  a <- as.character(skip3[i, 1])
		  colWanted <- paste(a, "MinorAlleleFreq", sep = "_")
		  b <- fpsummary[colWanted]
		  c <- title[title\$Barcode == a, ]["Sample_ID"]
		  e <- title[title\$Barcode == a, ]["Class"]
		  d <- data.frame(Source = paste(a, c\$Sample_ID, e\$Class, sep =" "), FPsum = b[, 1])
		  fpsum <- rbind(fpsum, d)
		}
		if(i <= 4){
		  ncol = 2
		}else if(i > 4 && i <= 9){
		  ncol = 3
		}else if(i > 9){
		  ncol = 4
		}
		ggplot(fpsum, aes(x = FPsum)) + geom_histogram(binwidth=0.001, fill = "black", color = "white") + my_theme(xAngle=0, hJust=0.5, base_size=24) + scale_x_continuous(labels = percent, breaks=c(0, 0.02, 0.04, 0.06, 0.08, 0.1), limits = c(0,0.1)) + theme(strip.text.x = element_text(size = 30))+ facet_wrap(~Source, ncol = ncol, scales="free") + xlab("Minor Allele Frequency")
	}else{
		cat("There are no samples with MAF > 0.01 in this pool.")
	}
	
     end.rcode-->
     <br><br>
     
     
      <h4><a name="PoolNormal">Pool Normal Genotyping</a></h4>
      <p>This plot shows the correlation of pre-determined SNP variant frequencies in the pool of normal sample. </p>
      <!--begin.rcode poolNormal, warning=F, echo=FALSE , fig.width=25, fig.height=25, out.extra='style="display:block;width:70%"', dev=c('svg', 'pdf'), 
      a <- "_ALL_GenotypePooledNormal.txt"
      poolNormal <- read.delim(paste(prefix, a, sep = ""))
      corTest <- with(poolNormal, cor.test(ExpectedVF, ObservedVF))
      ggplot(poolNormal, aes(x = ExpectedVF, y = ObservedVF)) + geom_point(size=6, alpha=0.7) + my_theme(xAngle = 0, hJust = 0.5, base_size = 26) + labs(title=paste("Expected and Observed variant frequencies\n for SNPs present in pool normal. Cor:", round(corTest\$estimate, digits=3), sep = " "))
      end.rcode-->
      <br><br>
      <h4><a name="GCcontent">GC Content</a></h4>
      <p>This plot shows the normalized coverage of each region plotted against the GC content of said region</p>
      <!--begin.rcode gc_content, fig.width=30, fig.height=20, out.extra='style="display:block;width:95%"', dev=c('svg', 'pdf'), echo=FALSE, results='asis'
    GCP<-"_ALL_gcbias.txt"
    x<- paste(prefix,GCP,sep="")
    xt<-read.delim(x)
    yvar <- paste(title\$Barcode, title\$Sample_ID, sep = ": ")
    colnames(xt) <- c("GC_Content", yvar)
    GC <- gvisLineChart(xt, xvar ="GC_Content", yvar = yvar, options = list(width = 1500, height = 1000, fontName='Helvetica', chartArea="{left:50,top:30,width:\\"75%\\",height:\\"75%\\"}"))
    print(GC, 'chart')
    #xt.m <- melt(xt, id.vars="X")
    #xt.m <- xt.m[grep("NTC", xt.m\$variable, invert = T), ]
    #ggplot(xt.m, aes(x = X, y = value, color = variable)) + geom_line(size=1.2) +my_theme(base_size=26, xAngle=0, legendPosition="right") + labs(title="Normalized Coverage vs GC-Content") + xlab("GC-Content") + ylab("Normalized Coverage") + scale_color_manual(name="Samples", values = rep(c48, 2)) + guides(colour = guide_legend(override.aes = list(size=5)))
    end.rcode-->
      <br><br>
      <h4><a name="capturespec">Capture Specificity</a></h4>
      <!--begin.rcode echo=FALSE
    hsmetricsfile<-"_ALL_HSmetrics.txt"
    x<- paste(prefix,hsmetricsfile,sep="")
    capture <- read.delim(x)
    baits <- data.frame(Sample =  paste(capture\$Sample, barcode, sep = " :"), OnBait = capture\$ON_BAIT_BASES, NearBait = capture\$NEAR_BAIT_BASES, OffBait = capture\$OFF_BAIT_BASES)
    baits.m <- melt(baits, id.vars="Sample")
    on_total <- sum(as.numeric(baits\$OnBait))
    off_total <- sum(as.numeric(baits\$OffBait))
    near_total <- sum(as.numeric(baits\$NearBait))
    total <- on_total+near_total+off_total
    on_near_perc <- (on_total+near_total)/total*100
    on_bait_perc <- on_total/total*100
    on_target_total <- sum(as.numeric(capture\$ON_TARGET_BASES))
    on_target_perc <- on_target_total/total*100
    end.rcode-->
      <p> This plot shows the capture specificty of the probes.Some of the details are as follows: <br>
      - Average % selected on/near bait = <b><!--rinline I(round(on_near_perc, digits=1))-->%</b><br>
      - Average % on bait = <b><!--rinline I(round(on_bait_perc, digits = 1))-->%</b><br>
      - Average % on target = <b><!--rinline I(round(on_target_perc, digits = 1))-->%</b></p>
      <!--begin.rcode capture_spec, fig.width=30, fig.height=20, out.extra='style="display:block;width:95%"', dev=c('svg', 'pdf'), echo=FALSE
    ggplot(baits.m, aes(x = Sample, y = value/1000000, fill = variable)) + geom_bar(stat="identity", width=0.7, color="black")+ my_theme(base_size=24, xAngle=45) + scale_fill_manual(name="Type", values = c3, labels = c("On Bait Bases", "Near Bait Bases", "Off Bait Bases")) + labs(title="Capture Specificity") + xlab("") + ylab("Total Bases (millions)")
    ggplot(baits.m, aes(x = Sample, y = value/1000000, fill = variable)) + geom_bar(stat="identity", position="fill",width=0.7, color="black")+ my_theme(base_size=24, xAngle=45) + scale_fill_manual(name="Type", values = c3, labels = c("On Bait Bases", "Near Bait Bases", "Off Bait Bases")) + labs(title="Capture Specificity") + xlab("") + ylab("Percent bases") + scale_y_continuous(labels=percent)
    end.rcode-->
      <br><br>
      <h4><a name="duplication">Duplication Rate</a></h4>
      <!--begin.rcode echo=FALSE
    duplication <- data.frame(Samples = paste(capture\$Sample, barcode, sep = " :"), DupRate = capture\$PERCENT_DUPLICATION)
    end.rcode-->
      <p>This chart shows the duplication rates.Average duplication rate is <b> <!--rinline I(round(mean(duplication\$DupRate)*100))-->%</b></p>
      <!--begin.rcode duplication_rate, fig.width=30, fig.height=20, out.extra='style="display:block;width:95%"', dev=c('svg', 'pdf'), echo=FALSE
    ggplot(duplication, aes(x = Samples, y = DupRate)) + geom_bar(stat = "identity") + my_theme(base_size=24, xAngle=45) + labs(title="Estimated Duplication Rate") + xlab("") + ylab("Duplication Rate")+ scale_y_continuous(labels=percent)
    end.rcode-->
      <br><br>
      <h4><a name="librarysize">Estimated Library Size</a></h4>
      <!--begin.rcode echo=FALSE
    complex <- data.frame(Samples = paste(capture\$Sample, barcode, sep = " :"), Comp = capture\$ESTIMATED_LIBRARY_SIZE)
    end.rcode-->
      <p> The following chart shows the estimated total library size. Average library size is <b><!--rinline I(round(mean(complex\$Comp)/1000000))--></b> million and total library size is <b><!--rinline I(round(sum(complex\$Comp)/1000000))--> million</b></p>
      <!--begin.rcode complexity, fig.width=30, fig.height=20, out.extra='style="display:block;width:95%"', dev=c('svg', 'pdf'), echo=FALSE
    ggplot(complex, aes(x = Samples, y = Comp/1000000)) + geom_bar(stat = "identity") + my_theme(base_size=24, xAngle=45) + labs(title="Estimated Library Size") + xlab("") + ylab("xMillion")
    end.rcode-->
      <br><br>
      
  

    </div><!-- span8-->
  </div><!--row-fluid-->
</div><!-- container-->
</body>
</html>
ENDSCRIPT

close RSCRIPT;

my $output2 = $sample."_All_Metrics.Rnw";
open RSCRIPT, ">", $output2 || $!;
print RSCRIPT <<ENDSCRIPT; 
\\documentclass[landscape]{article}
\\usepackage[sc]{mathpazo}
\\usepackage[T1]{fontenc}
\\renewcommand{\\familydefault}{\\sfdefault}
\\usepackage{helvet}
\\usepackage{geometry}
\\geometry{verbose,tmargin=0.4cm,bmargin=0.4cm,lmargin=0.4cm,rmargin=0.4cm}
\\setcounter{secnumdepth}{2}
\\setcounter{tocdepth}{2}
\\usepackage{url}
\\usepackage{graphicx}
\\usepackage{caption}
\\usepackage{subcaption}
\\usepackage{longtable}
\\usepackage[unicode=true,pdfusetitle, bookmarks=true,bookmarksnumbered=true,bookmarksopen=true,bookmarksopenlevel=2,breaklinks=false,pdfborder={0 0 1},backref=false,colorlinks=false]{hyperref}
\\hypersetup{pdfstartview={XYZ null null 1}}
\\usepackage{breakurl}
\\begin{document}


<<setup, include=FALSE, echo=FALSE, cache = FALSE>>=
#opts_chunk\$set(fig.width=4, fig.height=3, par=TRUE, out.width='2in', fig.pos='H', fig_path="figure/")
my_theme <- function(base_size = 16, base_family="Helvetica", dkpanel=FALSE, legendPosition = "right", xAngle = 45, hJust = 1) {
  thme <- theme(
      line = element_line(colour = "#838281", size = 0.5,
                          linetype = 1, lineend = "round"),
      rect = element_rect(fill = "white",
                          colour = NA, size = 0.5, linetype = 1),
      text = element_text(family = base_family, face = "plain",
                          colour = "black", size = base_size, hjust = 0.5, vjust = 0.5,
                          angle = 0, lineheight = 1),
      axis.text = element_text(size = rel(1), face = "italic"),
      axis.text.x = element_text(vjust = 1, hjust=hJust, angle = xAngle),
      axis.text.y = element_text(hjust = 1),
      axis.line = element_line(colour = "#d8d6d5"),
      axis.line.y = element_line(),
      axis.line.x = element_line(),
      axis.ticks = element_line(),
      axis.title = element_text(size = rel(1)),
      axis.title.x = element_text(),
      axis.title.y = element_text(angle = 90),
 
      axis.ticks.length = unit(base_size * 0.2 , "points"),
      axis.ticks.margin = unit(base_size * 0.5, "points"),
      
      legend.background = element_rect(color = "grey", linetype=1, size = rel(0.5)),
      legend.margin = unit(base_size * 2.5, "points"),
      legend.key = element_rect(),
      legend.key.size = unit(2, "lines"),
      legend.key.height = NULL,
      legend.key.width = NULL,
      legend.text = element_text(size = rel(1.1)),
      legend.text.align = NULL,
      legend.title = element_text(size = rel(1),  hjust = 1),
      legend.title.align = NULL,
      legend.position = legendPosition,
      legend.direction = NULL,
      legend.justification = "left",
      legend.box = "horizontal",
      
      panel.background = element_rect(linetype=0),
      panel.border = element_rect(fill = NA, color = "grey"),
      panel.grid.major = element_line(colour = "#d8d6d5", size=rel(1.25), linetype = 3),
      panel.grid.minor = element_line(colour = "#edecec", size = rel(1), linetype = 3),
      panel.margin = unit(c(1, 1, 1, 1), "lines"),
      
      strip.background = element_rect(fill = "white", colour = "#838281", linetype=1, size=1),
      strip.text = element_text(size = rel(0.8)),
      strip.text.x = element_text(),
      strip.text.y = element_text(),
      
      plot.background = element_rect(fill = "white", colour="white"),
      plot.title = element_text(face = "bold", size = rel(1.5), hjust=0.5),
      plot.margin = unit(c(6, 20, 9, 50) * 2, "points"),
      complete = TRUE)

  
  thme
}
c48 <- c("#1d915c","#5395b4",
                "#964a48",
                "#2e3b42",
                "#b14e72",
                "#402630","#f1592a",
                "#81aa90","#f79a70", # lt pink
                "#b5ddc2",
                "#8fcc8b", # lt purple
                "#9f1f63", # lt orange
                "#865444", "#a7a9ac",
                "#d0e088","#7c885c","#d22628","#343822","#231f20",
                "#f5ee31","#a99fce","#54525e","#b0accc",
                "#5e5b73","#efcd9f", "#68705d", "#f8f391", "#faf7b6", "#c4be5d", "#764c29", "#c7ac74", "#8fa7aa", "#c8e7dd", "#766a4d", "#e3a291", "#5d777a", "#299c39", "#4055a5", "#b96bac", "#d97646", "#cebb2d", "#bf1e2e", "#d89028", "#85c440", "#36c1ce", "#574a9e")
c3 <- c("#5395b4", "#f1592a", "#85c440")
prefix <- \"$sample\"
titlefile<-"_title.txt"
title <- paste(prefix,titlefile,sep="")
title <- read.delim(title)
title<-title[complete.cases(title),] #Removes cases with null values
barcode <- title\$Sample_ID #create a list of the Sample ID's to replace BC's later 
#Set up to have both barcode and sample_ID
combinedname<-vector()
for(i in 1:nrow(title)){
combinedname<-append(combinedname,paste(title\$Barcode[i]," :",barcode[i],sep=""))}
basequalityfile<-"_ALL_orgbasequalities.txt"
x<- paste(prefix,basequalityfile,sep="")
base2 <- read.delim(x)
readLength <- (nrow(base2)-1)/2
@
\\title{Run: \\Sexpr{prefix}}
\\author{Memorial Sloan Kettering Cancer Center\\\\
        Diagnostics Molecular Pathology Lab,\\\\
        IMPACT Test with bait version: \\Sexpr{title\$Bait_version[1]},\\\\
        Analysis Date: \\Sexpr{as.character(Sys.time())},
        Metrics Script Version: 1.3}
\\maketitle
This document shows the QC metrics results for run: \\Sexpr{prefix}. This sequencing run was performed on a HiSeq 2500 with \\textbf{\\Sexpr{readLength}}pb read pairs. We can add more information here if we want to. 
\\listoffigures
\\listoftables
\\clearpage
\\section{Title File}
<<titleFile, echo=FALSE, results="asis", message=FALSE>>=
library(xtable)
hsmetricsfile<-"_ALL_HSmetrics.txt"
x<- paste(prefix,hsmetricsfile,sep="")
capture <- read.delim(x)
canonical <- "_ALL_Canonical_exoncoverage.txt"
x<- paste(prefix,canonical,sep="")
canonical <- read.delim(x)
a <- colwise(median)(canonical[, c(3:ncol(canonical))])
canonical <- melt(a)
title\$MedianCoverage <- canonical[, 2]
title <- title[, c(1, 3:12)]
titleTable <- xtable(title, label="title_file", caption="Title File")
print.xtable(titleTable, hline.after=c(-1, 0, 1:nrow(title)), scalebox=0.8,include.rownames=FALSE)
@

\\clearpage
\\section{Software Versions}
<<versions, echo=FALSE, results="asis">>=
library(xtable)
PCF<-"RunConfigration.txt"
base <- dirname(prefix)
x<-paste(base,PCF,sep="/")
xdata<-read.delim(PCF)
xtable(xdata, label="versions", caption="Software Versions",include.rownames=FALSE)
@

\\clearpage

\\begin{figure}[hl]
  \\centering
    \\includegraphics[width = 0.8\\textwidth]{figure/trimmed_reads.pdf}
    \\caption{Percent read trimming}
  \\label{read_trimming}
\\end{figure}


\\clearpage
<<density, echo=FALSE, warning=FALSE, message=FALSE>>=
hsmetricsfile<-"_ALL_HSmetrics.txt"
x<- paste(prefix,hsmetricsfile,sep="")
base <- read.delim(x)
Period<-combinedname #x-axis values
a <- data.frame(Samples = Period, Read1 = base\$PerRead1Trimmed, Read2 = base\$PerRead2Trimmed)
a.m <- melt(a)
a.m\$value <- a.m\$value/100
density <- data.frame(Samples = Period, BothAlign = base\$both.reads.align, OneAlign = base\$one.read.aligns, NeitherAlign = base\$neither.read.aligns)
density.m <- melt(density)
density.m\$value <- density.m\$value/1000000
@

\\begin{figure}[hl]
  \\centering
    \\includegraphics[width = 0.9\\textwidth]{figure/cluster_density1.pdf}
    \\caption[Cluster Density]{Cluster Density: This plot shows the numbers of alignment of the read pairs. There were a total of \\textbf{\\Sexpr{round(sum(density.m\$value))}} million clusters on the lane. And \\textbf{\\Sexpr{round(sum(density.m[density.m\$variable == "BothAlign", ]["value"])/sum(density.m\$value)*100, digits=2)}}\\% of these had both reads align to the reference genome}
  \\label{cluster_density1}
\\end{figure}
\\begin{figure}[hl]
  \\centering
    \\includegraphics[width = 0.9\\textwidth]{figure/cluster_density2.pdf}
    \\caption[Cluster Density]{Cluster Density: This plot shows the numbers of alignment of the read pairs. There were a total of \\textbf{\\Sexpr{round(sum(density.m\$value))}} million clusters on the lane. And \\textbf{\\Sexpr{round(sum(density.m[density.m\$variable == "BothAlign", ]["value"])/sum(density.m\$value)*100, digits=2)}}\\% of these had both reads align to the reference genome}
  \\label{cluster_density2}
\\end{figure}


\\clearpage

\\begin{figure}[hl]
  \\centering
    \\includegraphics[width = 0.9\\textwidth]{figure/basequalities.pdf}
    \\caption{Base Qualities}
  \\label{base_qualities}
\\end{figure}


\\clearpage

\\begin{figure}[hl]
  \\centering
    \\includegraphics[width = 0.9\\textwidth]{figure/insert_size1.pdf}
    \\caption{Insert Size Distribution}
  \\label{insert_sizes}
\\end{figure}

\\clearpage

\\begin{figure}[hl]
  \\centering
    \\includegraphics[width = 0.9\\textwidth]{figure/insert_size2.pdf}
    \\caption{Peak Insert Size Values}
  \\label{insert_sizes}
\\end{figure}

\\clearpage
The following table shows the presence of hotspot mutations with VF > 1\\%.
<<hotspots, echo=FALSE, results="asis">>=
FPR<-"_ALL_genotypehotspotnormals.txt"
x<- paste(prefix,FPR,sep="")
hotspots <- read.delim(x, header = F)
if(length(hotspots\$V1) > 1){
colnames(hotspots) <- c("SampleType", "Sample", "Hotspot Cosmic Information", "Chr", "Position", "Ref", "Alt", "DP", "AD", "VF")
hotspotTable <- xtable(hotspots, caption=" Possible Hotspot Mutations in Normals", label="hotspots")
print.xtable(hotspotTable, scalebox=0.8,include.rownames=FALSE)
}else{
	cat("There are no hotspot mutations detected in the normals sequenced in this pool.")
}
@

\\clearpage

\\begin{figure}[hl]
  \\centering
    \\includegraphics[width = 0.65\\textwidth]{figure/samplemixups.pdf}
    \\caption{Sample Mixups: The Value below the key refers to the fraction of discordant homozygous alleles. Finding a low score between unrelated samples is a red flag}
  \\label{sample_mixups}
\\end{figure}

\\clearpage
<<mixups, echo=F, results="asis">>=
FPR1<-"_ALL_FPCResultsUnMatch.txt"
x1<- paste(prefix,FPR1,sep="")
unmatch<-read.delim(x1, header=T, comment.char="#")
pt <- title[title\$Class == "PoolTumor",]["Sample_ID"]
pn <- title[title\$Class == "PoolNormal",]["Sample_ID"]
ntc <- title[title\$Class == "NTC",]["Sample_ID"]
if(nrow(unmatch) == 0){
    cat("There are no unexpected matches in this pool.")
}else{
	if(length(pn\$Sample_ID) > 0){
    	unmatch <- unmatch[grep(pn\$Sample_ID, unmatch\$Sample1, invert=T), ]
    	unmatch <- unmatch[grep(pn\$Sample_ID, unmatch\$Sample2, invert=T), ]
    }
   	if(length(ntc\$Sample_ID) > 0){
   		unmatch <- unmatch[grep(ntc\$Sample_ID, unmatch\$Sample1, invert=T), ]
   		unmatch <- unmatch[grep(ntc\$Sample_ID, unmatch\$Sample2, invert=T), ]
   	}
    if(length(pt\$Sample_ID) > 0){
    		unmatch <- unmatch[grep(pt\$Sample_ID, unmatch\$Sample1, invert=T), ]
    		unmatch <- unmatch[grep(pt\$Sample_ID, unmatch\$Sample2, invert=T), ]
   	}
	mixupTable1 <- xtable(unmatch, caption="Unexpected matches.", label="mixups1")
	print.xtable(mixupTable1,tabular.environment='longtable',floating=FALSE,include.rownames=FALSE)
    }
@
\\clearpage
<<mixups2, echo=F, results="asis">>=
FPR2<-"_ALL_FPCResultsUnMismatch.txt"
x2<- paste(prefix,FPR2,sep="")
mismatch<-read.delim(x2, header=T, comment.char="#")
if(nrow(mismatch) == 0){
    cat("There are no unexpected mismatches in this pool.")
}else{
    if(length(pn\$Sample_ID) > 0){
    		mismatch <- mismatch[grep(pn\$Sample_ID, mismatch\$Sample1, invert=T), ]
    		mismatch <- mismatch[grep(pn\$Sample_ID, mismatch\$Sample2, invert=T), ]
    	}
    	if(length(ntc\$Sample_ID) > 0){
    		mismatch <- mismatch[grep(ntc\$Sample_ID, mismatch\$Sample1, invert=T), ]
    		mismatch <- mismatch[grep(ntc\$Sample_ID, mismatch\$Sample2, invert=T), ]
   		 }
    	if(length(pt\$Sample_ID) > 0){
    		mismatch <- mismatch[grep(pt\$Sample_ID, mismatch\$Sample2, invert=T), ]
    		mismatch <- mismatch[grep(pt\$Sample_ID, mismatch\$Sample1, invert=T), ]
    	}
    mixupTable2 <- xtable(mismatch, caption="Unexpected mis-matches.", label="mixups2")
	print.xtable(mixupTable2,,tabular.environment='longtable',floating=FALSE,include.rownames=FALSE)
}

@
\\clearpage
<<major, echo=FALSE>>=
FPH<-"_ALL_FPhet.txt"
      x<- paste(prefix,FPH,sep="")
      skip <- read.delim(x)
@
\\begin{figure}[hl]
  \\centering
    \\includegraphics[width = 0.9\\textwidth]{figure/major_conatmination.pdf}
    \\caption{Major Contamination Check: The mean fraction of positions that are hetrozygous is \\textbf{\\Sexpr{round(mean(skip\$PerHetrozygosPos, na.rm=TRUE), digits=5)}}}
  \\label{major_cont}
\\end{figure}

\\clearpage
<<minor, echo=FALSE>>=
FPH<-"_ALL_FPavgHom.txt"
      x<- paste(prefix,FPH,sep="")
      skip2 <- read.delim(x)
@
\\begin{figure}[hl]
  \\centering
    \\includegraphics[width = 0.9\\textwidth]{figure/minor_conatmination.pdf}
    \\caption{Minor Contamination Check: Average minor allele frequency is \\textbf{\\Sexpr{round(mean(skip2\$AvgMinorHomFreq, na.rm=TRUE), digits=5)}}}
  \\label{minor_cont}
\\end{figure}

ENDSCRIPT

my $file = glob("*FPavgHom.txt");
chomp($file);
open IN, "<", $file || die $!;
my $i = 0;
while(<IN>){
	chomp;
	next if(/Sample/);
	my @line = split("\t");
	#print "$line[1]\n";
	if($line[1] > 0.01){
		$i++;
	}
}
close IN;
#print "$i\n";
if($i > 0){
	print RSCRIPT "\\clearpage\n";
	print RSCRIPT "\\begin{figure}[hl]\n";
  	print RSCRIPT "\\centering\n";
    print RSCRIPT "\\includegraphics[width = 0.9\\textwidth]{figure/minor_conatmination_fpsummary.pdf}\n";
    print RSCRIPT "\\caption{Samples with average minor allele frequency > 1\\%.}\n";
  	print RSCRIPT "\\label{minor_allele_freq}\n";
	print RSCRIPT "\\end{figure}\n";
}

my $GenotypePoolNormalFile = $sample. "_ALL_GenotypePooledNormal.txt";
if(-e $GenotypePoolNormalFile){
	print RSCRIPT "\\clearpage\n";
	print RSCRIPT "<<poolNormal, echo=FALSE, warning=FALSE>>=\n";
	print RSCRIPT "a <- \"_ALL_GenotypePooledNormal.txt\"\n";
	print RSCRIPT "poolNormal <- read.delim(paste(prefix, a, sep = \"\"))\n";
	print RSCRIPT "corTest <- with(poolNormal, cor.test(ExpectedVF, ObservedVF))\n";
	print RSCRIPT "@\n";
	print RSCRIPT "\n";
	print RSCRIPT "\\begin{figure}[hl]\n";
	print RSCRIPT "\\centering\n";
	print RSCRIPT "   \\includegraphics[width = 0.7\\textwidth]{figure/poolNormal.pdf}\n";
	print RSCRIPT "   \\caption{Pool Normal: The observed variant frequencies of pre-determined SNPs in a pool of 10 normal samples vs the expected variant frequencies. Correlation: \\textbf{\\Sexpr{round(corTest\$estimate, digits = 3)}}.}\n";
	print RSCRIPT " \\label{GC content}\n";
	print RSCRIPT "\\end{figure}\n";

}

print RSCRIPT <<ENDSCRIPT; 
\\clearpage
\\begin{figure}[hl]
  \\centering
    \\includegraphics[width = 0.9\\textwidth]{figure/gc_content.pdf}
    \\caption{GC Content}
  \\label{GC content}
\\end{figure}


\\clearpage
<<capture_spec, echo=FALSE>>=
hsmetricsfile<-"_ALL_HSmetrics.txt"
x<- paste(prefix,hsmetricsfile,sep="")
capture <- read.delim(x)
baits <- data.frame(Sample = capture\$Sample, OnBait = capture\$ON_BAIT_BASES, NearBait = capture\$NEAR_BAIT_BASES, OffBait = capture\$OFF_BAIT_BASES)
baits.m <- melt(baits, id.vars="Sample")
on_total <- sum(as.numeric(baits\$OnBait))
off_total <- sum(as.numeric(baits\$OffBait))
near_total <- sum(as.numeric(baits\$NearBait))
total <- on_total+near_total+off_total
on_near_perc <- (on_total+near_total)/total*100
on_bait_perc <- on_total/total*100
on_target_total <- sum(as.numeric(capture\$ON_TARGET_BASES))
on_target_perc <- on_target_total/total*100
@

\\begin{figure}[hl]
  \\centering
    \\includegraphics[width = 0.9\\textwidth]{figure/capture_spec1.pdf}
    \\caption{Capture Specificity: Average \\% selected on/near bait = \\textbf{\\Sexpr{round(on_near_perc, digits=1)}}\\%. Average \\% on bait = \\textbf{\\Sexpr{round(on_bait_perc, digits = 1)}}\\%. Average \\% on target = \\textbf{\\Sexpr{round(on_target_perc, digits = 1)}}\\%}
  \\label{capture_spec3}
\\end{figure}

\\begin{figure}[hl]
  \\centering
    \\includegraphics[width = 0.9\\textwidth]{figure/capture_spec2.pdf}
    \\caption[Capture Specificity]{Capture Specificity: Average \\% selected on/near bait = \\textbf{\\Sexpr{round(on_near_perc, digits=1)}}\\%. Average \\% on bait = \\textbf{\\Sexpr{round(on_bait_perc, digits = 1)}}\\%. Average \\% on target = \\textbf{\\Sexpr{round(on_target_perc, digits = 1)}}\\%}
  \\label{capture_spec3}
\\end{figure}

\\clearpage
<<dupRate, echo=FALSE>>=
duplication <- data.frame(Samples = paste(capture\$Sample, barcode, sep = " :"), DupRate = capture\$PERCENT_DUPLICATION)
@
\\begin{figure}[hl]
  \\centering
    \\includegraphics[width = 0.9\\textwidth]{figure/duplication_rate.pdf}
    \\caption{Duplication Rate: Average duplication rate is \\textbf{\\Sexpr{round(mean(duplication\$DupRate)*100)}}\\%.}
  \\label{dup_rate}
\\end{figure}

\\clearpage
<<complexity, echo=FALSE>>=
complex <- data.frame(Samples = paste(capture\$Sample, barcode, sep = " :"), Comp = capture\$ESTIMATED_LIBRARY_SIZE)
@
\\begin{figure}[hl]
  \\centering
    \\includegraphics[width = 0.9\\textwidth]{figure/complexity.pdf}
    \\caption{Estimated Library Size: Average library size is \\textbf{\\Sexpr{round(mean(complex\$Comp)/1000000)}} million and total library size is \\textbf{\\Sexpr{round(sum(complex\$Comp)/1000000)}} million}
  \\label{library_size}
\\end{figure}

\\clearpage
<<coverage, echo=FALSE, message=FALSE>>=
canonical <- "_ALL_Canonical_exoncoverage.txt"
x<- paste(prefix,canonical,sep="")
canonical <- read.delim(x)
a <- colwise(median)(canonical[, c(3:ncol(canonical))])
canonical <- melt(a)
canonical\$variable <- paste(capture\$Sample, barcode, sep = " :")
colnames(canonical) <- c("Samples", "Cov")
a<-canonical[title\$Class == "Normal", ]["Cov"]
a<-canonical[title\$Class == "Normal", ]["Cov"]
normal_cov <- median(a\$Cov)
a <- canonical[title\$Class == "Tumor", ]["Cov"]
tumor_cov <- median(a\$Cov)
avrg_cov <- median(canonical\$Cov)
@
\\begin{figure}[hl]
  \\centering
    \\includegraphics[width = 0.9\\textwidth]{figure/coverage.pdf}
    \\caption{Median Target Coverage: Median canonical exon coverage across all samples is \\textbf{\\Sexpr{round(avrg_cov, digits = 0)}}X. Median coverage across normal samples is \\textbf{\\Sexpr{round(normal_cov, digits = 0)}}X and tumor samples is \\textbf{\\Sexpr{round(tumor_cov, digits = 0)}}X}
  \\label{library_size}
\\end{figure}




%
\\end{document}
ENDSCRIPT

open(OUT, ">", "metrics.r") or die $!;
print OUT "library(knitr)\n";
print OUT "library(ggplot2)\n";
print OUT "library(plyr)\n"; 
print OUT "library(grid)\n"; 
print OUT "library(reshape)\n";
print OUT "library(scales)\n";
print OUT "library(xtable)\n";
print OUT "knit(\"$output\")\n";
print OUT "knit2pdf(\"$output2\")\n";
close OUT;

`$R CMD BATCH metrics.r`;

my $toBeRemoved = $sample."_All_Metrics.tex";
`rm $toBeRemoved`;
$toBeRemoved =~ s/tex/aux/;
`rm $toBeRemoved`;
$toBeRemoved =~ s/aux/out/;
`rm $toBeRemoved`;
$toBeRemoved =~ s/out/lof/;
`rm $toBeRemoved`;
$toBeRemoved =~ s/lof/lot/;
`rm $toBeRemoved`;
$toBeRemoved =~ s/lot/log/;
`rm $toBeRemoved`;
