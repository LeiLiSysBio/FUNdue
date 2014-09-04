#!/usr/bin/env Rscript

#########################
# installing PIANO
# source("http://bioconductor.org/biocLite.R")
# biocLite("piano")
suppressPackageStartupMessages(library("piano"))
suppressPackageStartupMessages(library("optparse"))
myParseArgs <- function(usage,option_list,filenames.exist=NULL,multiple.options=NULL,mandatory=NULL) {
  parser <- OptionParser(usage = usage, option_list=option_list)
  opt <- parse_args(parser)
}
perror <- function(...) {
  cat(paste("[ERROR] ",...,"\n",sep=""),file=stderr())
}
pinfo <- function(...) {
  cat(paste("[INFO] ",...,"\n",sep=""))
}
pdebug.save.state <- function(file,stage="undef") {  
  if (pdebug.enabled) {
	assign("pdebug.stage",stage,envir = .GlobalEnv)
	save.image(file=paste(file,".Rdata",sep=""))
  }
}
pwarning <- function(...) {
  cat(paste("[WARNING] ",...,"\n",sep=""),file=stderr())
}
pmissing <- function(...) {
  cat(paste("[MISSING FILE] ",...,"\n",sep=""))
}

perror <- function(...) {
  cat(paste("[ERROR] ",...,"\n",sep=""),file=stderr())
}
write.tsv <- function(x,file,header=TRUE,rownames.label=NULL) {
  if (!is.null(x)) {
	if (!is.matrix(x)) {
	  x <- as.matrix(x)
	}
	x <- apply(x,c(1,2),gsub,pattern="'",replacement="")
	if ( !is.null(rownames.label) ) {
	  y <- cbind(rownames(x),x)
	  colnames(y) <- append(rownames.label,colnames(x))
	  x <- y
	}
  }
  write.table(x,file,sep="\t",row.names=F,col.names=header,quote=F)
  return(1)
}

gen.plot2report <- function(filename=NULL,
					 dir=NULL,
					 bg="white",
					 width=400,
					 height=400,
					 size.unit="px",
					 html=TRUE,
					 ps=TRUE,
					 data.table=NA,
					 caption="",
					 to.plot=NULL) {

  if ( is.null(dir) ) {
	if ( grepl("^/",filename) ) {
	  dir <- ""
	} else {
	  dir <- "."
	}
  }
  if ( dir != "" ) {
	dir <- paste(dir,"/",sep="")
  }  
  # automatic filename  
  if ( is.null(filename)) {
	time.label <- format(Sys.time(), "%d%m%Y%H%M%S")
	filename <- paste("graph_",time.label,".png",sep="")
  }
  ps.file <- NA
  data.file <- NA
  png.file <- paste(dir,filename,".png",sep="")
  png.file <- sub(".png.png$",".png",png.file)
  # width and height is in pixels
  # convert to in
  if ( size.unit == "px" ) {
	# width: 600px=20cm=7.8in
	width.in <- width*7.8/600
	height.in <- height*7.8/600 
  } else {
	width.in <- width
	height.in <- height
  }
  # PNG
  png(filename=png.file,width=width.in, height=height.in,bg=bg,res=300,units="in")
  # Catch the errors
  err <- try(to.plot())
  dev.off()
  if ( class(err)=="try-error") {
	pwarning("Failed to generate plot ",png.file)
	return(NULL);
  }
  # PS
  if (ps==TRUE) {
	ps.file <- paste(png.file,".eps",sep="")
	#ps(file=ps.file)
	postscript(file=ps.file,fonts=c("serif", "sans"))
	err <- try(to.plot())
	dev.off()
	if ( class(err)=="try-error") {
	  pwarning("Failed to generate plot ",png.file)
	  return(NULL);
	}
  }
  # DATA
  if ( sum(!is.na(data.table)) > 0 ) {
	data.file <- paste(png.file,".csv",sep="")
	write.tsv(data.table,data.file)
  }
  png.file
  # HTML
  if ( html == TRUE ) {
	html <- paste("<DIV class=\"dplot\"><IMG src='",basename(png.file),"' border =0 width=",width," height=",height,"><BR/>",sep="")
	files <- c()
	if (ps==TRUE) {
	  files <- append(files,basename(ps.file))
	  names(files)[length(files)] <- "PS"
	}
	if (!is.na(data.file)) {
	  files <- append(files,basename(ps.file))
	  names(files)[length(files)] <- "TSV"
	}
	#
	html <- paste(html,html.download.bar(files),"</DIV>",sep="")
	html <- paste(html,"<p class='caption'>",caption,"</p>",sep="")
  } else {
	HTML=FALSE
  }
  return(list(png.file=png.file,html=html,data.file=data.file,ps.file=ps.file))
}
#
#############################################################
#main program
usage <- "Rprogram --tsv file --go file [options]"
option_list <- list(
  make_option(c("-i", "--tsv"), type="character", dest="express.file", default=NULL,help="TSV file (gene ids should appear in the first column)"),
  make_option(c("-p", "--pvalue-col"), type="numeric", dest="pval.col", default=8,help="Column with the p-values [default %default]"),
  make_option(c("-f", "--foldchange-col"), type="numeric", dest="fold.col", default=5,help="Column with the foldchange values [default %default]"),
  make_option(c("-l", "--pvalue"), type="numeric", dest="pvalue", default=0.05,help="p-value [default %default]"),
  make_option(c("-s", "--minsize"), type="numeric", dest="min.size", default=5,help="Minimum size of a cluster [default %default]"),
  make_option(c("-x", "--maxsize"), type="numeric", dest="max.size", default=+Inf,help="Maximum size of a cluster [default %default]"),
  make_option(c("-m", "--method"), type="character", dest="method", default="median",help="Method [default: %default]"),
  make_option(c("-e", "--minedge"), type="numeric", dest="minedge", default=3,help="Minimum number of genes (>1) shared between two gene sets (used only in the plots) [default: %default]"),
  make_option(c("-t", "--top"), type="numeric", dest="top", default=Inf,help="Top gene sets (>1) used in the plot [default: %default]"),
  make_option(c("--go"), type="character", dest="go.file", default=NULL,help="TSV file with the mapping of the genes to GO terms (two columns file: geneid goterm)"),
  make_option(c("--title"), type="character", dest="title", default="",help="Main title to include in the plot"),
  make_option(c("-o", "--out"), type="character",default=NULL,help="Output file name prefix.")
)
##### Parse multiple options #####
multiple.options = list(method=c("mean","median","sum","fisher","stouffer","tailStrength","wilcoxon","reporter","page"))
filenames <- c("express.file","go.file") ;

mandatory <- c("express.file","out")
opt <- myParseArgs(usage=usage, option_list=option_list,filenames.exist=filenames,multiple.options=multiple.options,mandatory=mandatory)

pinfo("Parameters parsed.")

##############################
pinfo("Reading ",opt$express.file)
df <- read.csv(opt$express.file,sep="\t",header=T,quote="\"")
# Filter duplicate entries
if ( length(unique(df[,1]))!=nrow(df) ) {
  pinfo("Found duplicated gene entries in ",opt$tsv.file,".")
  dup.genes <- unique(as.character(df[duplicated(df[,1]),1]))
  pinfo("Picking the best statistic for each set of entries related to a gene")
  pinfo("Starting with ",nrow(df)," entries")
  for ( g in dup.genes ) {
	to.exclude <- which(df[,1]==g)
	pick <- which(df[to.exclude,opt$pval.col]==min(df[to.exclude,opt$pval.col]))[1]
	to.exclude <- to.exclude[-pick]
	df <- df[-to.exclude,] 
  }
  pinfo("No more duplicates")
  pinfo("",nrow(df)," unique entries")
}
fc <- df[,opt$fold.col]
names(fc) <- df[,1]
pval <- as.numeric(df[,opt$pval.col])
names(pval) <- df[,1]
pinfo("Data in ",opt$tsv.file," loaded: ",nrow(df)," entries")
##############################
if (!is.null(opt$go.file)) {
  pinfo("Reading ",opt$go.file)
  gene2go <- read.table(opt$go.file,sep="\t",header=F,comment.char="!",quote="\"")
  pinfo("Data in ",opt$go.file," loaded: ",nrow(gene2go)," entries")
  colnames(gene2go) <- c("gene","go.id")
}
out.tsv.file <- paste(opt$out,".csv",sep="")

##############################################
# print options
pinfo("Options")
pinfo("method=",opt$method)
pinfo("p-value=",opt$pvalue)
pinfo("top=",opt$top)
pinfo("minedge=",opt$minedge)
pinfo("minimum size=",opt$min.size)
pinfo("maximum size=",opt$max.size)
pinfo("annotation col=",opt$annotation_col)
pinfo("...............................")
##############################################

pval <- pval[!is.na(pval)]
fc <- fc[!is.na(fc)]
pinfo("FC mean ",round(mean(abs(fc),na.rm=T),2),"+-",round(sd(fc,na.rm=T),2))
pval[is.infinite(pval)] <- 1
pval.lim <- opt$pvalue
pinfo("p-value mean ",round(mean(abs(pval),na.rm=T),2),"+-",round(sd(pval,na.rm=T),2))
# set of genes
set1 <- names(pval[pval<=pval.lim])

set1.annot <- gene2go[gene2go[,1] %in% set1,]
print.empty <- FALSE

minSizeLim <- opt$min.size
maxSizeLim <- opt$max.size
pinfo("#Genes=",length(pval)," #genes (pvalue<=",pval.lim,")=",length(set1))

if ( length(set1)>=minSizeLim && print.empty==FALSE )  {
  print.empty <- FALSE
  myGSC <- loadGSC(set1.annot)
  gsaRes <- try(runGSA(geneLevelStats=pval,directions=fc,geneSetStat=opt$method,
					   gsc=myGSC,gsSizeLim=c(minSizeLim,maxSizeLim)))
  if ( class(gsaRes)=="try-error") {
	perror("Failed runGSA(",opt$method,")")
	pinfo("Probably no gene sets were found.")
	print.empty <- TRUE
  } else {
	 gsaRes2plot <- gsaRes
	 sel <- NULL
	   if ( !is.infinite(opt$top) ) {
		 v1 <- as.numeric(gsaRes2plot$pAdjDistinctDirDn)
		 v2 <- as.numeric(gsaRes2plot$pAdjDistinctDirUp)
		 names(v1) <- names(gsaRes2plot$gsc)
		 names(v2) <- names(gsaRes2plot$gsc)
		 v <- apply(cbind(v1,v2),1,min,na.rm=TRUE)
		 # select the top gene sets (lowest p-values, larger number of genes)
		 v2 <- as.numeric(gsaRes2plot$nGenesTot)
		 m <- cbind(v,v2)       
		 rownames(m) <- names(gsaRes2plot$gsc)
		 sorted <- m[order(v,-v2),]
		 sel <- head(rownames(sorted),opt$top)
	   }
	   if ( ! opt$method %in% c("fisher","stouffer","reporter","tailStrength")) {
		 gen.plot2report(filename=paste(opt$out,"_class_distinct_dir_both.png",sep=""),
						 width=800,
						 height=800,
						 html=FALSE,
						 ps=TRUE,
						 to.plot = function () {
						 	nw1 <- networkPlot(gsaRes2plot,class="distinct",direction="both",edgeWidth=c(2,20),adjusted=T,significance=opt$pvalue,geneSets=sel,overlap=opt$minedge,cexLegend=0.8,main=opt$title)
						 }                   
				)
	   }
	   sel <- NULL
	   if ( !is.infinite(opt$top) ) {
		# select the top gene sets (lowest p-values, larger number of genes)
		 v <- as.numeric(gsaRes2plot$pAdjNonDirectional)
		 v2 <- as.numeric(gsaRes2plot$nGenesTot)
		 m <- cbind(v,v2)       
		 rownames(m) <- names(gsaRes2plot$gsc)
		 sorted <- m[order(v,-v2),]
		 sel <- head(rownames(sorted),opt$top)
	   }
	   gen.plot2report(filename=paste(opt$out,"_class_non_dir_both.png",sep=""),
					   width=550,
					   height=550,
					   html=FALSE,
					   ps=TRUE,
					   to.plot = function () {
						 nw1 <- networkPlot(gsaRes2plot,class="non",edgeWidth=c(2,20),adjusted=T,significance=opt$pvalue,geneSets=sel,overlap=opt$minedge,cexLegend=0.8,main=opt$title)
					   }                     
					   )
	 write.tsv(GSAsummaryTable(gsaRes),file=out.tsv.file,header=T)
	 pinfo("Created ",out.tsv.file)
	 print.empty <- FALSE
   }
}

if(print.empty) {
  gsaRes=NA
  write.tsv(matrix(),file=out.tsv.file)
  pinfo("Created empty ",out.tsv.file)
}
save(gsaRes,file=paste(opt$out,".Rdata",sep=""))

quit(save="no",status=0)

