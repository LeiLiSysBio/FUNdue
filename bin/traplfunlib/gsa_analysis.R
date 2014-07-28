#!/usr/bin/env Rscript

#########################
# PIANO
# source("http://bioconductor.org/biocLite.R")
# biocLite("piano")
suppressPackageStartupMessages(library("piano"))

suppressPackageStartupMessages(library("optparse"))

myParseArgs <- function(usage,option_list,filenames.exist=NULL,multiple.options=NULL,mandatory=NULL) {

  # get command line options, if help option encountered print help and exit,
  # otherwise if options not found on command line then set defaults,
  parser <- OptionParser(usage = usage, option_list=option_list)
  opt <- parse_args(parser)

  for ( m in mandatory ) {
	if ( is.null(opt[[m]]) ) {
		perror("Parameter ",m," needs to be defined")
		q(status=1)
	}
  }

  for ( p in filenames.exist ) {
	if (! is.null(opt[[p]]) ) {
	  if (! file.exists(opt[[p]]) ) {
		perror("File ",opt[[p]]," not found")
		q(status=1)
	  }
	}
  }

  for ( op in names(multiple.options) ) {
	if ( ! opt[[op]] %in% multiple.options[[op]] ) {
	  perror("Invalid value ",opt[[op]]," for option ",op)
	  q(status=1)
	}
  }
  return(opt)
}

perror <- function(...) {
  cat(paste("[ERROR] ",...,"\n",sep=""),file=stderr())
}
pinfo <- function(...) {
  cat(paste("[INFO] ",...,"\n",sep=""))
}
read.tsv <- function(file,header=T) {
  read.csv(file,sep="\t",header,quote="\"")
}
pdebug <- function(...) {
  if ( pdebug.enabled) {
	cat(paste("[DEBUG] ",...,"\n",sep=""),file=stderr())
  }
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
  #
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
# generates a png file and a PDF file
# return the HTML to display the image
# with (optionally) links to the pdf
# and data (TSV) file
# returns a list with the filenames for png, pdf, data, and
# the html code
# (png.file=,pdf.file=,data.file=,html=)
# to.plot - function to make the plot
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
	  # abspath
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
  #print(err)
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
  make_option(c("-a", "--annotation"), type="character",default=NULL,help="Annotation file (TSV format)"),
  make_option(c("--annotation_col"), type="character",default="GOterm",help="Column in the annotation file [default %default]"),
  make_option(c("-i", "--tsv"), type="character", dest="tsv.file", default=NULL,help="TSV file (gene ids should appear in the first column)"),
  make_option(c("-p", "--pvalue-col"), type="numeric", dest="pval.col", default=8,help="Column with the p-values [default %default]"),
  make_option(c("-f", "--foldchange-col"), type="numeric", dest="fold.col", default=5,help="Column with the foldchange values [default %default]"),
  make_option(c("-l", "--pvalue"), type="numeric", dest="pvalue", default=0.05,help="p-value [default %default]"),
  make_option(c("-s", "--minsize"), type="numeric", dest="min.size", default=5,help="Minimum size of  a cluster [default %default]"),
  make_option(c("-x", "--maxsize"), type="numeric", dest="max.size", default=+Inf,help="Maximum size of  a cluster [default %default]"),
  make_option(c("-m", "--method"), type="character", dest="method", default="median",help="Method [default: %default]"),
  make_option(c("-e", "--minedge"), type="numeric", dest="minedge", default=3,help="Minimum number of genes shared between two gene sets (used only in the plots) [default: %default]"),
  make_option(c("-t", "--top"), type="numeric", dest="top", default=Inf,help="Top gene sets used in the plot [default: %default]"),
  make_option(c("--annot-only"),action="store_true",dest="annot.only",default=FALSE,help="Consider only genes with annotation [default %default]"),
  make_option(c("--plot-annot-only"),action="store_true",dest="plot.annot.only",default=FALSE,help="Plot  only genes with annotation [default %default]"),
  make_option(c("--dup-use-best"),action="store_true",dest="dup.use.best",default=TRUE,help="Pick the best statistic (lowest p-value) for duplicated genes [default %default]"),
  make_option(c("--debug"),action="store_true",dest="debug",default=FALSE,help="Debug mode (save a Rdata file for debugging purposes)"),
  make_option(c("--go"), type="character", dest="go.file", default=NULL,help="TSV file with the mapping of the genes to GO terms (two columns file: geneid goterm)"),
  make_option(c("--title"), type="character", dest="title", default="",help="Main title to include in the plot"),
  make_option(c("-o", "--out"), type="character",default=NULL,help="Output file name prefix.")
)
multiple.options = list(method=c("mean","median","sum","fisher","stouffer","tailStrength","wilcoxon","reporter","page"))
# required filenames
filenames <- c("tsv.file","go.file","annotation") ;

# check multiple options values
mandatory <- c("tsv.file","out")
opt <- myParseArgs(usage = usage, option_list=option_list,filenames.exist=filenames,multiple.options=multiple.options,mandatory=mandatory)

pdebug.enabled <- opt$debug
pinfo("Parameters parsed.")

if ( opt$minedge < 1 ) {
  perror("--minedge value should be greater or equal than 1")
  q(status=1)
}
if ( opt$top < 1 ) {
  perror("--top value should be greater or equal than 1")
  q(status=1)
}
if ( is.null(opt$annotation) && is.null(opt$go.file) ) {
  perror("--go or --annotation needs to be defined.")
  q(status=1)
}
##############################

pinfo("Reading ",opt$tsv.file)
df <- read.tsv(opt$tsv.file,header=T)
if ( opt$pval.col>ncol(df) || opt$pval.col<1 ) {
  perror("Invalid pvalue-col: ",opt$pval.col)
  quit(status=1)
}
if ( opt$fold.col>ncol(df) || opt$fold.col<1 ) {
  perror("Invalid foldchange-col: ",opt$fold.col)
  quit(status=1)
}
#save.image("debug.Rdata")
#load("debug.Rdata")
#####################$
# check for duplicates
if ( length(unique(df[,1]))!=nrow(df) ) {
  pinfo("Found duplicated gene entries in ",opt$tsv.file,".")
  if (!opt$dup.use.best) {
	perror("Please enable --dup-use-best to proceed")
	q(status=1)
  }
  # filter the duplicated genes
  dup.genes <- unique(as.character(df[duplicated(df[,1]),1]))
  pinfo("Picking the best statistic for each set of entries related to a gene")
  pinfo("Starting with ",nrow(df)," entries")
  for ( g in dup.genes ) {
	pdebug("Processing duplicated gene ",g)
	to.exclude <- which(df[,1]==g)
	pdebug("gene has ",length(to.exclude)," entries")
	pick <- which(df[to.exclude,opt$pval.col]==min(df[to.exclude,opt$pval.col]))[1]
	pdebug("Picked ",df[to.exclude[pick],])
	to.exclude <- to.exclude[-pick]
	df <- df[-to.exclude,] 
  }
  pinfo("No more duplicates")
  pinfo("",nrow(df)," unique entries")
}
# 
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
  if (ncol(gene2go)!=2) {
	perror("Error in ",opt$go.file,": expected 2 columns (id goterm)")
	q(status=2)
  }
  colnames(gene2go) <- c("gene","go.id")
} else {
  # load annot file
  pinfo("Loading annotation file ",opt$annotation,"...")
  annot.table <- load.annot(opt$annotation)
  pinfo("Loading annotation file ",opt$annotation,"...done.")
  # check if columns exist
  if ( ! opt$annotation_col %in% colnames(annot.table) ) {
	perror("Selected column ",opt$annotation_col," not found in ",opt$annotation)
	q(status=2)
  }
  pinfo("Generating mapping from gene to goterm/pathway",opt$go.file)
  tmp.m <- as.data.frame(list(ID=annot.table$ID,term=annot.table[,opt$annotation_col]))
  # which terms to use?
  # GOterm
  # KEGG
  #s.tmp.m <- head(tmp.m)
  s.tmp.m <- tmp.m
  x <- strsplit(as.character(s.tmp.m$term),split=";")
  lx <- lapply(x,length)
  names(lx) <- s.tmp.m$ID
  myrep <- function(x,n) {
	rep(x,n[x])
  }

  gene2go <- data.frame(list(gene=unlist(lapply(as.character(s.tmp.m$ID),myrep,lx)),
							 go.id=unlist(x)),stringsAsFactors = F)
}
# output filename
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
pdebug.save.state("irap_GSE_piano","p0")

pinfo("Starting to play piano...")
# create the genesets
# map gene to different 
# sets: 
# exclude NAs
pval <- pval[!is.na(pval)]
fc <- fc[!is.na(fc)]
pinfo("FC mean ",round(mean(abs(fc),na.rm=T),2),"+-",round(sd(fc,na.rm=T),2))
# fc.dir <- fc
# fc.dir[fc.dir<0] <- "down"
# fc.dir[fc.dir>0] <- "up"
# TODO: replace inf by a big value?
pval[is.infinite(pval)] <- 1
pval.lim <- opt$pvalue
pinfo("p-value mean ",round(mean(abs(pval),na.rm=T),2),"+-",round(sd(pval,na.rm=T),2))
# set of genes
set1 <- names(pval[pval<=pval.lim])

# opt$plot.annot.only

if ( opt$annot.only ) {
  # expand set 1
  # genes without annotation
  #print(head(gene2go))
  no.value <- is.na(gene2go[,2]) | gene2go[,2]==""
  to.exclude <- as.character(gene2go[no.value,1])
  #print(length(to.exclude))
  #print(head(set1))
  n1 <- length(set1)
  set1 <- set1[!set1 %in% to.exclude]
  n2 <- length(set1)
  pinfo("Excluded ",n1-n2," genes with no annotation")
}

set1.annot <- gene2go[gene2go[,1] %in% set1,]
print.empty <- FALSE

minSizeLim <- opt$min.size
maxSizeLim <- opt$max.size
pinfo("#Genes=",length(pval)," #genes (pvalue<=",pval.lim,")=",length(set1))

if(nrow(set1.annot)==0 ) {
  err.file <- opt
  ifile <- opt$annotation
  if ( !is.null(opt$go) ) {
	ifile <- opt$go
  }
  pwarning("No gene annotation found...giving up. Please check if you have GO/Pathway information in ",ifile,".")
  print.empty <- TRUE
  q(status=2)
}
#save.image("debug.Rdata")
#pinfo(minSizeLim)
pdebug.save.state("irap_GSE_piano","p1")
if ( length(set1)>=minSizeLim && print.empty==FALSE )  {
  print.empty <- FALSE
  myGSC <- loadGSC(set1.annot)
  # stats - vector fold change with names()=gene name
  #  pass the pval and foldchange

   #if ( opt$method %in% c("fisher","stouffer","reporter","wilcoxon","page") ) {
   #  fc <- NULL
   #}adjMethod="none
   #save.image("debug.Rdata")
  gsaRes <- try(runGSA(geneLevelStats=pval,directions=fc,geneSetStat=opt$method,
					   gsc=myGSC,gsSizeLim=c(minSizeLim,maxSizeLim)))
  pdebug.save.state("irap_GSE_piano","p2")
  if ( class(gsaRes)=="try-error") {
	perror("Failed runGSA(",opt$method,")")
	pinfo("Probably no gene sets were found.")
	print.empty <- TRUE
  } else {
	 # Create the plots
	 # 
	 # which sets to include
	 gsaRes2plot <- gsaRes
	 #opt$plot.annot.only <- T 
	 if (opt$plot.annot.only ) {
	   #Do not plot gene sets without annotation
	   gsaRes2Plot.sel <- (names(gsaRes$gsc)=="" )
	   #gsaRes2plot$gsc[[gsaRes2Plot.sel]] <- NULL
	   #gsaRes2plot$gsc <- gsaRes2plot$gsc[!is.null(gsaRes2plot$gsc)]
	   if ( sum(gsaRes2Plot.sel)==length(names(gsaRes2plot$gsc))) {
		 gsaRes2plot <- NULL
	   } else {
		 if (sum(gsaRes2Plot.sel)>0 ) {
		   gsaRes2plot$gsc <- gsaRes2plot$gsc[!gsaRes2Plot.sel]                
		   for (v in append(grep("^p" ,names(gsaRes2plot),value=T),"nGenesTot") ) {
			 gsaRes2plot[[v]] <- matrix(gsaRes2plot[[v]][!gsaRes2Plot.sel,],ncol=1)
		   }
		 }
	   }
	 }
	 sel <- NULL
	 if (is.null(gsaRes2plot) || length(names(gsaRes2plot$gsc))<3) {
	   pwarning("Unable to generate plots, not enough gene sets")
	 } else {
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
	   #If geneSetStat is set to"fisher","stouffer","reporter" or"tailStrength" only p-values are allowed as geneLevelStats.
	   if ( ! opt$method %in% c("fisher","stouffer","reporter","tailStrength")) {
		 gen.plot2report(filename=paste(opt$out,"_class_distinct_dir_both.png",sep=""),
						 width=450,
						 height=450,
						 html=FALSE,
						 ps=TRUE,
						 to.plot = function () {
						   nw1 <- networkPlot(gsaRes2plot,class="distinct",direction="both",edgeWidth=c(2,20),adjusted=T,significance=opt$pvalue,geneSets=sel,overlap=opt$minedge,cexLegend=0.8,main=opt$title)
						 }                     
						 )
	   }
	 # which sets to include
	   sel <- NULL
	   if ( !is.infinite(opt$top) ) {
										# select the top gene sets (lowest p-values, larger number of genes)
		 v <- as.numeric(gsaRes2plot$pAdjNonDirectional)
		 v2 <- as.numeric(gsaRes2plot$nGenesTot)
										#length(v2)
										#length(names(gsaRes2plot$gsc))
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
	 }
	 # save the tsv file and the rdata file
	 #save(gsaRes,file=paste(opt$out,".Rdata",sep=""),envir=sys.frame())
	 write.tsv(GSAsummaryTable(gsaRes),file=out.tsv.file,header=T)
	 pinfo("Created ",out.tsv.file)
	 print.empty <- FALSE
   }
} else {
  pwarning("No significant entries found for p-value of ",opt$pvalue)
  print.empty <- TRUE
}
pdebug.save.state("irap_GSE_piano","p3")
if(print.empty) {
  gsaRes=NA
  write.tsv(matrix(),file=out.tsv.file)
  pinfo("Created empty ",out.tsv.file)
}
save(gsaRes,file=paste(opt$out,".Rdata",sep=""))

quit(save="no",status=0)

