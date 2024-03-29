#!/usr/bin/env Rscript

library("KEGGREST")
library("getopt")
library("RCurl")
library("png")
suppressPackageStartupMessages(library("pathview"))
suppressPackageStartupMessages(library("gage"))
opt = getopt(matrix( c('help', 'h', 0, "logical", 
                       'verbose', 'v', 0, "integer",
                       'organism', 'o', 1, "character",
                       'database', 'd', 1, "character",
                       'path_id', 'p', 1, "character",
                       'result_path', 'r', 1, "character",
                       'folder', 'f', 1, "character",
                       'go_id', 'g', 1, "character",
                       'color', 'c', 1, "character"
), ncol=4, byrow=TRUE ) );

org <- opt$organism
## perform retrival pathway
if(! is.null(opt$organism) &&  ! is.null(opt$database)){
	getPathwayInfo <- function(pathlist){
		pathnames <- c()
		for( j in 1:length(pathlist)){
			query <- keggGet(pathlist[[j]])
			pathnames <- c(pathnames, query[[1]]$NAME)
		}
		return(cbind(pathlist, pathnames))
	}
	paths <- keggLink("pathway", org)	
	print("Getting the detailed pathway information, normally need 1 min")
	id2name <- getPathwayInfo(unique(paths))	
	pathtable <- cbind(sub(paste(org,":",sep=''), '', names(paths)), sub("path:",'',paths), sapply(paths, function(x) id2name[id2name[,1] == x,2]))	
	write.table(pathtable, file=opt$database, append=FALSE, quote=FALSE, sep="\t", row.names=FALSE, col.names=c('gene_id','path_id','path_name'))
}
## visualization pathway
if(! is.null(opt$organism) && ! is.null(opt$path_id) && ! is.null(opt$result_path) && ! is.null(opt$folder)){
	path <- opt$path_id
	path_result <- opt$result_path
	folder <- opt$folder
	print("Read differetial expression files")
	input <- readExpData(path_result,row.names=NULL)
	cnames=input[,1]
	input=input[,2]
	names(input)=cnames
	input_base <- basename(path_result)
	path_suffix <- paste(input_base,path,sep="_")
	setwd(folder)
	print("Pathway shown")
	pv.out <- pathview(gene.data=input,gene.idtype="kegg",pathway.id=path,species=org,out.suffix=path_suffix,kegg.native=T)
}