#!/usr/bin/env Rscript

library("KEGGREST")
library("getopt")
suppressPackageStartupMessages(library("pathview"))
suppressPackageStartupMessages(library("gage"))
opt = getopt(matrix( c('help', 'h', 0, "logical", 
                       'verbose', 'v', 0, "integer",
                       'organism', 'o', 1, "character",
                       'database', 'd', 1, "character",
                       'path_id', 'p', 1, "character",
                       'path_result', 'r', 1, "character",
                       'folder', 'f', 1, "character"
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
	#getGeneInfo <- function(pathlist){
	#	pathnames <- c()
	#	for( j in 1:length(pathlist)){
	#		query <- keggConv("ncbi-geneid", pathlist[[j]])
	#		pathnames <- c(pathnames, query[[1]])
	#	}
	#	return(cbind(pathlist,pathnames))
	#}	
	print("Getting the detailed pathway information, normally need 1 min")
	id2name <- getPathwayInfo(unique(paths))	
	#print("Getting the gene information, normally need 10-15 min")
	#id2geneid <- getGeneInfo(names(paths))	
	#print(id2geneid)
	pathtable <- cbind(sub(paste(org,":",sep=''), '', names(paths)), sub("path:",'',paths), sapply(paths, function(x) id2name[id2name[,1] == x,2]))	
	write.table(pathtable, file=opt$database, append=FALSE, quote=FALSE, sep="\t", row.names=FALSE, col.names=c('gene_id','path_id','path_name'))
}

if(! is.null(opt$organism) && ! is.null(opt$path_id) && ! is.null(opt$path_result) && ! is.null(opt$folder)){
	path <- opt$path_id
	path_result <- opt$path_result
	folder <- opt$folder
	input <- readExpData(path_result,row.names=1)
	#print(input)
	setwd(folder)
	pv.out <- pathview(gene.data=input,gene.idtype="kegg",pathway.id=path,species=org,out.suffix=path,kegg.native=T)
}
