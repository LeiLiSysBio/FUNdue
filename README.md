# README #


### What is this repository for? ###

* Quick summary
Mainly used for functional annotation of transcriptome data, but it can be broadly used in any other high-throughput sequencing data. 

### Waiting list ###

Welcome for your comments

## Version 0.5.2 ###
Add the function for visualizing gene ontology trees

### Version 0.5 ###
Add the function for pathway visulization
fix the bugs caused by pathway enrichment 

### Version 0.4.2 ###
Add gene set enrichment function for the pathway and ontology analysis

### Version 0.4.1 ###
Add a few parameters for the clustering and gene set analysis 

### Version 0.4 ###
Add Clustering analysis function
input(expression folder)

### Version 0.3 ###
Add function for Gene Set Analysis

### Version 0.2 ###
Add the output for FDR (false discovery rate)

Remove top gene ontology nodes

Add the ontologies column for output

fix visualisation issues 

### Version 0.15 ###
Allow input multiple target files

Add the go viz result for displaying results in http://revigo.irb.hr/ 

New blast2go function

For the gene ontology analysis, you need paste a list of NCBI RefSeq IDs
(such as YP\_005179941.1)

For the pathway analysis, you need paste a list of locus\_tag IDs
(such as SL1344\_0003)

### How do I get set up? ###

* Dependencies
** Python modules (scipy)
** R modules (KEGGREST,getopt, piano, optparse, gsge, pathview)
** install methods
source("http://bioconductor.org/biocLite.R")
biocLite()


### Who do I talk to? ###

* Talk with me directly
