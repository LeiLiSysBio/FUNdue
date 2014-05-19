main(){
	PYTHON_PATH=python3.3
	TRAPL_FUN_PATH=bin/traplfun
	TRAPL_FUN_FOLDER=TRAPL_FUN_analysis
	FTP_SOURCE=ftp://ftp.ncbi.nih.gov/genomes/Bacteria/Salmonella_enterica_serovar_Typhimurium_SL1344_uid86645/
	GO_OBO_FILE=/storage2/lei/project/TRAPL_FUN/gene_ontology.1_2.obo
	UNIPROT_IDMAPPING=/storage2/lei/project/TRAPL_FUN/UniprotKB_GO/idmapping_selected.tab

	# create_folders
	# set_up_analysis_folder
	# get_gff_files
	# run_retrieveGO
	 run_go_enrichment_analysis
	# generate_package_to_send
 	# compress_files
}

create_folders(){
    for FOLDER in bin
    do
	if ! [ -d $FOLDER ]
	then
	    mkdir -p $FOLDER
	fi
    done
}

set_up_analysis_folder(){
    if ! [ -d $TRAPL_FUN_FOLDER ]
    then
        $PYTHON_PATH $TRAPL_FUN_PATH create $TRAPL_FUN_FOLDER
    fi
}

get_gff_files(){
    echo "Downloading GFF annotation files"
    wget -cP $TRAPL_FUN_FOLDER/input/annotations/ \
	$FTP_SOURCE/*gff
    echo "----------------------------------------"
    echo "Content of" $TRAPL_FUN_FOLDER/input/annotations/
    ls -l $TRAPL_FUN_FOLDER/input/annotations/
}

run_retrieveGO(){
	$PYTHON_PATH $TRAPL_FUN_PATH \
	retrieve \
	$TRAPL_FUN_FOLDER
}

run_go_enrichment_analysis(){
	$PYTHON_PATH $TRAPL_FUN_PATH \
	go_stat \
	$TRAPL_FUN_FOLDER	
}

generate_package_to_send(){
	SEND_FOLDER=XXX-Functional-annotation-analysis
	TRAPL_FUN=TRAPL_FUN_analysis/
	mkdir ${SEND_FOLDER}
	cp -r ${TRAPL_FUN}/Gene_ontology/ ${SEND_FOLDER}
	zip -r ${SEND_FOLDER}.zip ${SEND_FOLDER}
}

compress_files(){
	find ${SEND_FOLDER} -type f -print0 | xargs -n1 -0 -P24 bzip2
}

main
