import concurrent.futures
import os
import sys
import json
from traplfunlib.paths import Paths
from traplfunlib.projectcreator import ProjectCreator
from traplfunlib.retrievego import RetrieveGo
from traplfunlib.gff3 import Gff3Parser
from traplfunlib.go_enrich import goenrichanalysis
from traplfunlib.retrievepathway import Retrievepathway
from traplfunlib.pathway_enrich import Pathway_analysis
from traplfunlib.blast2go import blast2go_analysis
from traplfunlib.go_viz import Goviz
from traplfunlib.gene_set_analysis import GSEA_analysis
from traplfunlib.clustering import cluster_analysis
from traplfunlib.pathway_viz import Pathviz

class Controller(object):
    """Control all the commands."""

    def __init__(self, args):
        self._args = args
        self._paths = Paths(args.project_path)

    def create_project(self, version):
        project_creator = ProjectCreator()
        project_creator.create_root_folder(self._args.project_path)
        project_creator.create_subfolders(self._paths.required_folders())
        project_creator.create_version_file(self._paths.version_path,version)
        sys.stdout.write("Created folder \"%s\" and required subfolders. \n" %(
            self._args.project_path))
    
    def blast2go(self):
        self._test_folder_existance(
            self._paths.required_input_folders())
        blast2go_file = self._paths.get_blast2go_files()
        blast2go_path = self._paths.set_blast2go_id_paths(blast2go_file)
        blast2go = blast2go_analysis()
        for blast2go_single_file in blast2go_path:
            blast2go.run_blast2go(self._paths.blast2go_path,
                    blast2go_single_file, self._paths.go_background_list_path)
        
    def go_terms(self):
        self._test_folder_existance(
            self._paths.required_go_folders())
        uniprot = self._args.uniprot
        annotation_files = self._paths.get_annotation_files()
        annotation_paths = self._paths.set_annotation_paths(annotation_files)
        go_creator = RetrieveGo()
        name_list = set()
        gffparser = Gff3Parser()
        for annotation_single_file in annotation_paths:
            for entry in gffparser.entries(open(annotation_single_file)):
                if entry.feature == 'CDS':
                    try:
                        ref_name = entry.attributes["Name"]
                        if ref_name not in name_list:
                            name_list.add(ref_name)
                    except:
                        ref_name = ''

        open(self._paths.go_background_list_path,"w").close()
        go_creator.retrieve_go(name_list,
                               uniprot,
                               self._paths.go_background_list_path)

    def kegg_terms(self):
        self._test_folder_existance(
            self._paths.required_kegg_folders())
        retrieve_pathway = Retrievepathway()
        retrieve_pathway.buildpathway(self._paths.common_r_path,
                        self._args.code,
                        self._paths.kegg_background_list_path)
        
        
    def go_stat(self):
        target_id_files = self._paths.get_target_id_files()
        target_id_paths = self._paths.set_target_id_paths(target_id_files)
        self._test_folder_existance(target_id_paths)
        gsea = self._args.gsea
        fdr = self._args.gsea
        go_obo = self._args.gobo
        goterm_analysis = goenrichanalysis(gsea,fdr)
        for target_id_each_file in target_id_paths:
            target_id_filename = os.path.basename(target_id_each_file)
            go_enrich_each_file = self._paths.go_enrich_folder + '/' + \
                                        target_id_filename + '_' + \
                                        "go_stat.csv"
            go_gsea_each_file = self._paths.go_enrich_folder + '/'+ \
                                        target_id_filename + '_' + \
                                        "go_gsea_stat.csv"
            open(go_enrich_each_file,"w").close()
            goterm_analysis.go_enrichment(target_id_each_file,
                                    self._paths.go_background_list_path,
                                    go_obo,
                                    go_enrich_each_file,go_gsea_each_file)


    def pathway_stat(self):
        target_id_file = self._paths.get_target_id_files()
        target_id_paths = self._paths.set_target_id_paths(target_id_file)
        self._test_folder_existance(target_id_paths)
        gsea_option = self._args.gsea
        fdr_option = self._args.fdr
        kegg_enrichment = Pathway_analysis(gsea_option,fdr_option)
        for target_id_each_file in target_id_paths:
            target_id_filename = os.path.basename(target_id_each_file)
            kegg_enrich_each_file = self._paths.kegg_enrich_folder + "/" + \
                                    target_id_filename + "_" + \
                                    "path_stat.csv"
            kegg_gsea_each_file = self._paths.kegg_enrich_folder + "/" + \
                                    target_id_filename + "_" + \
                                    "path_gsea_stat.csv"
            open(kegg_enrich_each_file,"w").close()
            open(kegg_gsea_each_file,"w").close()
            kegg_enrichment.pathway_enrichment(target_id_each_file,
                                    self._paths.kegg_background_list_path,
                                    kegg_enrich_each_file, 
                                    kegg_gsea_each_file)
        
    def go_viz(self):
        target_id_files = self._paths.get_target_id_files()
        target_id_paths = self._paths.set_target_id_paths(target_id_files)
        gene_slim = self._args.slim
        go_obo = self._args.gobo
        for target_id_each_file in target_id_paths:
            target_id_filename = os.path.basename(target_id_each_file)
            go_enrich_each_file = self._paths.go_enrich_folder + '/' + \
                                    target_id_filename + '_' + \
                                    "go_stat.csv"
            go_viz_each_file = self._paths.go_viz_folder + '/' + \
                                    target_id_filename + '_' + \
                                    "go_viz_revigo.txt"
            go_viz_file_tag = self._paths.go_viz_folder + '/' + \
                                    target_id_filename
            open(go_viz_each_file,"w").close()
            go_viz_object = Goviz(go_obo, \
                        gene_slim, \
                        self._paths.common_r_path)
            go_viz_object.go_viz(go_enrich_each_file,
                            go_viz_each_file, go_viz_file_tag)
    
    def path_viz(self):
        target_id_file = self._paths.get_target_id_files()
        target_id_paths = self._paths.set_target_id_paths(target_id_file)
        self._test_folder_existance(target_id_paths)
        for target_id_each_file in target_id_paths:
            target_id_filename = os.path.basename(target_id_each_file)
            kegg_gsea_each_file = self._paths.kegg_enrich_folder + "/" + \
                                    target_id_filename + "_" + \
                                    "path_gsea_stat.csv"
            path_viz_object=Pathviz(self._paths.common_r_path, \
                kegg_gsea_each_file, target_id_each_file, self._args.code, \
                self._paths.kegg_viz_folder)
            path_viz_object.path_viz()
        
        
    def clustering(self):
        expression_files = self._paths.get_expression_files()
        expression_paths = self._paths.set_expression_paths(expression_files)
        row_method = self._args.rowmethod
        row_matrix = self._args.rowmatrix
        col_method = self._args.colmethod
        col_matrix = self._args.colmatrix
        color = self._args.color
        for expression_single_file in expression_paths:
            expression_filename = os.path.basename(expression_single_file)
            expression_filepath = self._paths.clustering_base_folder + "/" + \
                                expression_filename
            cluster_object=cluster_analysis(expression_single_file, \
                            self._paths.clustering_base_folder, row_method, \
                            row_matrix, col_method, col_matrix, color)
            cluster_object.cluster()
                
    def gsea(self):
        expression_files = self._paths.get_expression_files()
        expression_paths = self._paths.set_expression_paths(expression_files)
        method = self._args.method
        minsize = self._args.minsize
        minedge = self._args.minedge
        pvalue = self._args.pcol
        foldchange = self._args.fcol
        gsea_object=GSEA_analysis(pvalue,foldchange,method,minsize,minedge)
        for expression_single_file in expression_paths:
            expression_filename = os.path.basename(expression_single_file)
            expression_filepath = self._paths.gsea_base_folder + "/" + \
                                expression_filename
            gsea_object.gsea(self._paths.gsa_path, expression_single_file, \
                self._paths.go_background_list_path, \
                self._paths.go_ontology_obo_path, expression_filepath)
            
    def _test_folder_existance(self, task_specific_folders):
        for folder in (
            self._paths.required_base_folders() + task_specific_folders):
            if not os.path.exists(folder):
                self._write_err_msg_and_quit(
                    "Error! Folder '%s' does not exist! Is the given project "
                    "folder name correct?\n" % folder)


