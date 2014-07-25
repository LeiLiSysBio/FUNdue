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
from traplfunlib.gsea import GSEA_analysis


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
                               self._paths.uniprot_id_mapping_path,
                               self._paths.go_background_list_path)

    def kegg_terms(self):
        self._test_folder_existance(
            self._paths.required_kegg_folders())
        retrieve_pathway = Retrievepathway()
        retrieve_pathway.buildpathway(self._paths.build_pathway_path,
                        self._args.code,
                        self._paths.kegg_background_list_path)
        
        
    def go_stat(self):
        target_id_files = self._paths.get_target_id_files()
        target_id_paths = self._paths.set_target_id_paths(target_id_files)
        self._test_folder_existance(target_id_paths)
        goterm_analysis = goenrichanalysis()

        open(self._paths.go_enrich_list_path,"w").close()
        for target_id_each_file in target_id_paths:
            target_id_filename = os.path.basename(target_id_each_file)
            go_enrich_each_file = self._paths.go_enrich_folder + '/' + \
                                        target_id_filename + '_' + \
                                        "go_stat.txt"
            open(go_enrich_each_file,"w").close()
            goterm_analysis.go_enrichment(target_id_each_file,
                                    self._paths.go_background_list_path,
                                    self._paths.go_ontology_obo_path,
                                    go_enrich_each_file)


    def pathway_stat(self):
        target_id_file = self._paths.get_target_id_files()
        target_id_paths = self._paths.set_target_id_paths(target_id_file)
        self._test_folder_existance(target_id_paths)
        kegg_enrichment = Pathway_analysis()
        
        open(self._paths.kegg_enrich_list_path,"w").close()
        kegg_enrichment.pathway_enrichment(target_id_paths,
                                            self._paths.kegg_background_list_path,
                                            self._paths.kegg_enrich_list_path)
        
    def go_viz(self):
        go_viz_object = Goviz()
        open(self._paths.go_viz_list_path,"w").close()
        target_id_files = self._paths.get_target_id_files()
        target_id_paths = self._paths.set_target_id_paths(target_id_files)
        for target_id_each_file in target_id_paths:
            target_id_filename = os.path.basename(target_id_each_file)
            go_enrich_each_file = self._paths.go_enrich_folder + '/' + \
                                    target_id_filename + '_' + \
                                    "go_stat.txt"
            go_viz_each_file = self._paths.go_viz_folder + '/' + \
                                    target_id_filename + '_' + \
                                    "go_viz_revigo.txt"
            open(go_viz_each_file,"w").close()
            go_viz_object.go_viz(go_enrich_each_file,
                            go_viz_each_file)
    
    def gsea(self):
        gsea_object=GSEA_analysis()
        expression_files = self._paths.get_expression_files()
        expression_paths = self._paths.set_expression_paths(expression_files)
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


