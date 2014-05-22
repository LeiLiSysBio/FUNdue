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
        target_id_file = self._paths.get_target_id_files()
        target_id_paths = self._paths.set_target_id_paths(target_id_file)
        self._test_folder_existance(target_id_paths)
        goterm_analysis = goenrichanalysis()

        open(self._paths.go_enrich_list_path,"w").close()
        goterm_analysis.go_enrichment(target_id_paths,
                                      self._paths.go_background_list_path,
                                      self._paths.go_ontology_obo_path,
                                      self._paths.go_enrich_list_path)


    def pathway_stat(self):
        target_id_file = self._paths.get_target_id_files()
        target_id_paths = self._paths.set_target_id_paths(target_id_file)
        self._test_folder_existance(target_id_paths)
        kegg_enrichment = Pathway_analysis()
        
        open(self._paths.kegg_enrich_list_path,"w").close()
        kegg_enrichment.pathway_enrichment(target_id_paths,
                                            self._paths.kegg_background_list_path,
                                            self._paths.kegg_enrich_list_path)
        
    
    def _test_folder_existance(self, task_specific_folders):
        for folder in (
            self._paths.required_base_folders() + task_specific_folders):
            if not os.path.exists(folder):
                self._write_err_msg_and_quit(
                    "Error! Folder '%s' does not exist! Is the given project "
                    "folder name correct?\n" % folder)


