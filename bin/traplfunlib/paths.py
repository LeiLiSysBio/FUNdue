import os
import sys

class Paths(object):
    def __init__(self, base_path="."):
        self.base_path = base_path
        self._set_folder_names()
        self._set_static_files()

    def _set_folder_names(self):
        self.input_folder = "%s/input" % (self.base_path)
        self.output_folder = "%s/output" % (self.base_path)
        self._set_input_folder_names()
        self._set_go_folder_names()
        self._set_kegg_folder_names()
        self._set_gsea_folder_names()
        self._set_clustering_folder_names()

    def _set_input_folder_names(self):
        self.annotation_folder = "%s/annotations" % self.input_folder
        self.target_id_folder = "%s/target_ids" % self.input_folder
        self.blast2go_folder = "%s/blast2go_xml" % self.input_folder
        #expression could be DESeq result
        self.expression_folder = "%s/expression" % self.input_folder

    def _set_go_folder_names(self):
        self.go_base_folder = "%s/gene_ontology" % self.output_folder
        self.go_background_folder = "%s/background_go/" % self.go_base_folder
        self.go_enrich_folder = "%s/go_enrichment/" % self.go_base_folder
        self.go_viz_folder = "%s/go_viz/" % self.go_base_folder
        self.go_background_list_path = (
            "%s/background_go.csv" %
            self.go_background_folder)
        self.go_enrich_list_path = (
            "%s/go_term_stat.csv" %
            self.go_enrich_folder)
        self.go_viz_list_path = (
            "%s/go_viz_revigo.txt" %
            self.go_viz_folder)

    def _set_kegg_folder_names(self):
        self.kegg_base_folder = "%s/pathway" % self.output_folder
        self.kegg_background_folder = "%s/background_path/" % self.kegg_base_folder
        self.kegg_enrich_folder = "%s/pathway_enrichment/" % self.kegg_base_folder
        self.kegg_viz_folder = "%s/path_viz/" %self.kegg_base_folder
        self.kegg_background_list_path = (
            "%s/background_kegg.csv" %
            self.kegg_background_folder)
        
    def _set_clustering_folder_names(self):
        self.clustering_base_folder = "%s/clustering_analysis/" % self.output_folder
             
    def _set_gsea_folder_names(self):
        self.gsea_base_folder = "%s/gene_set_analysis" % self.output_folder
         
    def _set_static_files(self):
        self.version_path = "%s/used_trapl_fun_version.txt" % (self.output_folder)
        self.uniprot_id_mapping_path = "/storage2/lei/Project/TRAPL_FUN/UniprotKB_GO/idmapping_selected.tab"
        self.go_ontology_obo_path = "/storage2/lei/Project/TRAPL_FUN/gene_ontology.1_2.obo"
        self.slim_go_ontology_obo_path = "/storage2/lei/Project/TRAPL_FUN/goslim_generic.obo"
        self.common_r_path = "bin/traplfunlib/common.R"
        self.blast2go_path = "bin/b2g4pipe"
        self.gsa_path = "bin/traplfunlib/gsa_analysis.R"

    def required_folders(self):
        return(self.required_base_folders() +
               self.required_input_folders() +
               self.required_go_folders() +
               self.required_kegg_folders() + 
               self.required_gsea_folders() +
               self.required_clustering_folder())

    def required_base_folders(self):
        return [self.input_folder, self.output_folder]

    def required_input_folders(self):
        return [self.annotation_folder, self.target_id_folder,
                self.blast2go_folder, self.expression_folder]

    def required_go_folders(self):
        return [self.go_base_folder,self.go_enrich_folder,
                self.go_viz_folder, self.go_background_folder]
    
    def required_kegg_folders(self):
        return [self.kegg_base_folder, self.kegg_background_folder,
                self.kegg_enrich_folder, self.kegg_viz_folder]
    
    def required_clustering_folder(self):
        return [self.clustering_base_folder]
        
    def required_gsea_folders(self):
        return [self.gsea_base_folder]
        
    def _get_sorted_folder_content(self, folder):
        "Return the sorted file list of a folder"
        return list(filter(lambda file:
                           not (file.endswith("~") or
                                os.path.basename(file).startswith(".")),
                           sorted(os.listdir(folder))))

    def get_annotation_files(self):
        """Read the names of annotation files."""
        return self._get_sorted_folder_content(self.annotation_folder)
    
    def get_expression_files(self):
        """Read the names of expression files."""
        return self._get_sorted_folder_content(self.expression_folder)
    
    def get_target_id_files(self):
        """Read the names of target id files."""
        return self._get_sorted_folder_content(self.target_id_folder)
    
    def get_blast2go_files(self):
        return self._get_sorted_folder_content(self.blast2go_folder)
    
    def set_annotation_paths(self, annotation_files):
        return self._path_list(
            self.annotation_folder, annotation_files)
            
    def set_expression_paths(self, expression_files):
        return self._path_list(
            self.expression_folder, expression_files)

    def set_target_id_paths(self, target_id_files):
        return self._path_list(
            self.target_id_folder,target_id_files)
            
    def set_blast2go_id_paths(self, blast2go_id_files):
        return self._path_list(
            self.blast2go_folder, blast2go_id_files)

    def _path_list(self, folder, files, appendix=""):
        return ["%s/%s%s" % (folder, file, appendix) for file in files]
