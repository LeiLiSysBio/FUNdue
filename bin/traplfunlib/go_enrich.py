import sys
from scipy import stats
import collections
from traplfunlib.obo_parser import GODag
import random
from traplfunlib.gsea import GSEA
import tempfile
import time
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector
import csv

class goenrichanalysis(object):
    def __init__(self,gsea_option, fdr_option):
        self._gsea_option = gsea_option
        self._fdr_option = fdr_option
        
    """Gene ontolgoy enrichment analysis"""
    def go_enrichment(self,target_id_file, background_file, obo_file_path, go_enrich_out,gsea_out):
        timer = time.clock()
        target_list = []
        background_list = []
        target_no = 0
        association_list = {}
        go_obo = GODag(obo_file_path)
        tmp_gofile = tempfile.NamedTemporaryFile(mode="a",delete=False)
        """gene set enrichment analysis"""
        print("Gene set enrichment analysis")
        if self._gsea_option == 'True':
            print("Processing the gene ontology file, Please wait...")
            for entry in open(background_file, "r"):
                uni_line = entry.rstrip().split("\t")
                if len(uni_line) > 3:
                    if ';' in uni_line[3]:
                        go_list = uni_line[3].replace(" ", "").split(";")
                        for each_go_list in go_list:
                            tmp_gofile.write(uni_line[0] + "\t" + \
                                each_go_list + "\t" + \
                                go_obo[each_go_list].name + "\n")
                    else:
                        tmp_gofile.write(uni_line[0] + "\t" + \
                            uni_line[3] + "\t" + \
                            go_obo[uni_line[3]].name + "\n")
        
            GSEA_analysis = GSEA(tmp_gofile.name, target_id_file, gsea_out,1)
            GSEA_analysis.gsea_analysis() 
        """Normal ontologies analysis"""
        print("Enrichment analysis using Fisher t test")
        for entry in open(''.join(target_id_file), "r"):
            uni_line = entry.rstrip("\n")
            uni_lines = uni_line.split("\t")
            if uni_lines not in target_list:
                target_list.append(uni_lines[0])

        for entry in open(background_file, "r"):
            uni_line = entry.rstrip("\n")
            uni_lines = uni_line.split("\t")
            if uni_lines[3]:
                go_split = set(uni_lines[3].replace(" ", "").split(";"))
                association_list[uni_lines[0]] = go_split
                if uni_lines[0] not in background_list:
                    background_list.append(uni_lines[0])
        
        for entry in target_list:
            entry_str = ''.join(entry)
            if entry_str in background_list:
                target_no = target_no + 1
        
        go_obo.update_association(association_list)
        count_obj = count()
        background_term = count_obj.count_terms(background_list,
                                      association_list, go_obo)
        target_term = count_obj.count_terms(target_list,
                                  association_list, go_obo)
        background_no = len(background_list)
        #target_no = len(target_list)
                
        ##iteratory each gene ontology, the related gene lists
        target_term_list = count_obj.list_terms(target_list,
                                    association_list, go_obo)
        num = 0
        summary = []
        assoc_target_list={}
        for term, target_count in target_term.items():
            if go_obo[term].name == "'cellular_component'|'biological_process'|'molecular_function'":
                pass
            else:
                num = num + 1
                print("the "+ str(num) +" term is processing:\t" + str(term))
                background_count = background_term[term]
                target_other = target_no - target_count
                background_other = background_no - background_count
                oddsratio, pvalue = stats.fisher_exact([[target_count,target_other],
                                        [background_count,background_other]])
                ratio_target = float(target_count) / float(target_no)
                ratio_background = float(background_count) / float(background_no)
                summary.append([term, go_obo[term].name, go_obo[term].namespace, target_count, \
                    target_no, ratio_target, background_count, background_no, \
                    ratio_background, pvalue,target_term_list[str(term)]])
        pval_sum = [x[9] for x in summary]
        adjust_stats = importr('stats')
        p_adjust_fdr = adjust_stats.p_adjust(FloatVector(pval_sum),method='fdr')
        p_adjust_bf = adjust_stats.p_adjust(FloatVector(pval_sum),method='bonferroni')
        for row_no in range(len(summary)):
            summary[row_no].extend([p_adjust_fdr[row_no],p_adjust_bf[row_no]])
        print("Writing to file")
        with open(go_enrich_out,'w',newline='') as csvfile:
            writer = csv.writer(csvfile,delimiter="\t")
            writer.writerow(['Gene ontology term', 'ontology description','ontologies', \
                'target number','total target numbers', 'ratio of  targets', 'background_number', \
                'total background numbers', 'ratio of background', 'pvalue','assoc_target_list','FDR', 'Bonferroni'])
            for line in summary:
                writer.writerow(line)
        print("# %0.2f seconds process time" % (time.clock() - timer))        
        

class count(object):
    def count_terms(self, geneset, assoc, obo_dag):
      term_cnt = collections.defaultdict(int)
      self.assoc = assoc
      for each_gene in geneset:
        each_gene_str = ''.join(each_gene)
        if each_gene_str in assoc:
            for x in assoc[each_gene_str]:
                if x in obo_dag:
                    term_cnt[obo_dag[x].id] += 1

      return term_cnt
    
    def list_terms(self, geneset, assoc, obo_dag):
        term_list = dict()
        self.assoc = assoc
        for each_gene in geneset:
          each_gene_str = ''.join(each_gene)
          if each_gene_str in assoc:
              for x in assoc[each_gene_str]:
                  if x in obo_dag:
                      if obo_dag[x].id in term_list:
                          term_list[obo_dag[x].id].append(each_gene_str)
                      else:
                          term_list[obo_dag[x].id]=[each_gene_str]
        
        return term_list

