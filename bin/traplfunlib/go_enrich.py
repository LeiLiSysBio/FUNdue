import sys
from scipy import stats
import collections
from traplfunlib.obo_parser import GODag
import random


class goenrichanalysis(object):
    """Gene ontolgoy enrichment analysis"""
    def go_enrichment(self,target_id_file, background_file, obo_file_path, go_enrich_out):
        timer = time.clock()
        target_list = []
        background_list = []
        target_no = 0
        association_list = {}
        go_out_file = open(go_enrich_out,"a")

        for entry in open(''.join(target_id_file), "r"):
            uni_line = entry.rstrip("\n")
            uni_lines = uni_line.split("\t")
            if uni_lines not in target_list:
                target_list.append(uni_lines)

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
        
        go_obo = GODag(obo_file_path)
        go_obo.update_association(association_list)

        count_obj = count()

        background_term = count_obj.count_terms(background_list,
                                      association_list, go_obo)

        target_term = count_obj.count_terms(target_list,
                                  association_list, go_obo)

        background_no = len(background_list)
        #target_no = len(target_list)
        go_out_file.write("Gene ontology term"+ "\t" + "ontology description" + "\t"+ "ontologies"\
                          + "\t" + "target number" + "\t" + "total target numbers" \
                          + "\t" + "ratio of  targets" + "\t" + "background_number" \
                          + "\t" + "total background numbers" + "\t" + "ratio of background" \
                          + "\t" + "pvalue" + "\t" + "FDR (False discovery rate)" + "\n")
        pvals = {}
        print("calculating the FDR (false discovery rate) may take for a while, please wait...")
        num = 0
        for term, target_count in target_term.items():
            if term == "GO:0003674|GO:0008150|GO:0005575":
                pass
            else:
                num = num + 1
                print("the "+ str(num) +" term is processing")
                background_count = background_term[term]
                target_other = target_no - target_count
                background_other = background_no - background_count
                oddsratio, pvalue = stats.fisher_exact([[target_count,target_other],
                                        [background_count,background_other]])
                ratio_target = float(target_count) / float(target_no)
                ratio_background = float(background_count) / float(background_no)
                #pvals[term] = pvalue
                         
                ## Calculate FDR (FalseDiscoveryRate) Random select the same number of clusters
                ## And calculate the probability of p less than original one
                distribution = []
                for i in range(500):
                    new_target = random.sample(background_list, target_no)
                    new_term_study = count_obj.count_terms(new_target, association_list, go_obo)
                    min_p = 1
                    if term in new_term_study:
                        new_target_other = target_no - new_term_study[term]
                        new_background_count = background_term[term]
                        new_background_other = background_no - new_background_count
                        new_odd, p = stats.fisher_exact([[new_term_study[term], new_target_other],
                                        [new_background_count, new_background_other]])
                        if p < min_p:
                            min_p = p
                        distribution.append(min_p)
                sum = 0
                for random_p in distribution:
                    if random_p < pvalue:
                        sum = sum + 1
                if len(distribution) == 0:
                    fdr_adjust_p = float(sum)
                else:
                    fdr_adjust_p = float(sum)/len(distribution)
                go_out_file.write(term + "\t" + go_obo[term].name + "\t"+ go_obo[term].namespace +"\t" + str(target_count)
                    + "\t" + str(target_no) + "\t" + str(ratio_target) + "\t"
                    + str(background_count) + "\t" + str(background_no) + "\t"
                    + str(ratio_background) + "\t" + str(pvalue) + "\t"
                    + str(fdr_adjust_p) + "\n")
        print("# %0.2f seconds process time" % (time.clock() - timer))        

class FDR(object):
    def __init__(self, p_val_distribution, results, a=.05):
        self.corrected_pvals = fdr = []
        for rec in results:
            q = (sum(1 for x in p_val_distribution if x < rec.p_uncorrected) \
                * 1.0 / len(p_val_distribution))
            fdr.append(q)
 
class Bonferroni(Correction):
    def __init__(self,pvals, a=0.05):
        self.pvals = self.corrected_pvals = np.array(pvals)
        self.n = len(self.pvals)
        self.a = a
    def set_correction(self):
        self.corrected_pvals *= self.n

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
