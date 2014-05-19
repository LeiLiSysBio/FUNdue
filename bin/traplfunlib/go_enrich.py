import sys
from scipy import stats
import collections
from traplfunlib.obo_parser import GODag


class goenrichanalysis(object):
    """Gene ontolgoy enrichment analysis"""
    def go_enrichment(self,target_id_file, background_file, obo_file_path, go_enrich_out):

        target_list = []
        background_list = []
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

        go_obo = GODag(obo_file_path)
        go_obo.update_association(association_list)

        count_obj = count()

        background_term = count_obj.count_terms(background_list,
                                      association_list, go_obo)

        target_term = count_obj.count_terms(target_list,
                                  association_list, go_obo)

        background_no = len(background_list)
        target_no = len(target_list)
        go_out_file.write("Gene ontology term"+ "\t" + "ontology description" \
                          + "\t" + "target number" + "\t" + "total target numbers" \
                          + "\t" + "ratio of  targets" + "\t" + "background_number" \
                          + "\t" + "total background numbers" + "\t" + "ratio of background" \
                          + "\t" + "pvalue" +"\n")
        for term, target_count in target_term.items():
            background_count = background_term[term]
            target_other = target_no - target_count
            background_other = background_no - background_count
            oddsratio, pvalue = stats.fisher_exact([[target_count,target_other],
                                                    [background_count,background_other]])
            ratio_target = float(target_count) / float(target_no)
            ratio_background = float(background_count) / float(background_no)
            go_out_file.write(term + "\t" + go_obo[term].name + "\t" + str(target_count)
                + "\t" + str(target_no) + "\t" + str(ratio_target) + "\t"
                + str(background_count) + "\t" + str(background_no) + "\t"
                + str(ratio_background) + "\t" + str(pvalue) + "\n")

class count(object):
    def count_terms(self, geneset, assoc, obo_dag):
      term_cnt = collections.defaultdict(int)
      self.assoc = assoc
      for g in geneset:
        g_str = ''.join(g)
        if g_str in assoc:
            for x in assoc[g_str]:
                if x in obo_dag:
                    term_cnt[obo_dag[x].id] += 1

      return term_cnt
