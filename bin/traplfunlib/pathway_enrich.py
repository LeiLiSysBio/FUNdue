import sys
from scipy import stats
import collections
from traplfunlib.gff3 import Gff3Parser
from traplfunlib.paths import Paths

class Pathway_analysis(object):
	# perform pathway enrichment analysis
	def pathway_enrichment(self,target_id_file, background_file, pathway_enrich_out):
		target_list = []
		background_list = []
		associate_id_map = {}
		associate_map_desc = {}
		pathway_enrich_file = open(pathway_enrich_out,"a")
		
		pathway_out_file = open(pathway_enrich_out,"a")
		
		for entry in open(''.join(target_id_file), "r"):
			uni_line = entry.rstrip("\n")
			uni_lines = uni_line.split("\t")
			if uni_lines not in target_list:
				target_list.append(uni_lines)
		
		with open(background_file,"r") as back:
			next(back)
			for entry in back:
				uni_line = entry.rstrip("\n")
				uni_lines = uni_line.split("\t")
				associate_id_map[uni_lines[0]] = uni_lines[1]
				associate_map_desc[uni_lines[1]] = uni_lines[2]
				if uni_lines[0] not in background_list:
					background_list.append(uni_lines[0])
		
		count_obj = count()
		background_term = count_obj.count_terms(background_list,
											associate_id_map)
		target_term = count_obj.count_terms(target_list,
											associate_id_map)
		
		background_no = len(background_list)
		target_no = len(target_list)
		pathway_enrich_file.write("KEGG pathway term"+ "\t" + "pathway description" \
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
			pathway_enrich_file.write(term + "\t" + associate_map_desc[term] + "\t" \
				+ str(target_count) + "\t" + str(target_no) + "\t" \
				+ str(ratio_target) + "\t" + str(background_count) + "\t" \
				+ str(background_no) + "\t" + str(ratio_background) + "\t" + str(pvalue) + "\n")
		
class count(object):
	def count_terms(self, geneset, assoc):
		term_cnt = collections.defaultdict(int)
		self.assoc = assoc
		for 	each_gene in geneset:
			each_gene_str = ''.join(each_gene)
			if each_gene_str in assoc:
				term_cnt[assoc[each_gene_str]] += 1
					
		return term_cnt
			
		
			
		
		
		