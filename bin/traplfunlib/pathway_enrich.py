from scipy import stats
import collections
from traplfunlib.gsea import GSEA
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector
import csv
import time

class Pathway_analysis(object):
	def __init__(self,gsea_option,fdr_option):
		self._gsea = gsea_option
		self._fdr = fdr_option
		
	def pathway_enrichment(self, target_id_file, background_file, pathway_enrich_out, path_gsea):
		timer = time.clock()
		target_list = []
		background_list = []
		target_no = 0
		associate_id_map = {}
		associate_map_desc = {}
		
		"""gene set enrichment analysis"""
		print("Gene set enrichment analysis")
		GSEA_analysis = GSEA(background_file, target_id_file, path_gsea,1)
		GSEA_analysis.gsea_analysis()
		
		"""Normal pathway analysis"""
		print("Enrichment analysis using Fisher t test")
		for entry in open(target_id_file, "r"):
			uni_line = entry.rstrip("\n")
			uni_lines = uni_line.split("\t")
			if uni_lines not in target_list:
				target_list.append(uni_lines[0])
		with open(background_file,"r") as back:
			next(back)
			for entry in back:
				uni_line = entry.rstrip("\n")
				uni_lines = uni_line.split("\t")
				associate_id_map[uni_lines[0]] = uni_lines[1]
				associate_map_desc[uni_lines[1]] = uni_lines[2]
				if uni_lines[0] not in background_list:
					background_list.append(uni_lines[0])

		for entry in target_list:
			entry_str =''.join(entry)
			if entry_str in background_list:
				target_no = target_no + 1
		count_obj = count()
		background_term = count_obj.count_terms(background_list,associate_id_map)
		target_term = count_obj.count_terms(target_list,associate_id_map)
		background_no = len(background_list)
		summary = []
		for term, target_count in target_term.items():
			background_count = background_term[term]
			target_other = target_no - target_count
			background_other = background_no - background_count
			oddsratio, pvalue = stats.fisher_exact([[target_count,target_other],
													[background_count,background_other]])
			ratio_target = float(target_count) / float(target_no)
			ratio_background = float(background_count) / float(background_no)
			summary.append([term, associate_map_desc[term],target_count, \
				target_no, ratio_target, background_count, \
				background_no, ratio_background, pvalue])
		pval_sum = [x[8] for x in summary]
		adjust_stats = importr('stats')
		p_adjust_fdr = adjust_stats.p_adjust(FloatVector(pval_sum),method='fdr')
		p_adjust_bf = adjust_stats.p_adjust(FloatVector(pval_sum),method='bonferroni')
		for row_no in range(len(summary)):
			summary[row_no].extend([p_adjust_fdr[row_no],p_adjust_bf[row_no]])
		print("Writing to file")
		with open(pathway_enrich_out,"a") as csvfile:
			writer = csv.writer(csvfile, delimiter="\t")
			writer.writerow(['pathway', 'pathname', 'target number', \
				'total target numbers','ratio of  targets', 'background_number', \
				'total background numbers', 'ratio of background', 'pvalue', \
				'FDR','Bonferroni'])
			for line in summary:
				writer.writerow(line)
		print("# %0.2f seconds process time" % (time.clock() - timer))        

		
class count(object):
	def count_terms(self, geneset, assoc):
		term_cnt = collections.defaultdict(int)
		self.assoc = assoc
		for each_gene in geneset:
			each_gene_str = ''.join(each_gene)
			if each_gene_str in assoc:
				term_cnt[assoc[each_gene_str]] += 1

		return term_cnt






