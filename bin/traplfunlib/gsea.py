from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector
import os
from scipy import stats
import math
import operator
from itertools import islice
from collections import defaultdict
import csv


class GSEA(object):
	"""Perform gene set analysis"""
	def __init__(self,database,genelist,pathout,fdr):
		self._db = database
		self._genelist = genelist
		self._fdr = fdr
		self._pathout = pathout
		
	def _take(self,n, iterable):
		"Return first n items of the iterable as a list"
		return list(islice(iterable, n))	
		
	def gsea_analysis(self):
		path_id = defaultdict(list)
		gene_list = {}
		path_desc = {}
		with open(self._genelist,"r") as gene_file:
			for gene_line in gene_file:
				gene_uni_line = gene_line.rstrip("\n").split("\t")
				if gene_uni_line[0] not in gene_list and gene_uni_line[1] != 'NA' :
					gene_list[gene_uni_line[0]] = float(gene_uni_line[1])
		length = len(gene_list)
		sort_gene_list = sorted(gene_list.items(),key=operator.itemgetter(1),reverse=True)
		
		"""Read database file
		Format: gene_id path_id path_name"""
		with open(self._db,"r") as db_file:
			next(db_file)
			for line in db_file:
				uni_line = line.rstrip("\n").split("\t")
				path_id[uni_line[1]].append(uni_line[0])
				if uni_line[1] not in path_desc:
					path_desc[uni_line[1]] = uni_line[2]	
								
		"""For each pathway, iterature all the ranked gene list"""
		sum = []
		count = 0
		for id_list in path_id:
			count = count + 1
			print("processing term id\t" + str(id_list) + "\t" + str(count) + "/" + str(len(path_id)))
			sort_total_gene_list = self._take(length, sort_gene_list)
			gene_in_this_pathway =len([x for x in sort_total_gene_list if x[0] in path_id[id_list]])
			gene_not_this_pathway = length - gene_in_this_pathway
			max_p = 0
			max_m = 0
			max_p_stat = 0
			max_m_stat = 0
			max_p_g_up = 0
			max_m_g_down = 0
			max_p_g = 0
			max_m_g = 0
			val_list = []
			for i in range(1,length+1):
				sort_n_gene_list = self._take(i, sort_gene_list)
				hit =len([x for x in sort_n_gene_list if x[0] in path_id[id_list]])
				#print(str(i) + "\t" + str(length) + "\t" + str(hit) + "\t" + str(gene_in_this_pathway))
				p_val_up = stats.hypergeom.sf(hit-1,length,gene_in_this_pathway,i)
				p_val_down = stats.hypergeom.sf(gene_in_this_pathway - hit - 1, \
							length, gene_in_this_pathway,length-i)
				if p_val_up == '1.0':
					log_p_up = '0.0'
				else:
					log_p_up = -math.log(p_val_up,10)
				if p_val_down == '1.0':
					log_p_down = '0.0'
				else:
					log_p_down = -math.log(p_val_down,10)
				## up-regulate genes
				if sort_gene_list[i-1][1] >= 0:
					if log_p_up > max_p:
						max_p = log_p_up
						max_p_stat = sort_gene_list[i-1][1]
						max_p_g = hit
						max_p_g_up = i
				elif sort_gene_list[i-1][1] < 0:
					if log_p_down > max_m:
						max_m = log_p_down
						max_m_stat = sort_gene_list[i-1][1]
						max_m_g = gene_in_this_pathway - hit
						max_m_g_down = length - i
				val_list.append([i,log_p_up,log_p_down])
			sort_gene_list_value_up = [float(str(x[1]).replace('-0.0','0.0')) for x in sort_total_gene_list]
			sort_gene_list_value_down = [-float(str(x[1]).replace('-0.0','0.0')) for x in sort_total_gene_list]
			val_list_up = [float(str(x[1]).replace('-0.0','0.0')) for x in val_list]
			val_list_down = [float(str(x[2]).replace('-0.0','0.0')) for x in val_list]
			up_cor, up_val = stats.spearmanr(sort_gene_list_value_up,val_list_up)
			down_cor, down_val = stats.spearmanr(sort_gene_list_value_down,val_list_down)
			sum.append([id_list, path_desc[id_list], gene_in_this_pathway, max_p, \
				up_cor, max_p_stat, max_p_g, max_p_g_up, max_m,down_cor,max_m_stat,max_m_g, \
				max_m_g_down])
		print("sum up the pathway result...")	
		adjust_stats = importr('stats')
		pval_up = [10**(-x[3]) for x in sum]
		pval_down = [10**(-x[8]) for x in sum]
		print("Calculating false discvoery rate")
		p_adjust_up = adjust_stats.p_adjust(FloatVector(pval_up),method='BH')
		p_adjust_down = adjust_stats.p_adjust(FloatVector(pval_down),method='BH')
		for row_no in range(len(sum)):
			sum[row_no][3] = -math.log(p_adjust_up[row_no],10)
			sum[row_no][8] = -math.log(p_adjust_down[row_no],10)
		# write to one file
		print("Writing to file")
		with open(self._pathout,'wt') as csvfile:
			writer = csv.writer(csvfile,delimiter="\t")
			writer.writerow(['pathway','pathname','path_genes', \
				'up_p','up_cor', 'up_max_stat','path_genes_up', \
				'genes_up', 'down_p','down_cor','down_max_stat', \
				'path_genes_down', 'genes_down'])
			for line in sum:
				writer.writerow(line)
					
							
			
	