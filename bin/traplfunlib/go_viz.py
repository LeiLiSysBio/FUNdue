from traplfunlib.paths import Paths
from traplfunlib.obo_parser import GODag
import os
import matplotlib
matplotlib.use("Agg")
import pylab as pl
import rpy2.robjects as rpy
from subprocess import call
import urllib
import urllib.request
import urllib.parse
from operator import itemgetter
from collections import OrderedDict


class Goviz(object):
	"""Uses mapping to detect the go terms"""
	def __init__(self,obo_path,slim_obo_path,viz_path):
		self._obo_path = obo_path
		self._slim_obo_path = slim_obo_path
		self._viz_path = viz_path
	
	def _mapslim(self, go_term, go_dag, goslim_dag):
		# check parameters
		if not isinstance(go_dag, GODag):
			raise TypeError("go_dag must be an instance of GODag")
		if not isinstance(goslim_dag, GODag):
			raise TypeError("goslim_dag must be an instance of GODag")
		if go_term not in go_dag:
			raise ValueError("go_term must be an accession that is in the go_dag")
		all_ancestors = set()
		covered_ancestors = set()
		# get all paths for the term in the go_dag
		paths = go_dag.paths_to_top(go_term)
		for path in paths:
			# the next loop needs to run bottom->up, i.e. from the go_term item to
			# the root, thus we need to reverse the list prior to iteration
			path.reverse()
			got_leaf = False
			for term in path:
				if term.id in goslim_dag:
					all_ancestors.add(term.name)
					if got_leaf:
						covered_ancestors.add(term.name)
					got_leaf = True
		# get the direct ancestors, i.e. those that are not covered by a earlier
		# ancestor of the GO-Slim in _any_ path (in bottom->top order)
		direct_ancestors = all_ancestors - covered_ancestors
		return direct_ancestors
	
	def _draw_go(self,go_list,go_viz_output):
		url = "http://amigo.geneontology.org/visualize"
		color='lightblue'
		new_go_list=[]
		for go_id in go_list:
			new_go_list.append('\"'+str(go_id)+'\":{\"fill\": \"'+str(color)+'\"}')
		term_data="{"+','.join(str(go_tmp) for go_tmp in new_go_list)+"}"
		term_data_type='json'
		postdata = urllib.parse.urlencode(
						{'mode': 'amigo',
						'inline': 'false',
						'format':'png',
						'term_data_type':'json',
						'term_data':term_data
						})
		postdata = postdata.encode('utf-8')
		outfile=open(go_viz_output,'wb')
		res=urllib.request.urlopen(url,postdata)
		outfile.write(res.read())
		outfile.close()
	
	def go_viz(self, enrichment_file, viz_revigo, viz_tag):
		#For REVIGO use
		viz_revigo_file = open(viz_revigo,"a")
		viz_revigo_file.write("% created by TRAPL_FUN version 0.2" + "\n")
		viz_revigo_file.write("% Enriched gene ontology list" + "\n")
		viz_revigo_file.write("% p-value represent the enrichment of gene ontology terms" + "\n")
		viz_revigo_file.write("% GeneGroup" + "\t" + "pValue" + "\n")
		gene_ontology_object=GODag(self._obo_path)
		slim_gene_ontology_object=GODag(self._slim_obo_path)
		#For slim use
		enrich_go_term_bp_no = {}
		enrich_go_term_bp_id = []
		enrich_go_term_mf_no = {}
		enrich_go_term_mf_id = []
		enrich_go_term_cc_no = {}
		enrich_go_term_cc_id = []
		#For go taxonmy view use
		enrich_go_term_bp = []
		enrich_go_term_mf = []
		enrich_go_term_cc = []
		with open(enrichment_file,"r") as go_enrich:
			next(go_enrich)
			for entry in go_enrich:
				uni_lines = entry.rstrip().split("\t")
				gene_list=uni_lines[10]
				if float(uni_lines[9]) <= 0.05 and uni_lines[1] != "biological_process" \
						and uni_lines[1] != "cellular_component" and \
						uni_lines[1] != "molecular_function":
					viz_revigo_file.write(uni_lines[0] + "\t" + uni_lines[9] + "\n")
					if uni_lines[2] == "biological_process" and float(uni_lines[5]) > float(uni_lines[8]):
						direct_term_set = self._mapslim(uni_lines[0], \
							gene_ontology_object, slim_gene_ontology_object)
						direct_term = ''.join(str(s) for s in direct_term_set)
						enrich_go_term_bp.append(uni_lines[0])
						if direct_term in enrich_go_term_bp_id and direct_term !=  "biological_process":
							sum_no = enrich_go_term_bp_no[direct_term] + int(uni_lines[3])
							enrich_go_term_bp_no[direct_term]= sum_no
						elif direct_term !=  "biological_process":
							enrich_go_term_bp_id.append(direct_term)
							enrich_go_term_bp_no[direct_term]=int(uni_lines[3])
					elif uni_lines[2] == "molecular_function" and float(uni_lines[5]) > float(uni_lines[8]):
						direct_term_set = self._mapslim(uni_lines[0], \
							gene_ontology_object, slim_gene_ontology_object)
						direct_term = ''.join(str(s) for s in direct_term_set)
						enrich_go_term_mf.append(uni_lines[0])
						if direct_term in enrich_go_term_mf_id and direct_term != "molecular_function":
							sum_no = enrich_go_term_mf_no[direct_term] + int(uni_lines[3])
							enrich_go_term_mf_no[direct_term]= sum_no
						elif direct_term != "molecular_function":
							enrich_go_term_mf_id.append(direct_term)
							enrich_go_term_mf_no[direct_term]=int(uni_lines[3])
					elif uni_lines[2] == "cellular_component" and float(uni_lines[5]) > float(uni_lines[8]):
						direct_term_set = self._mapslim(uni_lines[0], \
							gene_ontology_object, slim_gene_ontology_object)
						direct_term = ''.join(str(s) for s in direct_term_set)
						enrich_go_term_cc.append(uni_lines[0])
						if direct_term in enrich_go_term_cc_id and direct_term != "cellular_component":
							sum_no = enrich_go_term_cc_no[direct_term] + int(uni_lines[3])
							enrich_go_term_cc_no[direct_term]= sum_no
						elif direct_term != "cellular_component":
							enrich_go_term_cc_id.append(direct_term)
							enrich_go_term_cc_no[direct_term]=int(uni_lines[3])
		## Draw go taxonomy tree
		self._draw_go(enrich_go_term_bp,viz_tag + "_bp.png")
		self._draw_go(enrich_go_term_mf,viz_tag + "_mf.png")
		self._draw_go(enrich_go_term_cc,viz_tag + "_cc.png")
		## Draw go slim figure
		fig = pl.figure()
		ax = fig.add_subplot(111)
		pl.rcParams['font.size'] = 8.0
		fig.autofmt_xdate()
		enrich_go_term_bp_sort=OrderedDict(sorted(enrich_go_term_bp_no.items(),key=itemgetter(1),reverse=True))
		x_axis=[]
		y_axis=[]
		for line in enrich_go_term_bp_sort.keys():
			x_axis.append(line)
			y_axis.append(int(enrich_go_term_bp_sort[line]))
		pl.bar(range(len(enrich_go_term_bp_sort)), y_axis,color="black")
		enrich_go_term_bp=[]
		pl.xticks(range(len(enrich_go_term_bp_sort)),x_axis)
		pl.title('GO Slim category summary',color='black')
		pl.savefig(viz_tag + "_biological_process.png",format='png',dpi=300,bbox_inches='tight')
		pl.close()
		
		fig = pl.figure(figsize=(12,12))
		ax = fig.add_subplot(111)
		pl.rcParams['font.size'] = 8.0
		fig.autofmt_xdate()
		enrich_go_term_mf_sort=OrderedDict(sorted(enrich_go_term_mf_no.items(),key=itemgetter(1),reverse=True))
		x_axis=[]
		y_axis=[]
		for line in enrich_go_term_mf_sort.keys():
			x_axis.append(line)
			y_axis.append(int(enrich_go_term_mf_sort[line]))
		pl.bar(range(len(enrich_go_term_mf_sort)), y_axis,color="black")
		enrich_go_term_mf=[]
		pl.xticks(range(len(enrich_go_term_mf_sort)),x_axis)
		pl.title('GO Slim category summary',color='black')
		pl.savefig(viz_tag + "_molecular_function.png",format='png',dpi=300,bbox_inches='tight')
		pl.close()
		
		fig = pl.figure(figsize=(1,1))
		ax = fig.add_subplot(111)
		pl.rcParams['font.size'] = 8.0
		fig.autofmt_xdate()
		enrich_go_term_cc_sort=OrderedDict(sorted(enrich_go_term_cc_no.items(),key=itemgetter(1),reverse=True))
		x_axis=[]
		y_axis=[]
		for line in enrich_go_term_cc_sort.keys():
			x_axis.append(line)
			y_axis.append(int(enrich_go_term_cc_sort[line]))
		pl.bar(range(len(enrich_go_term_cc_sort)), y_axis,color="black")
		enrich_go_term_cc=[]
		pl.xticks(range(len(enrich_go_term_cc_sort)),x_axis)
		pl.title('GO Slim category summary',color='black')
		pl.savefig(viz_tag + "_cellular_component.png",format='png',dpi=300,bbox_inches='tight')
		pl.close()
		
		