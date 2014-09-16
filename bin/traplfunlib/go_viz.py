from traplfunlib.paths import Paths
from traplfunlib.obo_parser import GODag
import os
import matplotlib
matplotlib.use("Agg")
import pylab as pl

class Goviz(object):
	"""Uses mapping to detect the go terms"""

	def go_viz(self, enrichment_file, viz_revigo, viz_tag,obo_path):
		viz_revigo_file = open(viz_revigo,"a")
		viz_revigo_file.write("% created by TRAPL_FUN version 0.2" + "\n")
		viz_revigo_file.write("% Enriched gene ontology list" + "\n")
		viz_revigo_file.write("% p-value represent the enrichment of gene ontology terms" + "\n")
		viz_revigo_file.write("% GeneGroup" + "\t" + "pValue" + "\n")
		gene_ontology_object=GODag(obo_path)
		num = 0
		enrich_go_term_bp = []
		enrich_go_term_bp_id = []
		enrich_go_term_mf = []
		enrich_go_term_mf_id = []
		enrich_go_term_cc = []
		enrich_go_term_cc_id = []
		with open(enrichment_file,"r") as go_enrich:
			next(go_enrich)
			for entry in go_enrich:
				num = num + 1
				uni_line = entry.rstrip("\n")
				uni_lines = uni_line.split("\t")
				if float(uni_lines[9]) <= 0.05 and uni_lines[1] != "biological_process" \
						and uni_lines[1] != "cellular_component" and \
						uni_lines[1] != "molecular_function":
					viz_revigo_file.write(uni_lines[0] + "\t" + uni_lines[9] + "\n")
					if uni_lines[2] == "biological_process" and float(uni_lines[5]) > float(uni_lines[8]):
						enrich_go_term_bp.append(uni_lines[3])
						enrich_go_term_bp_id.append(uni_lines[1])
					elif uni_lines[2] == "molecular_function" and float(uni_lines[5]) > float(uni_lines[8]):
						enrich_go_term_mf.append(uni_lines[3])
						enrich_go_term_mf_id.append(uni_lines[1])
					elif uni_lines[2] == "cellular_component" and float(uni_lines[5]) > float(uni_lines[8]):
						enrich_go_term_cc.append(uni_lines[3])
						enrich_go_term_cc_id.append(uni_lines[1])
					"""Test function, pygraphviz need be installed in the right python version"""
					#print("output for visulization")
					#term = gene_ontology_object.query_term(uni_lines[0],verbose=True)
					#print(viz_tag)
					#print(term)
					#output = viz_tag+'_'+str(num)+'png'
					#gene_ontology_object.draw_lineage([term],lineage_img=output)
		fig = pl.figure(figsize=[12,12])
		ax = fig.add_subplot(111)
		pl.rcParams['font.size'] = 8.0
		pl.pie(enrich_go_term_bp,labels=enrich_go_term_bp_id,autopct='%1.f%%',startangle=90)
		pl.title('Enriched Biological process',color='black')
		pl.savefig(viz_tag + "_biological_process.png",format='png',dpi=400)
		pl.close()
		fig = pl.figure(figsize=[12,12])
		ax = fig.add_subplot(111)
		pl.rcParams['font.size'] = 8.0
		pl.pie(enrich_go_term_mf,labels=enrich_go_term_mf_id,autopct='%1.f%%',startangle=90)
		pl.title('Enriched Molecular function',color='black')
		pl.savefig(viz_tag + "_molecular_function.png",format='png',dpi=400)
		pl.close()
		fig = pl.figure(figsize=[12,12])
		ax = fig.add_subplot(111)
		pl.rcParams['font.size'] = 8.0
		pl.pie(enrich_go_term_cc,labels=enrich_go_term_cc_id,autopct='%1.f%%',startangle=90)
		pl.savefig(viz_tag + "_cellular_component.png",format='png',dpi=400)
		pl.title('Enriched Cellular component',color='#afeeee')
		pl.close()