from traplfunlib.paths import Paths

class Goviz(object):
	"""Uses mapping to detect the go terms"""
	"""For the gene not found in uniprot, using blast2go"""

	def go_viz(self, enrichment_file, viz_revigo):
		viz_revigo_file = open(viz_revigo,"a")
		viz_revigo_file.write("% created by TRAPL_FUN version 0.2" + "\n")
		viz_revigo_file.write("% Enriched gene ontology list" + "\n")
		viz_revigo_file.write("% p-value represent the enrichment of gene ontology terms" + "\n")
		viz_revigo_file.write("% GeneGroup" + "\t" + "pValue" + "\n")
		with open(enrichment_file,"r") as go_enrich:
			next(go_enrich)
			for entry in go_enrich:
				uni_line = entry.rstrip("\n")
				uni_lines = uni_line.split("\t")
				if float(uni_lines[8]) <= 0.1:
					viz_revigo_file.write(uni_lines[0] + "\t" + uni_lines[8] + "\n")
		