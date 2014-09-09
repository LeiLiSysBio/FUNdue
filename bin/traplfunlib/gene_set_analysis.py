from scipy import stats
import collections
import tempfile
from traplfunlib.obo_parser import GODag
from subprocess import call
import os

class GSEA_analysis(object):
	"""Perform gene set analysis"""
	def __init__(self,pvalue,foldchange,method,minsize,minedge):
		self._pvalue = pvalue
		self._foldchange = foldchange
		self._method = method
		self._minsize = minsize
		self._minedge = minedge
	
	def gsea(self, gsea_script, expression_file,go_background_file, obo_file_path, prefix):
		"""Expression files format gene_id log2foldchange pval"""
		go_obo = GODag(obo_file_path)
		tmp_gofile = tempfile.NamedTemporaryFile(mode="a",delete=False)
		print("Processing the gene ontology file, Please wait...")
		for entry in open(go_background_file, "r"):
			uni_line = entry.rstrip().split("\t")
			if len(uni_line) > 3:
				if ';' in uni_line[3]:
					go_list = uni_line[3].replace(" ", "").split(";")
					for each_go_list in go_list:
						tmp_gofile.write(uni_line[0] + "\t" + \
							go_obo[each_go_list].name + "\n")
				else:
					tmp_gofile.write(uni_line[0] + "\t" + \
						go_obo[uni_line[3]].name + "\n")
				
		prefix_go = prefix + "_go"
		print("Running Gene Set Analysis, Please wait...")
		with open(os.devnull, "w") as devnull:
			call(["Rscript", gsea_script, \
			"--tsv" , expression_file,"-p",str(self._pvalue),"-f",str(self._foldchange), \
			"-s",str(self._minsize),"-e", str(self._minedge), "-m", str(self._method), \
			"--go", tmp_gofile.name, "-o", prefix_go])
		tmp_gofile.close()
		
