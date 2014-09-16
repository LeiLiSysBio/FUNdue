from traplfunlib.paths import Paths
import rpy2.robjects as rpy
from subprocess import call
import os
import matplotlib
matplotlib.use("Agg")
import pylab as pl
import numpy as np
import operator
from itertools import islice


class Pathviz(object):
	"""Uses mapping to detect the go terms"""
	def __init__(self,pathway_script, path_result, target_id, code, output):
		self._path_result = path_result
		self._target_id = target_id
		self._code = code
		self._output = output
		self._pathway_script = pathway_script

	def path_viz(self):
		enrich_path_up_list = {}
		enrich_path_up = []
		enrich_path_id_up = []
		enrich_path_down_list = {}
		enrich_path_down = []
		enrich_path_id_down = []
		with open(self._path_result,"r") as Path:
			print(self._path_result)
			next(Path)
			for entry in Path:
				line = entry.rstrip().split("\t")
				if float(line[3]) >= 2:
					enrich_path_up_list[str(line[0])]=float(line[3])
					print("processing Upregulated pathway\t" + str(self._code) + str(line[0][3:]))
					with open(os.devnull, "w") as devnull:
						call(["Rscript", self._pathway_script, \
						"-o", self._code, "-p", str(line[0][3:]), \
						"-r", self._target_id, "-f", self._output])
				elif float(line[8]) >= 2:
					enrich_path_down_list[str(line[0])]=float(line[8])
					print("processing Downregulated pathway\t" + str(self._code) + str(line[0][3:]))
					with open(os.devnull, "w") as devnull:
						call(["Rscript", self._pathway_script, \
						"-o", self._code, "-p", str(line[0][3:]), \
						"-r", self._target_id, "-f", self._output])
		sort_enrich_path_up_list = sorted(enrich_path_up_list.items(),key=operator.itemgetter(1),reverse=False)
		sort_enrich_path_down_list = sorted(enrich_path_down_list.items(),key=operator.itemgetter(1),reverse=False)
		for up_list in sort_enrich_path_up_list:
			enrich_path_id_up.append(up_list[0])
			enrich_path_up.append(up_list[1])
		for down_list in sort_enrich_path_down_list:
			enrich_path_id_down.append(down_list[0])
			enrich_path_down.append(down_list[1])
		target_id_filename = os.path.basename(self._target_id)
		fig = pl.figure()
		width=0.5
		ax = fig.add_subplot(111)
		pl.rcParams['font.size'] = 7.0
		pl.barh(range(len(enrich_path_id_up)),enrich_path_up,color = 'grey')
		ax.set_yticks(np.arange(len(enrich_path_id_up)) + width/2)
		ax.set_yticklabels(enrich_path_id_up)
		pl.title('Upregulated pathways',color='black')
		pl.xlabel("-log10(P)")
		pl.ylabel("Pathway")
		pl.savefig(self._output + "/"+ target_id_filename + "_upregulated_pathway.png",format='png',dpi=400)
		fig.autofmt_xdate()
		pl.close()
		fig = pl.figure()
		width=0.5
		ax = fig.add_subplot(111)
		pl.rcParams['font.size'] = 7.0
		pl.barh(range(len(enrich_path_id_down)),enrich_path_down,color = 'grey')
		ax.set_yticks(np.arange(len(enrich_path_id_down)) + width/2)
		ax.set_yticklabels(enrich_path_id_down)
		pl.xlabel("-log10(P)")
		pl.ylabel("Pathway")
		pl.title('Downregulated pathways',color='black')
		pl.savefig(self._output + "/"+ target_id_filename + "_downregulated.png",format='png',dpi=400)
		fig.autofmt_xdate()
		pl.close()