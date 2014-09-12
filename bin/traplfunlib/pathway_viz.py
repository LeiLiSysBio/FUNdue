from traplfunlib.paths import Paths
import rpy2.robjects as rpy
from subprocess import call
import os

class Pathviz(object):
	"""Uses mapping to detect the go terms"""
	def __init__(self,pathway_script, path_result, target_id, code, output):
		self._path_result = path_result
		self._target_id = target_id
		self._code = code
		self._output = output
		self._pathway_script = pathway_script

	def path_viz(self):
		with open(self._path_result,"r") as Path:
			print(self._path_result)
			next(Path)
			for entry in Path:
				line = entry.rstrip().split("\t")
				if float(line[3]) >= 1 or float(line[8]) >= 1:
					print(self._code)
					print(line[0][3:])
					with open(os.devnull, "w") as devnull:
						call(["Rscript", self._pathway_script, \
						"-o", self._code, "-p", str(line[0][3:]), \
						"-r", self._target_id, "-f", self._output])
			