from traplfunlib.gff3 import Gff3Parser
from traplfunlib.paths import Paths
from subprocess import call
import os

class Retrievepathway(object):
	## First build a pathway db
	# calling R in python
	# Rscript buildpathwayDB.R -o sey -d sey.txt
	def buildpathway(self, pathway_script, organism_code, background_id):
		with open(os.devnull, "w") as devnull:
			call(["Rscript", pathway_script, \
			"-o" , organism_code, "-d", background_id])
			
	 