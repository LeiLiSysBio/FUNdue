import collections
import string
import sys
import argparse
import os
import subprocess
import tempfile

class blast2go_analysis(object):
	
	def run_blast2go(self, blast2go_path, xml_file, background_file):
		
		# transform the xml files
		in_handle = open(xml_file)
		footer = "  </BlastOutput_iterations>\n</BlastOutput>\n"
		header = ""
		while True:
			line = in_handle.readline()
			if not line:
				print("Problem with XML file?")
				sys.exit()
			if line.strip() == "<Iteration>":
				break
			header += line
		
		if "<BlastOutput_program>blastx</BlastOutput_program>" in header:
			print("BLASTX xml detected")
		elif "<BlastOutput_program>blastp</BlastOutput_program>" in header:
			print("BLASTP xml detected")
		else:
			in_handle.close()
			print("Expect BLASTP or BLASTX output")
			sys.exit()
		
		out_handle = tempfile.NamedTemporaryFile(delete=False)
		out_handle.write(header.encode())
		out_handle.write(line.encode())
		count = 1
		while True:
			line = in_handle.readline()
			if not line:
				break
			elif line.strip() == "<Iteration>":
			   #Insert footer/header
			   out_handle.write(footer.encode())
			   out_handle.write(header.encode())
			   count += 1
			out_handle.write(line.encode())
		
		
		in_handle.close()
		print("Input has %i queries" % count)

		# run commands
		cmd = ["java", "-cp", blast2go_path + '/blast2go.jar' + ":" + \
		blast2go_path + "/ext/*", "es.blast2go.prog.B2GAnnotPipe",
		"-in", out_handle.name,
		"-prop", blast2go_path + '/b2gPipe.properties',
		"-out", background_file, #Used as base name for output files
		"-annot", # Generate *.annot tabular file
		"-img", # Generate images, feature not in v2.3.5
		]
		
		try:
			child = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		except Exception:
			print("Error invoking command:\n%s\n\n%s\n" % (" ".join(cmd), err))
			sys.exit()

		stdout, stderr = child.communicate()
		#return_code = child.returncode
		
		if stdout:
			#sys.stderr.write("Standard out:\n%s" % stdout)
			sys.stderr.write("Blast2go finished")
		if stderr:
			sys.stderr.write("Standard error:\n%s\n\n" % stderr)

		#out_handle.close()
		