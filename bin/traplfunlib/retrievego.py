from traplfunlib.gff3 import Gff3Parser
from traplfunlib.paths import Paths

class RetrieveGo(object):
    """Uses mapping to detect the go terms"""
    """For the gene not found in uniprot, using blast2go"""

    def retrieve_go(self, name_list, id_mapping, background_id):
        """ retrieve go terms  """

        idmapping = open(id_mapping,"r")
        background_id_file = open(background_id,"a")
        print("mapping your file, normally need 10-30min to finish")
        num = 0
        for uni_id in idmapping:
            uni_line = uni_id.rstrip("\n")
            uni_lines = uni_line.split("\t")
            #print("You file is mapping..." + uni_lines[3])
            #for ref_name in name_list:
            if uni_lines[3] in name_list:
                num = num + 1
                print("Getting your %s hit" % num )
                background_id_file.write(uni_lines[3] + "\t" + uni_lines[0] + "\t" + \
                                         uni_lines[2] + "\t" + uni_lines[6] + "\n")

        background_id_file.close()

