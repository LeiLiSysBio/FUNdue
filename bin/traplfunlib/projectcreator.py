import json
import os
import sys

class ProjectCreator(object):

    def create_root_folder(self, project_name):
        if not os.path.exists(project_name):
            os.mkdir(project_name)
        else:
            sys.stderr.write("Cannot create folder \"%s\"! File/folder with "
                             "the same name exists already.\n" % project_name)
            sys.exit(2)

    def create_subfolders(self, subfolders):
        for folder in subfolders:
            if not os.path.exists(folder):
                os.mkdir(folder)

    def create_version_file(self, version_file_path, version):
        with open(version_file_path,"w") as fh:
            fh.write("TRAPL_FUN version %s" % version)
