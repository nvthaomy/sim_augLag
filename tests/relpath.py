import os, sys

def setrelpath(file):
    dir_path = os.path.dirname(os.path.realpath(file))
    newpath = "/".join(dir_path.split("/")[:-3])
    sys.path.insert(0,newpath)
    print("expected sim path: {}".format(newpath))
