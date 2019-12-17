import os, sys
if sys.version_info[0] < 3:
    raise Exception("Must be using Python 3 or higher")
dirmod=os.path.dirname(os.path.realpath(__file__))
sys.path.append(dirmod)
with open(os.path.join(dirmod, "influx_version.txt"), "r") as f:
    __version__ = f.read().rstrip()
