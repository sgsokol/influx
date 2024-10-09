import sys
from pathlib import Path
if sys.version_info[0] < 3:
    raise Exception("Python 3 or higher required")
#print("sys.path1=", sys.path)
dirmod=Path(__file__).resolve().parent
sys.path.insert(0, str(dirmod))
sys.path.insert(0, str(dirmod/"bin"))
#print("sys.path2=", sys.path)
__version__ = (dirmod/"influx_version.txt").read_text().rstrip()
del(dirmod)
