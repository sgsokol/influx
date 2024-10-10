import sys
from pathlib import Path
from platform import python_version
from packaging.version import Version
if Version(python_version()) < Version("3.5"):
    raise Exception("Python 3.5 or higher required")
#print("sys.path1=", sys.path)
dirmod=Path(__file__).resolve().parent
sys.path.insert(0, str(dirmod))
sys.path.insert(0, str(dirmod/"bin"))
#print("sys.path2=", sys.path)
__version__ = (dirmod/"influx_version.txt").read_text().rstrip()
del(dirmod)
