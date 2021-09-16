import sys
from pathlib import Path
if sys.version_info[0] < 3:
    raise Exception("MPython 3 or higher required")
dirmod=Path(__file__).resolve().parent
sys.path.append(str(dirmod))
__version__ = (dirmod/"influx_version.txt").read_text().rstrip()
del(dirmod)
