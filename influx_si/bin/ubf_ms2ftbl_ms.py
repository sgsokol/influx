#! /usr/bin/python3

"""
convert ubf ms data to ftbl MASS_SPECTROMETRY
"""

import sys
import re

print("MASS_SPECTROMETRY")
print("	META_NAME	FRAGMENT	WEIGHT	VALUE	DEVIATION")

ls=sys.stdin.readlines()
for met, data in zip(ls[::2], ls[1::2]):
   data = data.split(" ")
   met = re.sub("_(C\d)-", "\t\\1-", met)
   m, f = re.match("(_?[a-zA-Z0-9_-]*\w+)(\tC\d-C\d fragment)?", met).group(1,2)
   m = m.strip()
   if f:
      b, e = re.match("\tC(\d)-C(\d) fragment", f).group(1, 2)
   else:
      b, e = "1", str(len(data)-1)
   print(f"\t{m}\t{b}~{e}\t0\t{data[0].strip()}\t0.01")
   for i, d in enumerate(data[1:]):
      print(f"\t\t\t{i+1}\t{d.strip()}\t0.01")
