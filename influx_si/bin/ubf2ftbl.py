#! /usr/bin/env python3

"""
usage: ubf2txt.py < model.ubf > model.txt
"""
# translate University of Barcelona format to plain txt format for labeling flux system
import sys
import re
from string import ascii_lowercase as lett
from datetime import datetime
import argparse
import subprocess

def ubf2txt(lines):
   """Read content from lines list till reaction section is found ("// Reactions" case insensitive) then proceed
Return a list of strings formatted as lines of FTBL file
"""
   sres = [] # list strings: name<tab>reac
   res = [] # tuples, 1 per reaction
   reacfound = False
   for l in lines:
      if not reacfound:
         if l[:12].lower() == "// reactions":
            reacfound = True
         continue
         
      # here reactions go
      ## get names
      nms, parts = l.split(":")
      nmf, nms = re.match("([^,]*,)?(.+)", nms).group(1,2) # get name full, name short
      nmf = (nmf[:-1] if nmf else "").strip()
      nms = nms.strip()
      ## split reac and atom mapping
      sep = re.findall(" *(<?[-=]+>?) *", parts)
      if len(sep) != 2:
         raise Exception(f"Expected two reaction sings, instead got '{sep}' in '{l}'")
      if sep[0] != sep[1]:
         tmp = "' and '".join(sep)
         raise Exception(f"Expected two identical reaction sings, instead got '{tmp}' in '{l}'")
      sep = sep[0]
      try:
         p1, p2, p3 = re.split(" *<?[-=]+>? *", parts)
      except:
         sys.stderr.write(f"line '{l}'\n")
      rl = re.findall("[^+ ]+", p1)
      ar = re.findall("[0-9]+", p3)
      rr, al = zip(*(re.findall("([_a-zA-Z][^+ ]*)|([0-9]+)", p2)))
      rr = list(filter(None, rr))
      al = list(filter(None, al))
      # translate '123' to 'abc'
      for a in (ar, al):
         a[:] = ["".join(lett[int(l)-1] for l in it) for it in a]
      # complete metab list
      for r,a in [(rl, al), (rr, ar)]:
         if len(r) < len(a) and all(len(ai) == len(a[0]) for ai in a):
            # repeat last item in rl to match the len(al)
            r += [r[-1] for i in range(len(a)-len(r))]
      res.append((nmf, nms, rl, sep, rr, al, ar))

   # format txt as
   # reac_name:<tab>a (laba) + b (labb) -> c (labc)
   sres += [f"# automatically produced by ubf2ftbl.ubf2txt()\n# at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n"]
   for r in res:
      if r[0]:
         sres += [f"# {r[0]}"] # full reac name in comment
      sres += [f"{r[1]}:\t{' + '.join(m+' ('+l+')' for m,l in zip(r[2], r[5]))} {r[3]} {' + '.join(m+' ('+l+')' for m,l in zip(r[4], r[6]))}"]
   return(sres)
def ubf2ms(lines):
   """Convert ubf ms data to ftbl MASS_SPECTROMETRY section
Return list of strings
"""
   sres=[]
   sres.append("MASS_SPECTROMETRY")
   sres.append("	META_NAME	FRAGMENT	WEIGHT	VALUE	DEVIATION")

   for met, data in zip(lines[::2], lines[1::2]):
      data = data.split(" ")
      met = re.sub("_(C\d)-", "\t\\1-", met)
      m, f = re.match("(_?[a-zA-Z0-9_-]*\w+)(\tC\d-C\d fragment)?", met).group(1,2)
      m = m.strip()
      if f:
         b, e = re.match("\tC(\d)-C(\d) fragment", f).group(1, 2)
      else:
         b, e = "1", str(len(data)-1)
      sres.append(f"\t{m}\t{b}~{e}\t0\t{data[0].strip()}\t0.01")
      for i, d in enumerate(data[1:]):
         sres.append(f"\t\t\t{i+1}\t{d.strip()}\t0.01")
   return(sres)

parser = argparse.ArgumentParser(description="Convert and compile University of Barcelona format (UBF) files to FTBL format")
parser.add_argument("network", help="file name having network reactions and labeling in UBF format")
parser.add_argument("--ms", help="file name having MS measurements in UBF format")
args = parser.parse_args()
# convert ubf to txt
with open(args.network, "r") as f:
   txt = "\n".join(ubf2txt(f.readlines()))
# convert txt to FTBL
ftbl = subprocess.run(["txt2ftbl"], input=txt, capture_output=True, text=True)
ftbl = ftbl.stdout.split("\n")
# convert ms to ftbl
if args.ms:
   with open(args.ms, "r") as f:
      ms = ubf2ms(f.readlines())
   # where ms starts in ftbl?
   ims = ftbl.index("MASS_SPECTROMETRY")
   ftbl = ftbl[:ims+2] + ms[2:] + ftbl[ims+2:]

print("\n".join(ftbl))
