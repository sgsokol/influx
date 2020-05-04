#! /usr/bin/env python3

"""
usage: ubf2ftbl.py < model.ubf > model.ftbl
"""
# translate University of Barcelona format to plain txt format for labeling flux system
import sys
import re
from string import ascii_lowercase as lett
from datetime import datetime
import argparse
import subprocess

def ubf2txt(lines, fname=""):
   """Read content from lines list till reaction section is found ("// Reactions" case insensitive) then proceed
Return a list of strings formatted as lines of FTBL file
"""
   sres = [] # list strings: name<tab>reac
   res = [] # tuples, 1 per reaction
   metabs = {} # dict like "npyrc": "pyr"
   reacfound = False
   for l in lines:
      l=l.strip()
      # skip empty lines and comments
      if not l or (len(l) > 1 and l[:2] == "//"):
         continue
      # detect reaction start
      if not reacfound:
         if l[:10].lower() == "reactions:":
            reacfound = True
            continue
         # here metabs go
         li = re.split(" +", l)
         if li[0] != "Metab":
            continue
         p,lng,mi,mn = li
         metabs[mi] = mn.strip('" ')
         continue
      # here reactions go
      li = re.split(" +", l)
      if "$" not in li and "->" not in li:
         continue # not reaction string
      idol = li.index("$")
      iarr = li.index("->")
      ## get reaction parts
      nms = li[0]
      nmf = ""
      sep = "->"
      rl = [metabs.get(it.strip("y[]./"), it.strip("y[]./")) for it in li[1:iarr]]
      rr = [metabs.get(it.strip("y[]./"), it.strip("y[]./")) for it in li[iarr+1:idol]]
      if not rl and not rr:
         continue
      al,ar = [a.split("+") for a in li[-1].strip("()").split(">")]
      # translate '123' to 'abc'
      for a in (ar, al):
         a[:] = ["".join(lett[int(l)-1] for l in it) for it in a]
      # complete metab list
      if not rl and rr:
         rl=[s+"_in" for s in rr]
         al=ar
      if not rr and rl:
         rr=[s+"_sink" for s in rl]
         ar=al
      for r,a in [(rl, al), (rr, ar)]:
         if len(r) < len(a) and all(len(ai) == len(a[0]) for ai in a):
            # repeat last item in rl to match the len(al)
            r += [r[-1] for i in range(len(a)-len(r))]
      # complete lacking carbon in/out
      sl="".join(al)
      sr="".join(ar)
      ndiff=len(sl) - len(sr)
      if ndiff < 0:
         # add _cX_in
         rl.append("_c"+str(-ndiff)+"_in")
         al.append("".join(sorted(set(sr)-set(sl))))
      elif ndiff > 0:
         # add _cX_sink
         rr.append("_c"+str(ndiff)+"_sink")
         ar.append("".join(sorted(set(sl)-set(sr))))
      # make 'Mal' symmetric
      for (r,a) in ((rl,al),(rr,ar)):
         try:
            i=r.index("Mal")
         except:
            continue
         #import pdb; pdb.set_trace()
         a[i]=a[i]+"+"+a[i][::-1]
      res.append((nmf, nms, rl, sep, rr, al, ar))

   # format txt as
   # reac_name:<tab>a (laba) + b (labb) -> c (labc)
   sres += [f"# automatically produced "+(f"from '{fname}' " if fname else "")+f"by ubf2ftbl.ubf2txt()\n# at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n"]
   for r in res:
      if r[0]:
         sres += [f"# {r[0]}"] # full reac name in comment
      sres += [f"{r[1]}:\t{' + '.join(m+' ('+l+')' for m,l in zip(r[2], r[5]))} {r[3]} {' + '.join(m+' ('+l+')' for m,l in zip(r[4], r[6]))}"]
   return(sres)
def ubf2ms(lines, fname='', ms1=False):
   """Convert ubf ms data to ftbl MASS_SPECTROMETRY section
Return list of strings
"""
   sres = []
   sres.append("MASS_SPECTROMETRY")
   sres.append("	META_NAME	FRAGMENT	WEIGHT	VALUE	DEVIATION")
   sres.append(f"	// automatically produced "+(f"from '{fname}' " if fname else "")+f"by ubf2ftbl.ubf2ms()\n\t// at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
   onelineformat = len(lines)%2 == 1 or (lines and ":" in lines[0])
   source = (l.split(":") for l in lines if l.strip()) if onelineformat else zip(lines[::2], lines[1::2])
   for met, data in source:
      data = [d.strip() for d in data.strip().split(" ")]
      if ms1:
         # norm data to sum up to 1
         data = [float(d) for d in data]
         s = sum(data)
         data = [d/s for d in data]
      met = re.sub("_?(fragment )?[Cc](\d)-?[Cc](\d)( fragment)?\s*", "\t\\2\t\\3", met)
      li = met.split("\t")
      m, b, e = li if len(li) > 1 else (met, "1", str(len(data)-1))
      m = m.strip()
      sres.append(f"\t{m}\t{b}~{e}\t0\t{data[0]}\t0.01")
      for i, d in enumerate(data[1:]):
         sres.append(f"\t\t\t{i+1}\t{d}\t0.01")
   return(sres)

parser = argparse.ArgumentParser(description="Convert and compile University of Barcelona format (UBF) files to FTBL format")
parser.add_argument("network", help="file name having network reactions and labeling in UBF format")
parser.add_argument("--ms", help="file name having MS measurements in UBF format")
parser.add_argument("--textonly", help="stop treatment at text stage, no ftbl is produced, --ms arguments is ignored", action="store_true")
parser.add_argument("--ms1", help="norm MS data to sum up to 1", action="store_true")
args = parser.parse_args()
# convert ubf to txt
with open(args.network, "r") as f:
   txt = "\n".join(ubf2txt(f.readlines(), args.network))
if args.textonly:
   print(txt)
   exit(0)
# convert txt to FTBL
ftbl = subprocess.run(["txt2ftbl"], input=txt, capture_output=True, text=True)
if ftbl.returncode != 0:
   raise Exception(ftbl.stderr)
ftbl = ftbl.stdout.split("\n")
# convert ms to ftbl
if args.ms:
   with open(args.ms, "r") as f:
      ms = ubf2ms(f.readlines(), args.ms, args.ms1)
   # where ms starts in ftbl?
   ims = ftbl.index("MASS_SPECTROMETRY")
   ftbl = ftbl[:ims+2] + ms[2:] + ftbl[ims+2:]

print("\n".join(ftbl))
