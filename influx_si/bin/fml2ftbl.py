#! /usr/bin/env python3

import xmltodict
from datetime import datetime
import re

def is_list(l): return(type(l)==type([]))

def is_num(x):
   try:
      float(x)
      return(True)
   except ValueError:
      return(False)

# inverse inequality sign
def revineq(op): return("<=" if op == ">=" else ">=")

# get default value from list
def get_dl(obj, i, dflt): return(obj[i] if i < len(obj) else dflt)

# get xml to dict
with open('Ecoli3.fml') as fd:
   doc = xmltodict.parse(fd.read())

# check that we are in statioanry mode
if doc['fluxml']['configuration']['@stationary'].lower() != "true":
   raise Exception("instationary mode is not treated yet")

# extract reactions
reacs={}
rout=set()
for rxml in doc['fluxml']['reactionnetwork']['reaction']:
   reac=rxml['@id']
   key='reduct'
   l=[(p['@id'], '#'+p['@cfg']) for p in (rxml[key] if is_list(rxml[key]) else [rxml[key]])] if key in rxml else []
   key='rproduct'
   r=[(p['@id'], '#'+p['@cfg']) for p in (rxml[key] if is_list(rxml[key]) else [rxml[key]])] if key in rxml else []
   if not r:
      # if nothing on rhs, put inputs with prefix "_sink_" and add reac to rout set
      key='reduct'
      r=[("_sink_"+p['@id'], '#'+p['@cfg']) for p in (rxml[key] if is_list(rxml[key]) else [rxml[key]])] if key in rxml else []
      rout.add(reac)
   reacs[reac]=[l, r] # id: (list of 2-tuples educt, list of 2-tuples product)

# extract constraints
cnx={"net": {"c": {}, "eq": {}, "ineq": {}}, "xch": {"c": {}, "eq": {}, "ineq": {}}}
pcnx=re.compile("([^<>=]*)([<>=]+)([^<>=]*)")
for nx in ("net", "xch"):
   s=doc['fluxml']['constraints'][nx]['textual']
   s=[[t.strip() for t in pcnx.match(ex).group(1,2,3)] for ex in s.split(";")]
   for (l,op,r) in s:
      if op == "=":
         if is_num(l) or is_num(r):
            v,e=(l,r) if is_num(l) else (r,l)
            if e in reacs:
               # constrained flux
               cnx[nx]["c"][e]=v
            else:
               # eq constraint
               cnx[nx]["eq"][e]=v
         else:
            # pass all to rhs=0
            cnx[nx]["eq"][r+" - "+l.replace("-", "{:}").replace("+", "-").replace("{:}", "+")]="0"
      elif op == "<=" or op == ">=":
         if is_num(l) or is_num(r):
            # value at left leaves op as is
            v,e,o=(l,r,op) if is_num(l) else (r,l,revineq(op))
            cnx[nx]["ineq"][e]=(v, o)
         else:
            cnx[nx]["ineq"][r+" - "+l.replace("-", "{:}").replace("+", "-").replace("{:}", "+")]=("0", op)
for r in rout:
   if r not in cnx["xch"]["c"]:
      cnx["xch"]["c"][r]="0"

# extract flux starting values
fst={'net': {}, 'xch': {}}
for r in doc['fluxml']['configuration']['simulation']['variables']['fluxvalue']:
   nx=r['@type']
   fl=r['@flux']
   if fl in cnx[nx]['c']:
      raise Exception(f"{nx} flux '{fl}' is both in constrained (fluxml/constraints/{nx}) and free fluxes (fluxml/configuration/simulation/variables/fluxvalue)")
   fst[nx][fl]=r['#text']

# extract measurements
## data
mdoc=doc['fluxml']['configuration']['measurement']['data']['datum']
mdoc=[mdoc] if  not is_list(mdoc) else mdoc
mdata={}
for d in mdoc:
   mdata[d['@id']]=mdata.get(d['@id'], [])
   mdata[d['@id']].append({'val': d['#text'], 'sd': d['@stddev'], 'ms+': d.get('@weight')})

## flux
meas={"f": {}}
l=doc['fluxml']['configuration']['measurement']['model']['fluxmeasurement']['netflux']
l=[l] if not is_list(l) else l
for fm in l:
   meas["f"][fm['textual']]=mdata[fm['@id']][0]

## ms
meas['ms']=[]
l=doc['fluxml']['configuration']['measurement']['model']['labelingmeasurement']['group']
l=[l] if not is_list(l) else l
for gr in l:
   meas['ms'].append((gr['textual'], gr['@scale'], mdata[gr['@id']]))

# extract input labeling
ldoc=doc['fluxml']['configuration']['input']
ldoc=[ldoc] if not is_list(ldoc) else ldoc
inlab={}
for dlab in ldoc:
   if dlab['@type'].lower() != 'isotopomer':
      raise Exception(f"for pool '{dlab['@pool']}' type of input labeling is not isotopomer")
   inlab[dlab['@pool']]=[(it['@cfg'], it['#text']) for it in dlab['label']]

#-----------------------------------------------------
# make output
# print FTBL header
print(f"""PROJECT
	NAME	VERSION	FORMAT	DATE	COMMENT
	fml2ftbl	1		{datetime.now().strftime("%Y-%m-%d %H:%M:%S")}	converted from fml by fml2ftbl.py

NETWORK
	FLUX_NAME	EDUCT_1	EDUCT_2	PRODUCT_1	PRODUCT_2""")
# print reactions
for rnm in reacs:
   rt=reacs[rnm]
   me,le=zip(*rt[0]) if rt[0] else ([], [])
   mp,lp=zip(*rt[1])if rt[1] else ([], [])
   for i in range((max(len(me), len(mp))+1)//2):
      print(f"\t{rnm}\t{get_dl(me, i, '')}\t{get_dl(me, i+1, '')}\t{get_dl(mp, i, '')}\t{get_dl(mp, i+1, '')}")
      print(f"\t\t{get_dl(le, i, '')}\t{get_dl(le, i+1, '')}\t{get_dl(lp, i, '')}\t{get_dl(lp, i+1, '')}")

# print flux definitions
print("\nFLUXES")
for nx,NX in (('net', 'NET'), ('xch', 'XCH')):
   print(f"""	{NX}
		NAME	FCD	VALUE(F/C)	ED_WEIGHT	LOW(F)	INC(F)	UP(F)""")
   for rnm in reacs:
      if rnm in fst[nx]:
         val=fst[nx][rnm]
         if nx == 'xch':
            val=float(val)
            val=val/(1.+val)
         print(f"\t\t{rnm}\tF\t{val}")
      elif rnm in cnx[nx]['c']:
         val=cnx[nx]['c'][rnm]
         if nx == 'xch':
            val=float(val)
            val=val/(1.+val)
         print(f"\t\t{rnm}\tC\t{val}")
      else:
         print(f"\t\t{rnm}\tD")

# print equalities
print("\nEQUALITIES")
for nx in ('net', 'xch'):
   print("\t"+nx.upper()+"\n\t\tVALUE\tFORMULA")
   for e in cnx[nx]["eq"]:
      v=cnx[nx]["eq"][e]
      print(f"\t\t{v}\t{e}")

# print inequalities
print("\nINEQUALITIES")
for nx in ('net', 'xch'):
   print("\t"+nx.upper()+"\n\t\tVALUE\tCOMP\tFORMULA")
   for e in cnx[nx]["ineq"]:
      v, op=cnx[nx]["ineq"][e]
      print(f"\t\t{v}\t{op}\t{e}")

# print flux measurements
print("\nFLUX_MEASUREMENTS\n\tFLUX_NAME\tVALUE\tDEVIATION")
for rnm,m in meas["f"].items():
   print(f"\t{rnm}\t{m['val']}\t{m['sd']}")

# print label input
print("\nLABEL_INPUT\n\tMETA_NAME\tISOTOPOMER\tVALUE")
for met,li in inlab.items():
   print(f"\t{met}\t#{li[0][0]}\t{li[0][1]}")
   for lit in li[1:]:
      print(f"\t\t#{lit[0]}\t{lit[1]}")

# print ms measurements
print("\nMASS_SPECTROMETRY\n\tMETA_NAME\tFRAGMENT\tWEIGHT\tVALUE\tDEVIATION")
for frag,scale,msli in meas["ms"]:
   met,ali=frag.split("[") # metab and atom list
   ali=ali.split("]")[0]
   if scale == 'auto':
      vsum=0.
      for v in msli:
         v['val']=float(v['val'])
         vsum += v['val']
      for v in msli:
         v['val'] /= vsum
   print(f"\t{met}\t{ali}\t{msli[0]['ms+']}\t{msli[0]['val']}\t{msli[0]['sd']}")
   for it in msli[1:]:
      print(f"\t\t\t{it['ms+']}\t{it['val']}\t{it['sd']}")

# print options
print("""
OPTIONS
	OPT_NAME	OPT_VALUE
	//optctrl_history	1 // nlsic
	commandArgs	--ln --TIMEIT --noscale
	posttreat_R	save_all.R; plot_smeas.R
""")
