#!/usr/bin/env python3

"""
transform a series of TXT and TSV files into FTBL file.

Copyright 2021, INRAE, INSA, CNRS
Author: Serguei Sokol (sokol at insa-toulouse dot fr)
License: Gnu Public License (GPL) v2 http://www.gnu.org/licenses/gpl.html
"""

import re
import os
import sys
from pathlib import Path
import pandas as pa
import numpy as np
import stat
import argparse
import datetime as dt
from scipy import linalg
from numpy import diag
vsadd=np.core.defchararray.add # vector string add
import influx_si
from C13_ftbl import formula2dict

version="1.0"
#me=os.path.basename(sys.argv[0] or "txt2ftbl")
me="txt2ftbl"
LOCAL_TIMEZONE=dt.datetime.now(dt.timezone.utc).astimezone().tzinfo
invcomp={">=": "<=", "=>": "<=", "<=": ">=", "=<": ">="}

def natural_sort_key(s, _re=re.compile(r'(\d+)')):
    # last 2 fields are inverted for sorting
    return [int(t) if i & 1 else t.lower() for li in (_re.split(s)[:-1],) for i, t in enumerate(li[:-3]+li[-1:]+li[-2:-4:-1])]
def plain_natural_key(s, _re=re.compile(r'(\d+)')):
    return [int(t) if i & 1 else t.lower() for li in (_re.split(s)[:-1],) for i, t in enumerate(li)]
def dtstamp():
    "formatted date-time stamp"
    return dt.datetime.now(LOCAL_TIMEZONE).strftime('%Y-%m-%d %H:%M:%S %Z %z')
def werr(mes):
    raise Exception(f"{me}: {mes}")
def warn(mes):
    sys.stderr.write(f"Warning! {me}: {mes}\n")
def usage():
    sys.stderr.write(__doc__+"\n")
def itvl2li(v):
    "convert interval like '2-5'  to ['2', '3', '4', '5']"
    li=v.split("-")
    if len(li) == 1:
        return [str(int(li[0]))]
    else:
        return [str(i) for i in range(int(li[0]), int(li[1])+1)]
def dsec2out(dsec, fout):
    "write lines from dsec to fout"
    for k,li in dsec.items():
        for item in li:
            if type(item) == dict:
                dsec2out(item, fout)
            else:
                fout.write(item+"\n")
def tsv2df(f, sep="\t", comment="#", skip_blank_lines=True, append_iline="iline"):
    "Read file 'f' as TSV and return a DataFrame. Separator is 'sep', comment char is 'comment', blank lines are skept, header is in the first row, file line numbers are stored in a column 'line_nb' if there is no a column with this name"
    if type(f) == str or type(f) == type(Path()):
        li=Path(f).open().read().splitlines()
    else:
        li=f.readlines
    rows=np.array([], dtype=object)
    irow=0
    cnm=[] # col names
    for ili,row in enumerate(li):
        # remove comments
        if comment:
            row=re.sub("%s.*$"%comment, "", row)
        if skip_blank_lines and not row.strip():
            continue
        if not cnm:
            cnm=[v.strip() for v in row.split(sep)]
            ncol=len(cnm)
            if append_iline:
                if not any(v==append_iline for v in cnm):
                    cnm.append(append_iline)
                else:
                    append_iline=False
            rows=np.empty((0, len(cnm)), dtype=object)
            continue
        # append rows
        rli=[v.strip() for v in row.split(sep)]
        if len(rli) != ncol:
            werr("tsv2df: wrong column number %d, expected %d (%s: %d)"%(len(rli), ncol, f, ili+1))
        if append_iline:
            rli.append(ili+1)
        rows=np.vstack((rows, rli))
    df=pa.DataFrame(rows)
    df.columns=cnm
    return df
def try_ext(f, li):
    """See if file 'f' exists, if not try with extensions from 'li'.
    The first found is returned as Path() otherwise an exception is raised."""
    if type(f) in (type(Path()), str):
        # check if we need to add an extension
        pth=Path(f)
        if not pth.is_file():
            for suf in li:
                if not suf.startswith("."):
                    suf="."+suf
                npth=pth.with_suffix(suf)
                if npth.is_file():
                    return npth
                npth=pth.parent/(pth.name+suf)
                if npth.is_file():
                    return npth
            raise Exception(f"try_ext: not found file '{f}' neither literal nor with extensions: '"+"', '".join(li)+"'")
        else:
            return pth
    else:
        raise Exception("try_ext: unknown type of 'f'. Expecting 'str' or 'Path'.")
    

def txt_parse(ftxt, re_metab=re.compile(r"(?:(?P<coef>[\d\.]*)\s+)?(?:(?P<metab>[^() \t\r]+)\s*)(?:\(\s*(?P<carb>[^()]*)\s*\))?\s*"),
        re_labpat=re.compile(r"^[./*\d\s]*(?P<labpat>[a-zA-Z]*)\s*$")):
    """Parse txt file from fname which is in format:
### Glycolysis and OPP pathway
GLYC (abcdef) ->  G6P (abcdef)
G6P (abcdef)  <-> F6P (abcdef)
i.e.

comment: # blabla
non reversible reaction: [reac: ][N1] metab1 [(carb1)] [... + [N_i] metab_i [(carb_i)]] -> ...
where
 reac is an optional reaction name;
 N_i is optional stoechiometric coefficient
 metab_i is i-th metabolite name
 (carb_i) is optional 1-letter carbon names for carbon transition mapping

reversible reaction is represented by "<->" sign
non reversible reaction is represented by "->" sign.
reaction with imposed sens of reaction (from left to right) is representd with double ">>",
i.e. "->>" for non reversible reaction or "<->>" for reversible reaction. This
reaction with imposed sens will have an inequality "reac_name >= 0" in
FTBL/INEQUAITIES/NET section.

Retrun a list with following items:
- a list of carbon exchange reactions
- a list of non carbon echanging reactions
- a list of equalities net and xch
- a list of flux tuples (reac, rev, imposed_sens)
- a list of two lists: left and right metabolites [(metab, clen),]
- a carbon length dictionary {met: N}

the first item is a list of:
plain string == just a comments
list == reaction items: input, output: lists of tuples (metab, carb, coeff)
"""
    res=[] # main result. Each item can be a plain str (== comment) or a 4-tuple: reac, reversible?, imposed_sens?, list of in-metabo-tuples, list of out-metabo-tuples
    # each metabo-tuple is (coef, metab, labeling-pattern)
    
    eqs=[[], []] # lists (net & xch) of equalities: tuple (value, formula)
    ineqs=[[], []] # lists (net & xch) of inequalities: tuple (value, comp_sign, formula)
    fluxes=[] # list of tuples (nm_reac, rev)
    resnotr=[] # list of reaction without tracer (like biomass)
    dclen={}
    comment=""
    
    open_here=False
    if type(ftxt) == type(Path()):
        open_here=True
        fc=ftxt.open("r")
        fname=ftxt.name
    else:
        fc=ftxt
        fname=Path(ftxt.name).name
    m_left={} # metab sources
    m_right={} # metab products
    sto={} # stoechiometric matrix dictionary {flux:{metab: coef}}}
    mconn=[] # list of connected metabolite sets. Should be only one
    ireac=0
    ipath=1
    iline=0
    for l in fc.readlines():
        iline=iline+1
        l=l.strip()
        if len(l) == 0:
            continue
        if l[0]=="#":
            comment+="//"+l[1:]+"\n"
            if len(l) > 2 and l[:3] == "###":
                ipath=ipath+(ireac != 0)
                ireac=0
            continue
        # parse reaction
        ireac=ireac+1
        # strip possible end-line comment
        li=l.split("#", 1)
        if len(li)==2:
            l=li[0].strip()
        li=l.split(":", 1)
        if len(li)==2:
            nm_reac=li[0].strip()
            reac=li[1].strip()
        else:
            nm_reac="r"+str(ipath)+"."+str(ireac)
            reac=li[0].strip()
        in_out=[side.strip() for side in reac.split("->")]
        if len(in_out) != 2:
            werr("txt_parse: bad syntax. No reaction detected in '%s': %d"%(fname, iline))
        # True|False is for reversible or not, imposed_sens or not
        rev=False
        imposed_sens=False
        if in_out[1][0] == ">":
            in_out[1]=in_out[1][1:].strip()
            imposed_sens=True
        if in_out[0][-1:] == "<":
            in_out[0]=in_out[0][:-1].strip()
            rev=True
        # parse metabo-tuples
        # we cannot just split by "+" because of possible scrambling patterns having "+" in them
        lr=[[(c, m, [elab for lab in t.split("+") for elab in re_labpat.findall(lab)]) for (c, m, t) in re_metab.findall(io) if m != "+"] for io in in_out]
        tr_reac=len([1 for side in lr for (c, m, t) in side if t and t[0]]) # it exists carbon transition
        #import pdb; pdb.set_trace();
        if tr_reac:
            # gather carbon lengths
            for side in lr:
                for c,m,t in side:
                    if m in dclen:
                        if not dclen[m] == len(t[0]):
                            werr("txt_parse: metabolite '%s' has different label lengths %d and %d. The last is in '%s': %d"%(m, dclen[m], len(t[0]), fname, iline))
                    else:
                        dclen[m]=len(t[0])
                        
            [[[(c,m,t[0]) ] ]]
        if any(len(t) == 0 for side in lr for c,m,t in side):
            werr("txt_parse: wrong format for labeling pattern in '%s': %d."%(fname, iline))
        if not tr_reac:
            if comment:
                resnotr.append(comment)
                comment=""
            resnotr.append([(nm_reac, rev, imposed_sens)]+lr)
            fluxes.append((nm_reac, rev, imposed_sens, "F"))
            # add a row to stoechiometric matrix
            if nm_reac in sto:
                werr("txt_parse: reaction name '%s' was already used ('%s': %d)"%(nm_reac, fname, iline))
            d={}
            d.update((m, -(float(c) if c else 1.)+d.get(m, 0.)) for c,m,t in lr[0])
            d.update((m, (float(c) if c else 1.)+d.get(m, 0.)) for c,m,t in lr[1])
            sto[nm_reac]=d
            # get m_left and m_right metab dicts
            es=dict((m, 0) for c,m,t in lr[0])
            ps=dict((m, 0) for c,m,t in lr[1])
            if rev:
                es.update(ps)
                ps.update(es)
            m_left.update((m,clen) for m, clen in es.items() if m not in m_left)
            m_right.update((m,clen) for m, clen in ps.items() if m not in m_right)
            continue
        
        # get position and length of scrumbles then make twin reactions 1,2,... with equalities
        i_n_scr=[[(i,len(t)) for i,(c,m,t) in enumerate(side)] for side in lr]
        # prepare nested loops, one per scrumble
        ilr=[[0]*len(it) for it in i_n_scr]
        nlr=[[n for i,n in it] for it in i_n_scr]
        lrs=[]
        loop=["for ilr[%(ilr)d][%(isc)d] in range(%(nsc)d):"%{"ilr": il, "isc": isc, "nsc": nsc} for il in range(2) for isc,nsc in enumerate(nlr[il]) if nsc > 1]
        code="\n".join(" "*i+it for i,it in enumerate(loop))
        if loop:
            ire=0
            code+="\n"+" "*len(loop)+"lrs.append([[(c,m,(t[0] if (len(t) == 1) else t[ilr[isi][i]])) for (i, (c,m,t)) in enumerate(side)] for isi, side in enumerate(lr)])"
            exec(code, locals())
        else:
            lrs=[[[(c,m,t[0]) for c,m,t in side] for side in lr]]
        for ilr, lr in enumerate(lrs):
            nm_r=nm_reac+("_"+str(ilr+1) if len(lrs) > 1 else "")
            # add a row to stoechimetric matrix
            if nm_r in sto:
                werr("txt_parse: reaction name '%s' was already used ('%s': %d)"%(nm_r, fname, iline))
            d={}
            d.update((m, -(float(c) if c else 1.)+d.get(m, 0.)) for c,m,t in lr[0])
            d.update((m, (float(c) if c else 1.)+d.get(m, 0.)) for c,m,t in lr[1])
            sto[nm_r]=d
            # get m_left and m_right metab dicts
            es=dict((m, len(t)) for c,m,t in lr[0])
            ps=dict((m, len(t)) for c,m,t in lr[1])
            if rev:
                es.update(ps)
                ps.update(es)
            # update connected metabolite sets
            ise=[i for i,s in enumerate(mconn) if any(m in s for m in es)]
            if not ise:
                mconn.append(set())
                ise=len(mconn)-1
            else:
                ise=ise[0]
            mconn[ise].update(es.keys())
            mconn[ise].update(ps.keys())
            m_left.update((m,clen) for m, clen in es.items() if m not in m_left)
            m_right.update((m,clen) for m, clen in ps.items() if m not in m_right)
            rgr=[] # reaction group (for possible long reactions)
            r=[(nm_r, rev, imposed_sens), [], []] # place for elementary reaction (no more than 2 metabs on each side)
            fluxes.append((nm_r, rev, imposed_sens, "D" if ilr else "F"))
            if ilr:
                teq=(0, "%s - %s"%(nm_reac+"_1", nm_r))
                eqs[0].append(teq)
                if rev:
                    eqs[1].append(teq)
            if ilr == 0 and imposed_sens:
                ineqs[0].append((0, "<=", nm_r))
            # add sto row
            while len(lr[0]) or len(lr[1]):
                # split this reaction in many: max two metabolites on each side
                if len(lr[0]):
                    r[1].append(lr[0].pop(0))
                else:
                    r[1].append(("", "", ""))
                if len(lr[1]):
                    r[2].append(lr[1].pop(0))
                else:
                    r[2].append(("", "", ""))
                if (len(r[1]) == 2 or len(r[2]) == 2) and (len(lr[0]) or len(lr[1])):
                    # we have to split here
                    rgr.append(r[:])
                    r[1]=[]
                    r[2]=[]
            rgr.append(r[:]) # final flush
            if comment:
                res.append(comment)
                comment=""
            res+=rgr
        
    if open_here:
        fc.close()
    # collapse rightmost connected metab set
    rmc=range(len(mconn)-1, 0, -1)
    for i in rmc:
        for j in range(0, i):
            if mconn[j] & mconn[i]:
                mconn[j].update(mconn[i])
                del(mconn[i])
                break
            
    if len(mconn) > 1:
        mconn=sorted(mconn, key=lambda s: len(s), reverse=True)
        werr("txt_parse: detected %d disconnected sub-networks in '%s':\n\t%s"%(len(mconn)-1, fname, "\n\t".join(str(i+1)+": "+", ".join(s) for i,s in enumerate(mconn[1:]))))
    #import pdb; pdb.set_trace()
    return [res, resnotr, eqs, ineqs, fluxes, [m_left, m_right], sto, dclen]
def parse_miso(fmiso, clen, case_i=False):
    "Parse isotopic measurements TSV file. Return dict with keys: ms, lab, peak"
    if "name" in dir(fmiso):
        fname=fmiso.name
    else:
        fname=fmiso
    fname=os.path.basename(fname)
#    import pdb; pdb.set_trace()
    df=tsv2df(fmiso)
    if "Metabolite" in df.columns and "Species" in df.columns:
        df.rename({"Metabolite": "Specie", "Species": "Isospecies"}, inplace=True, axis=1)
    if case_i:
        if "Time" not in df or sum(df["Time"] != "") == 0:
            warn("parse_miso: instationary option is activated but 'Time' column is empty in '%s'. Only simulations can be run on result file, not fitting."%fname)
        else:
            df=df[df["Time"] != ""]
        df_kin=pa.DataFrame()
    else:
        if "Time" in df and sum(df["Time"] == "") == 0:
            werr("parse_miso: we are in stationary case but 'Time' column is not empty in '%s'"%fname)
        df=df[df["Time"] == ""] if "Time" in df else df
    res={"ms": [], "lab": []} # 'peak' and 'mean' are transformed into lab
    # split into kind of measurements: ms, peak, lab
    last_met=last_frag=last_dset=""
    cgr=1
    #import pdb; pdb.set_trace()
    meas_seen=set()
    for kgr, ligr in df.groupby(["Specie", "Fragment", "Dataset"]).groups.items():
        #print("gr=", kgr, ligr)
        met,frag,dset=kgr
        #if met == "M_accoa_c":
        #    import pdb; pdb.set_trace()
        if not met:
            werr("parse_miso: metabolite name is missing in '%s':%d\n%s"%(fname, ist, "\t".join(df.iloc[ligr[0], :])))
        ist=int(df.loc[ligr[0], "iline"])
        iend=int(df.loc[ligr[-1], "iline"])
        mets=np.array([v.strip() for v in met.split("+")]) # met can be A+B+C, take just the first name
        ibad=np.array([v not in clen for v in mets])
        if ibad.any():
            werr("parse_miso: metabolite '%s' was not seen in label transitions, '%s': %d "%(", ".join(mets[ibad]), fname, ist))
        # label length in metabolite(s)
        mlen=sorted(clen[v] for v in mets)
        if mlen[0] != mlen[-1]:
            werr("parse_miso: metabolites of different lengths (%s) are present in '%s', '%s': %d"%(", ".join(str(v) for v in np.unique(mlen)), met, fname, ist))
        mlen=mlen[0]
        # standard fragment form: 1,2,...
        if frag:
            sfrag=np.unique(np.sort(np.array([i for v in frag.split(",") for i in itvl2li(v)]).astype(int)))
            ifr=sfrag-1;
            sfrag=",".join(sfrag.astype(str))
        else:
            ifr=np.arange(mlen)
            sfrag=",".join((ifr+1).astype(str))
        # common sanity check
        #if len(ligr) == 1 and df.loc[ligr, "Isospecies"].iloc[0] != "mean":
        #    werr("parse_miso: a group %s of length 1 is not valid, '%s': %d"%(kgr, fname, ligr[0]+2))
        flen=len(ifr)
        if flen > mlen:
            werr("parse_miso: in group %s, fragment length %d is greater than metabolite length %d in '%s'"%(kgr, flen, mlen, fname))
        if not dset:
            werr("parse_miso: dataset name is missing in '%s':%d\n%s"%(fname, ist, "\t".join(df.iloc[ligr[0], :])))
        if ligr[-1]-ligr[0]+1 != len(ligr):
            werr("parse_miso: measurements %s are not contiguous in '%s'. They occupy rows: %s"%(kgr, fname, ", ".join((ligr+2).astype(str))))
        # ftbl frag
        if flen == mlen:
            ffrag=""
        else:
            ffrag=frag.replace("-", "~")
        val=df.loc[ligr, "Value"].to_numpy()
        sdv=df.loc[ligr, "SD"].to_numpy()
        # detect kind of species
        spec=df.loc[ligr, "Isospecies"]
        if all(spec.str.match("^ *M\d+ *$")):
            kind="ms"
        elif all(spec.str.match("^[ 01x+]+$")):
            kind="lab"
        elif all(spec.str.contains("->")):
            sep="->"
            kind="peak"
        elif all(spec.str.contains("→")):
            sep="→"
            kind="peak"
        elif len(ligr) == 1 and spec.iloc[0] == "mean":
            kind="mean"
        else:
            werr("parse_miso: unknown Isospecies '%s' in group %s in '%s': %d-%d"%(kgr, ", '".join(spec), fname, ist, iend))
        if case_i:
            dsp=dict() # {specie: times indexes}, e.g. "M0": vec("0.1", "0.2", ...)
            spli=[]
            ii0=[]
            dfgr=df.iloc[ligr,:]
            for sp, spi in dfgr.groupby(["Isospecies"]).groups.items():
                dsp[sp]=spi
                spli.append(sp)
                ii0.append(np.where(ligr == spi[0])[0][0])
                # check that all SD are the same for all time points
                u=np.unique(df.loc[spi, "SD"])
                if len(u) != 1:
                    werr(f"parse_miso: SD must be the same at all time points for {kgr}, {sp}: '{fname}': "+", ".join(df.loc[spi, "iline"]))
            ii0=sorted(ii0)
        if kind == "ms":
            # ms group here, like M0, M1
            #print("ms gr=", kgr)
            w=np.char.lstrip(df.loc[ligr, "Isospecies"].to_numpy().astype("str"), " M").astype(int)
            # ms sanity check
            if any(w > flen):
                werr("parse_miso: invalid MS weight '%s' in group %s in '%s':%d-%d"%(w[w>flen].astype(str)[0], kgr, fname, ist, iend))
            if not case_i and len(ligr) > flen+1:
                raise Excpetion("parse_miso: too many MS entries %d (max %d expected) in group %s in '%s':%d-%d"%(kgr, len(ligr) , flen+1, fname, ist, iend))
            # add this group to results
            #if frag == "":
            #    frag=",".join(str(i) for i in range(1,flen+1))
            if case_i:
                #if met == "Phe":
                #    import pdb; pdb.set_trace()

                res["ms"] += [f"\t{met}\t{ffrag}\t{w[0]}\tNA\t{sdv[0]}"+"   // %s: %d"%(fname, ist)]
                res["ms"] += [f"\t\t\t{w[i0]}\tNA\t{sdv[i0]}"+"   // %s: %s"%(fname, df.loc[ligr[i0], "iline"]) for i,i0 in zip(range(1, len(spli)), ii0[1:])]
                #import pdb; pdb.set_trace()
                if not dfgr[dfgr["Time"] != ""].empty:
                    for sp,spi in dsp.items():
                        df_kin=pa.concat([df_kin, pa.DataFrame(df.loc[spi, "Value"].to_numpy().reshape(1, -1), columns=df.loc[spi, "Time"], index=[f"m:{met}:{ffrag}:{sp[1:]}:{df.loc[spi[0],'iline']}"])])
            else:
                res["ms"] += [f"\t{met}\t{ffrag}\t{w[0]}\t{val[0]}\t{sdv[0]}"+"   // %s: %d"%(fname, ist)]
                res["ms"] += [f"\t\t\t{w[i]}\t{val[i]}\t{sdv[i]}"+"   // %s: %s"%(fname, df.loc[ligr[i], "iline"]) for i in range(1, len(ligr))]
        elif kind == "lab" or kind == "peak" or kind == "mean":
            # label group (like 01x+1x1)
            if kind == "lab":
                labs=[[vv.strip() for vv in v.split("+")] for v in df.loc[ligr, "Isospecies"]]
            elif kind == "peak":
                b=np.ones(mlen, str) # base where lab will be injected
                b.fill("x")
                b[ifr]="0"
                labs=[]
                for peak in spec:
                    li=np.array([v.strip() for v in peak.replace(sep, ",").split(",") if v.strip()]).astype(int)-1
                    tmp=b.copy()
                    tmp[li]="1"
                    labs.append(["".join(tmp[ifr])])
                #import pdb; pdb.set_trace()
            elif kind == "mean":
                b=np.ones(mlen, str) # base where lab will be injected
                b.fill("x")
                labs=[["x"*i+"1"+"x"*(flen-i-1) for i in range(flen)]]
                if flen > 1:
                    val[val != 'NA'] = vsadd(vsadd("(", val[val != 'NA']), ")*"+str(flen))
                    sdv = vsadd(vsadd("(", sdv), ")*"+str(flen))
            # label sanity check
            for i,li in enumerate(labs):
                for v in li:
                    if len(v) != flen:
                        werr("parse_miso: entry '%s' has length %d different from fragment length %d, '%s': %s"%(v, len(v), flen, fname, df.loc[ligr[i], "iline"]))
            # inject fragment species into full molecule
            if flen < mlen:
                b=np.ones(mlen, str) # base where lab will be injected
                b.fill("x")
                for li in labs:
                    for i in range(len(li)):
                        tmp=b.copy()
                        tmp[ifr]=list(li[i])
                        li[i]="".join(tmp)
            # normalize or not?
            collab=sorted(v for li in labs for v in li) # will be collapsed labels.
            # the group is normalizable if collapsed labs is composed of only "x"
            while True:
                collab=sorted(collab) # if there are two collapsible labs, they will be neighbors
                found=False
                for i in range(len(collab)):
                    if i == len(collab)-1:
                        break
                    s=np.array(list(collab[i])) if i == 0 else n # string
                    n=np.array(list(collab[i+1])) # next
                    cdif=s != n
                    if cdif.sum() == 1:
                        idif=np.where(cdif)[0][0]
                        found=s[idif] == "0" and n[idif] == "1"
                        if found:
                            s[idif]="x"
                            collab[i]="".join(s)
                            del(collab[i+1])
                            break
                if not found or len(collab) == 1:
                    break
            if len(collab) == 1 and collab[0] == "x"*len(collab) and not any(val == "NA"):
                norma=True
            else:
                norma=False
            # add this group to results
            #	META_NAME	CUM_GROUP	VALUE	DEVIATION	CUM_CONSTRAINTS
            if case_i:
                res["lab"] += [f"\t{met}\t1\t{val[0]}\t{sdv[0]}\t"+"+".join("#"+v for v in labs[0])+"   // %s: %d"%(fname, ist)]
                res["lab"] += [f"\t\t{i+1 if norma else 1}\t{val[i0]}\t{sdv[i0]}\t"+"+".join("#"+v for v in labs[i0])+"   // %s: %s"%(fname, df.loc[ligr[i0], "iline"]) for i,i0 in zip(range(1, len(spli)), ii0[1:])]
                if not dfgr[dfgr["Time"] != ""].empty:
                    for i,(sp,spi) in enumerate(dsp.items()):
                        #import pdb; pdb.set_trace()
                        df_kin=pa.concat([df_kin, pa.DataFrame(df.loc[spi, "Value"].to_numpy().reshape(1, -1), columns=df.loc[spi, "Time"], index=[f"l:{met}:{'+'.join('#'+v for v in labs[spi[0]])}:NA"])])
            else:
                if met != last_met or frag != last_frag:
                    last_met=met
                    last_frag=frag
                elif dset != last_dset:
                    last_dset=dset
                res["lab"] += [f"\t{met}\t{cgr}\t{val[0]}\t{sdv[0]}\t"+"+".join("#"+v for v in labs[0])+"   // %s: %d"%(fname, ist)]
                res["lab"] += [f"\t\t{i+cgr if norma else cgr}\t{val[i]}\t{sdv[i]}\t"+"+".join("#"+v for v in labs[i])+"   // %s: %s"%(fname, df.loc[ligr[i], "iline"]) for i in range(1, len(ligr))]
                cgr += len(ligr) if norma else 1
    return (res, df_kin) if case_i else res
def parse_linp(f, clen={}):
    "Parse label input TSV file. Return a list of lines to add to ftbl"
    if "name" in dir(f):
        fname=f.name
    else:
        fname=f
    fname=os.path.basename(fname)
    
    df=tsv2df(f)
    if "Metabolite" in df.columns:
        df.rename({"Metabolite": "Specie"}, inplace=True, axis=1)

    #print("df=", df)
    res=[]
    for met, ligr in df.groupby(["Specie"]).groups.items():
        mlen=clen.get(met)
        il=df.loc[ligr, "iline"].to_numpy() # strings
        # sanity check
        if clen:
            if mlen is None:
                werr("parse_linp: specie '%s' was not seen in label transitions, '%s': %s"%(met, fname, ", ".join(il)))
            ibad=np.where([len(v) != mlen for v in df.loc[ligr, "Isotopomer"]])[0]
            if len(ibad):
                werr("parse_linp: for specie '%s', isotopomer length(s) (%s) differ from its labeling atom length %d in '%s': %s"%(met, ", ".join(str(len(v)) for v in df.loc[ligr[ibad], "Isotopomer"]), mlen, fname, ", ".join(il)))
            if len(ligr) > 2**mlen:
                werr("parse_linp: too many isotopomers %d for metabolite '%s' in '%s':%s-%s"%(len(ligr), met, fname, il[0], il[-1]))
        if len(ligr) != ligr[-1]-ligr[0]+1:
            werr("parse_linp: for metabolite '%s', label input is not continguous in '%s': %s-%s"%(met, fname, il[0], il[-1]))
        res.append("\t%s\t#%s\t%s   // %s: %s"%(met, df.loc[ligr[0], "Isotopomer"], df.loc[ligr[0], "Value"], fname, il[0]))
        res += ["\t\t#%s\t%s   // %s: %s"%(df.loc[ligr[i], "Isotopomer"], df.loc[ligr[i], "Value"], fname, il[i]) for i in range(1, len(ligr))]
    return res
def parse_mflux(f, dfl={}):
    "Parse flux measurements TSV file. Return a list of lines to add to ftbl"
    if "name" in dir(f):
        fname=f.name
    else:
        fname=f
    fname=os.path.basename(fname)
    
    df=tsv2df(f)
    #print("df=", df)
    res=[]
    for flux,ligr in df.groupby(["Flux"]).groups.items():
        il=df.loc[ligr, "iline"].to_numpy() # strings
        # sanity check
        if len(ligr) > 1:
            werr("parse_mflux: flux '%s' has more than 1 measurement in '%s': %s"%(flux, fname, ", ".join(il)))
        if dfl and flux not in dfl:
            werr("parse_mflux: flux '%s' was not seen in metabolic network, '%s': %s"%(flux, fname, il[0]))
        res.append("\t%s\t%s\t%s   // %s: %s"%(flux, df.loc[ligr[0], "Value"], df.loc[ligr[0], "SD"], fname, il[0]))
    return res
def parse_opt(f):
    "Parse options TSV file. Return a list of lines to add to ftbl"
    if "name" in dir(f):
        fname=f.name
    else:
        fname=f
    fname=os.path.basename(fname)
    
    df=tsv2df(f).loc[:, ["Name", "Value"]].to_numpy().astype(str)
    res=vsadd(vsadd("\t", vsadd(df[:,0], "\t")), df[:,1]).tolist()
    return res
def parse_tvar(f, dfl={}, itnl_met=set()):
    "Parse variable type TSV file. Return a tuple of a dict and a list with lines to add to ftbl"
    if "name" in dir(f):
        fname=f.name
    else:
        fname=f
    fname=os.path.basename(fname)
    
    df=tsv2df(f)
    #print("df=", df)
    fl={}
    mets=[]
    for kgr,ligr in sorted(df.groupby(["Kind", "Name"]).groups.items(), key=lambda t: (t[0][0], plain_natural_key(t[0][1]))):
        il=df.loc[ligr, "iline"].to_numpy() # strings
        kind,nm=kgr
        #if kind == "NET":
        #    import pdb; pdb.set_trace();
        # sanity check
        if kind not in ("NET", "XCH", "METAB"):
            werr("parse_tvar: kind '%s' is unknown (expected one of NET, XCH or METAB) in '%s': %s"%(kind, ", ".join(il)))
        if len(ligr) > 1:
            werr("parse_tvar: groupe '%s' has more than 1 entry in '%s': %s"%(kgr, ", ".join(il)))
        val=df.loc[ligr[0], "Value"]
        ty=df.loc[ligr[0], "Type"]
        if ty not in ("F", "C", "D"):
            werr("parse_tvar: type '%s' is not valid (expected one of F, D or C) in '%s': %s"%(ty, il[0]))
        if ty != "D" and not val:
            werr("parse_tvar: type '%s' supposes non empty value in '%s': %s"%(ty, il[0]))
        if kind == "METAB":
            if itnl_met and nm not in itnl_met:
                werr("parse_tvar: metabolite '%s' was not seen in internal metabolites, '%s': %s"%(nm, il[0]))
            if ty != "C":
                val="-"+val
            mets.append("\t%s\t%s   // %s: %s"%(nm, val, fname, il[0]))
        else:
            if dfl and nm not in dfl:
                werr("parse_tvar: flux '%s' was not seen in metabolic network, '%s': %s"%(nm, il[0]))
            fl[kind]=fl.get(kind, [])
            fl[kind].append("\t\t%s\t%s\t%s   // %s: %s"%(nm, ty, val, fname, il[0]))
    return (fl, mets)
def parse_cnstr(f):
    "Parse constraint TSV file. Return a tuple of 2 dicts (eq, ineq) with lines to add to ftbl"
    if "name" in dir(f):
        fname=f.name
    else:
        fname=f
    fname=os.path.basename(fname)
    
    df=tsv2df(f)
    #print("df=", df)
    eq={}
    ineq={}
    for kgr,ligr in df.groupby(["Kind", "Formula", "Operator"]).groups.items():
        il=df.loc[ligr, "iline"].to_numpy() # strings
        kind,frml,op=kgr
        # sanity check
        if op != "==" and not op in invcomp:
            werr("parse_cnstr: operator '%s' is unknown in '%s': %s"%(op, ", ".join(il)))
        if kind not in ("NET", "XCH", "METAB"):
            werr("parse_cnstr: kind '%s' is unknown (expected one of NET, XCH or METAB) in '%s': %s"%(kind, ", ".join(il)))
        if len(ligr) > 1:
            werr("parse_cnstr: formula '%s' has more than 1 entry in '%s': %s"%(frml, ", ".join(il)))
        val=df.loc[ligr[0], "Value"]
        if op == "==":
            eq[kind]=eq.get(kind, [])
            eq[kind].append("\t\t%s\t%s   // %s: %s"%(val, frml, fname, il[0]))
        else:
            op=invcomp[op]
            ineq[kind]=ineq.get(kind, [])
            ineq[kind].append("\t\t%s\t%s\t%s   // %s: %s"%(val, op, frml, fname, il[0]))
    #import pdb; pdb.set_trace()
    return (eq, ineq, df)
def parse_mmet(f, smet=set()):
    "Parse metabolite concentration measurements TSV file. Return a list of lines to add to ftbl"
    if "name" in dir(f):
        fname=f.name
    else:
        fname=f
    fname=os.path.basename(fname)
    
    df=tsv2df(f)
    if "Metabolite" in df.columns:
        df.rename({"Metabolite": "Specie"}, inplace=True, axis=1)
    #print("df=", df)
    res=[]
    for met,ligr in df.groupby(["Specie"]).groups.items():
        il=df.loc[ligr, "iline"].to_numpy() # strings
        # sanity check
        if len(ligr) > 1:
            werr("parse_mmet: specie '%s' has more than 1 measurement in '%s': %s"%(met, fname, ", ".join(il)))
        if smet and any(m.strip() not in smet for m in met.split("+")):
            werr("parse_mmet: specie '%s' was not seen in internal metabolites of network, '%s': %s"%(met, fname, il[0]))
        res.append("\t%s\t%s\t%s   // %s: %s"%(met, df.loc[ligr[0], "Value"], df.loc[ligr[0], "SD"], fname, il[0]))
    return res

def compile(mtf, cmd, case_i=False, clen=None):
    "Compile FTBL content from mtf names: netw, miso etc. Return a dict of ftbl lines"
    # dict of ftbl sections. Contains list of lines to be completed by compilation
    dsec={
        "proj": [
            "PROJECT",
            "	NAME	VERSION	FORMAT	DATE	COMMENT",
        ],
        "netw": [
            "NETWORK",
            "	FLUX_NAME	EDUCT_1	EDUCT_2	PRODUCT_1	PRODUCT_2",
        ],
        "notr": [
            "NOTRACER_NETWORK",
            "	FLUX_NAME	EQUATION",
        ],
       "flux": [
            "FLUXES",
            {
                "NET": [
                    "	NET",
                    "		NAME	FCD	VALUE(F/C)	ED_WEIGHT	LOW(F)	INC(F)	UP(F)",
                ],
                "XCH": [
                    "	XCH",
                    "		NAME	FCD	VALUE(F/C)	ED_WEIGHT	LOW(F)	INC(F)	UP(F)",
                ],
            },
        ],
        "met_pool": [
            "METABOLITE_POOLS",
            "	META_NAME	META_SIZE",
        ],
        "eq": [
            "EQUALITIES",
            {
                "NET": [
                    "	NET",
                    "		VALUE	FORMULA",
                ],
                "XCH": [
                    "	XCH",
                    "		VALUE	FORMULA",
                ],
                "METAB": [
                    "	METAB",
                    "		VALUE	FORMULA",
                ],
            },
        ],
        "ineq": [
            "INEQUALITIES",
            {
                "NET": [
                    "	NET",
                    "		VALUE	COMP	FORMULA",
                ],
                "XCH": [
                    "	XCH",
                    "		VALUE	COMP	FORMULA",
                ],
                "METAB": [
                    "	METAB",
                    "		VALUE	COMP	FORMULA",
                ],
            },
        ],
        "linp": [
            "LABEL_INPUT",
            "	META_NAME	ISOTOPOMER	VALUE",
        ],
        "mflux": [
            "FLUX_MEASUREMENTS",
            "	FLUX_NAME	VALUE	DEVIATION",
        ],
        "meas_lab": [
            "LABEL_MEASUREMENTS",
            "	META_NAME	CUM_GROUP	VALUE	DEVIATION	CUM_CONSTRAINTS",
        ],
        "meas_peak": [
            "PEAK_MEASUREMENTS",
            "	META_NAME	PEAK_NO	VALUE_S	VALUE_D-	VALUE_D+	VALUE_DD	VALUE_T	DEVIATION_S	DEVIATION_D-	DEVIATION_D+	DEVIATION_DD/T",
        ],
        "meas_ms": [
            "MASS_SPECTROMETRY",
            "	META_NAME	FRAGMENT	WEIGHT	VALUE	DEVIATION",
        ],
        "mmet": [
            "METAB_MEASUREMENTS",
            "	META_NAME	VALUE	DEVIATION",
        ],
        "opt": [
            "OPTIONS",
            "	OPT_NAME	OPT_VALUE",
        ],
    }
    dsec_empty=dsec.copy()
    # Parse netw file if not empty
    if "netw" in mtf and mtf["netw"]:
        pth=try_ext(mtf["netw"], ["netw", "txt"])
        (netw, notr_netw, eqs, ineqs, fluxes, (m_left, m_right), sto, dclen)=txt_parse(pth)
        # build afl matrix: each row is a balance on an internal metabolite, each column is a flux values
        nb_flux=len(sto)
        nm_flux=sorted(sto.keys(), key=plain_natural_key)
        m_l=set(m_left.keys())
        m_r=set(m_right.keys())
        m_inp=m_l-m_r
        m_out=m_r-m_l
        itnal_met=(m_l|m_r) - m_inp - m_out
        nm_met=sorted(itnal_met, key=plain_natural_key)
        fl2i=dict((fl, i) for i, fl in enumerate(nm_flux))
        met2i=dict((m, i) for i, m in enumerate(nm_met))
        afl=np.zeros((len(nm_met)+len(eqs[0]), len(nm_flux)))
        for fl, row in sto.items():
            for m, c in row.items():
                if m not in met2i:
                    continue
                afl[met2i[m], fl2i[fl]]=c
        # add netw equations
        for i,eq in enumerate(eqs[0]):
            fls=eq[1].split(" - ")
            afl[len(nm_met)+i, fl2i[fls[0]]]=1.
            afl[len(nm_met)+i, fl2i[fls[1]]]=-1.
        if nb_flux > 1:
            # qr(afl)
            q,r,p=linalg.qr(afl, pivoting=True)
            d=diag(r)
            if d[0] == 0.:
                werr("Stoechiometrix matrix is 0")
            rank=sum(abs(d/d[0]) > 1.e-10)
            # free fluxes are the last in pivots
            ff=set(nm_flux[i] for i in p[rank:])
        else:
            ff=[]
        # prepare header
        base=mtf["netw"].name
        dsec["proj"] += ["	%s	%s			command='%s'"%(base, version, cmd)]
        # prepare the ftbl content
        for row in netw:
            if type(row) == str:
                dsec["netw"] += ["%s"%row.strip("\n\r")]
                continue
            dsec["netw"] += ["\t%s"%row[0][0]] # reac name
            # input metabs
            carbs="\t"
            imetab=0
            for (coef, metab, carb) in row[1]:
                imetab=imetab+1
                dsec["netw"][-1] += "\t%s%s"%(coef+"*" if coef and coef != "1" else "", metab)
                carbs=carbs+"\t%s"%(("#"+carb) if metab else "")
            if imetab==1:
                # complete by empty metabolite
                dsec["netw"][-1] += "\t"
                carbs=carbs+"\t"
            # output metabs
            imetab=0
            for (coef, metab, carb) in row[2]:
                imetab=imetab+1
                dsec["netw"][-1] += "\t%s"%metab
                carbs=carbs+"\t%s"%(("#"+carb) if metab else "")
            dsec["netw"] += ["%s"%carbs]
        for row in notr_netw:
            if type(row) == str:
                dsec["notr"] += ["%s"%row]
                continue
            dsec["notr"] += ["\t%s\t%s"%(row[0][0], " = ".join("+".join((c+"*" if c and c != "1" else "")+m  for c,m,t in side) for side in row[1:3]))]
        for e in eqs[0]:
            dsec["eq"][1]["NET"] += ["\t\t%s\t%s"%e]
        for e in eqs[1]:
            dsec["eq"][1]["XCH"] += ["\t\t%s\t%s"%e]
        for ine in ineqs[0]:
            dsec["ineq"][1]["NET"] += ["\t\t%s\t%s\t%s"%ine]
        if not "linp" in mtf:
            # add default full label input
            for m in sorted(m_inp, key=plain_natural_key):
                dsec["linp"] += ["\t%s\t#%s\t1.0"%(m, "1"*dclen[m])]
    else:
        dclen={} if clen is None else clen
        itnal_met=set()
        sto={}
        fluxes=[]
    sfl=set(sto.keys()) # set of fluxes
    if "cnstr" in mtf and mtf["cnstr"]:
        pth=try_ext(mtf["cnstr"], ["cntsr", "tsv", "txt"])
        ce,ci,df=parse_cnstr(pth)
        for k,v in ce.items():
            dsec["eq"][1][k] += v
        for k,v in ci.items():
            dsec["ineq"][1][k] += v
        # complete flux set by those from eq:net
        
        #import pdb; pdb.set_trace()
        for eq in df[(df["Kind"] == "NET") & (df["Operator"] == "==")]["Formula"]:
            sfl |= set(formula2dict(eq).keys())
    if "miso" in mtf and mtf["miso"]:
        pth=try_ext(mtf["miso"], ["miso", "tsv", "txt"])
        if case_i:
            meas, df_kin=parse_miso(pth, dclen, case_i)
            #import pdb; pdb.set_trace()
        else:
            meas=parse_miso(pth, dclen)
        dsec["meas_lab"] += meas["lab"]
        dsec["meas_ms"] += meas["ms"]
    else:
        if case_i:
            df_kin=pa.DataFrame() # return empty data frame for kinetic measurements
    # simple sections
    if "linp" in mtf and mtf["linp"]:
        pth=try_ext(mtf["linp"], ["linp", "tsv", "txt"])
        dsec["linp"] += parse_linp(pth, dclen)
    if "mflux" in mtf and mtf["mflux"]:
        pth=try_ext(mtf["mflux"], ["mflux", "tsv", "txt"])
        dsec["mflux"] += parse_mflux(pth, sfl)
    if "mmet" in mtf and mtf["mmet"]:
        pth=try_ext(mtf["mmet"], ["mmet", "tsv", "txt"])
        dsec["mmet"] += parse_mmet(pth, itnal_met)
    if "opt" in mtf and mtf["opt"]:
        pth=try_ext(mtf["opt"], ["opt", "tsv", "txt"])
        dsec["opt"] += parse_opt(pth)
    # with subsections NET/XCH/...
    if "tvar" in mtf and mtf["tvar"]:
        pth=try_ext(mtf["tvar"], ["tvar", "tsv", "txt"])
        tf,tm=parse_tvar(pth)
        stvar=dict((nx, set(v.split("\t")[2] for v in li)) for nx,li in tf.items())
        if "NET" in tf:
            dsec["flux"][1]["NET"] += tf["NET"]
        if "XCH" in tf:
            dsec["flux"][1]["XCH"] += tf["XCH"]
        dsec["met_pool"] += tm
    else:
        stvar={"NET": set(), "XCH": set()}
        dsec["met_pool"] += ["\t"+m+"\t0.1" for m in itnal_met]

    snrev=set() # set of non reversible fluxes
    for tpl in fluxes:
        f, rev, imposed_sens, fd=tpl
        #if f == "R_FORt":
        #    import pdb; pdb.set_trace()
        if not rev:
            dsec["flux"][1]["XCH"] += ["\t\t%s\tC\t0"%f]
            snrev.add(f)
        elif not "tvar" in mtf:
            dsec["flux"][1]["XCH"] += ["\t\t%s\tF\t0.01"%f]
    badf=stvar.get("XCH", set()) & snrev
    if badf:
        werr("following fluxes should not appear in '%s', with 'XCH' kind as they are non reversible in '%s':\n\t'%s'"%(mtf["tvar"], mtf["netw"], "'\n\t'".join(badf)))
    # f, rev, imposed_sens, fd=tpl
    dtnet=dict((tpl[0], "\t\t%s\t%s\t0.2E0"%(tpl[0], "F" if tpl[0] in ff else "D")) for tpl in fluxes)
    dtxch=dict((tpl[0], "\t\t%s\t%s"%(tpl[0], "%s\t0.01E0"%tpl[3])) for tpl in fluxes if not tpl[1])
    
    # fluxes that are not in tvar are completed from dtnet and dtxch
    dsec["flux"][1]["NET"] += [v for k,v in dtnet.items() if k not in stvar.get("NET", set())]
    dsec["flux"][1]["XCH"] += [v for k,v in dtxch.items() if k not in (stvar.get("XCH", set()) | snrev)]
    #import pdb; pdb.set_trace()
    return (dsec, dclen) if not case_i else (dsec, dclen, df_kin)
def main(argv=sys.argv[1:], res_ftbl=None):
    ord_args=[]
    class ordAction(argparse.Action):
        def __call__(self, parser, namespace, values, option_string=None):
            ord_args.append((self.dest, values))
            setattr(namespace, self.dest, values)
    parser=argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("--mtf", action=ordAction, help=
"""MTF is a coma separated list of files with following extensions/meanings:
 netw: a text file with stoichiometric reactions and label transitions (one per line)
    Comments starts with '#' but those starting with '###' introduce 
    pathways which are numbered as well as reactions in them. Reaction 
    name can precede the reaction itself and is separated by ":" If no 
    explicit name is given, reactions in FTBL file will be named 
    according a pattern 'rX.Y' where X is pathway number and Y is 
    reaction number in the pathway. But it is highly recommended to 
    give explicit names to reactions.
    Symbols "+", "(", ")" and ":" are not allowed in metabolite neither reaction names
    Example of reaction and label transition:
       edd: Gnt6P (ABCDEF) -> Pyr (ABC) + GA3P (DEF)
    Non reversible reactions are signaled with '->' (as in the example above).
    A sign '<->' can be used for reversible reactions.
    If 'netw' name is equal to '-', then its content is read from standard input, e.g.
      '--mtf netw=-'
 linp: label inputs (starting from this extensions, TSV files are assumed)
 miso: isotopic measurements (NMR (label, peak) and MS)
 mflux: flux measurements
 mmet: metabolite concentration measurements
 tvar: type of variables (NET or XCH , free or dependent, starting values, ...)
 cnstr: equality and inequality constraints on fluxes and concentrations
 opt: options
 ftbl: name of output FTBL file. If not given, it will be equal to 'netw'
   stem with '.ftbl' extension. If it is equal to '-', then the result 
   will be written to standard output, e.g.
     'txt2ftbl --mtf ftbl=-,ecoli.netw'
   Intermediate directories in ftbl path are silently created if non existent.
 vmtf: variable part of mtf approach.
   If a series of FTBL files has to be generated partially with 
   information common to all files (constant part) and partially with 
   sections proper to each FTBL (variable part) then files containing 
   variable sections (e.g. 'miso') can be given in a special file 
   having an extension (or prefix, cf. hereafter) 'vmtf'. In such a 
   way, each FTBL file will be produced from combination of MTF files 
   given directly in this option (constant part) and files given on a 
   corresponding row of 'vmtf' file. vmtf file is a TSV file with 
   columns using the same names: 'netw', 'linp', etc. Each row contains 
   file names that will be used to produce an FTBL file. Thus each row 
   must have 'ftbl' column with unique and non empty name. When 'vmtf' 
   is used, 'ftbl' cannot be present on the command line. If a file 
   type is present both in column names of 'vmtf' and in '--mtf' option 
   then the content of 'vmtf' file will take precedence. Empty values 
   in 'vmtf' file are ignored. All file paths in 'vmtf' file are 
   considered relative to the location of 'vmtf' file itself.
Only first 3 files are necessary to obtain a workable FTBL file, others 
are optional.
Example: 'txt2ftbl --mtf ecoli.netw,glu08C1_02U.linp,cond1.miso,cond1.mflux'
NB: no space is allowed around comas. If a file path has a spaces in 
its name, it must be enclosed into quotes or double quotes.
If an entry file cannot be renamed to have some of these extensions, then
they can be used as prefixes followed by a '=' sign, e.g.
  'txt2ftbl --mtf netw=ecoli.txt,linp=glu08C1_02U.tsv,cond1.miso,cond1.mflux'
As you can see from this example, both naming schemes can be mixed.
If for some reason, the same type of file is indicated several times
(no matter with extension or prefix), the last occurrence supersedes
all precedent ones.
""")
    parser.add_argument("--prefix", action=ordAction, help=
"""If all input files have the same name pattern and are different only 
in extensions then the pattern can be given as PREFIX, e.g.
  '--prefix somedir/ecoli'
Then in 'somedir', we suppose to have 'ecoli.netw', 'ecoli.linp' and 
other input files having names starting with 'ecoli' and ending with 
corresponding extensions.
NB. If some file is given in more than one option: '--prefix' and/or 
'--mtf' then the last occurrence overrides precedent ones.
""")
    parser.add_argument("--eprl", action="append", help=
"""Parallel experiments can be given with this option. It must
introduce a couple of linp/miso files and optional auxiliary ftbl name. These files
correspond to a given parallel experiment. This option can be repeated as
many times as there are additional parallel experiments, e.g.
  'txt2ftbl --mtf ec.netw,glc6.linp,glc6.miso --eprl glc1.linp,glc1.miso --eprl glc4.linp,glc4.miso'
This command will produce a main FTBL file 'ec.ftbl' including all necessary
sections (NETWORK, etc.) but also two auxiliary FTBL files: 'glc1.ftbl' and
'glc4.ftbl' having only label input/measurement sections. They will correspond
to 2 additional parallel experiments. If ftbl file is not given in --eprl
option, the name of miso file will be used for it. If intermediate
directories in ftbl path are non existent they will be silently created.
Auxiliary ftbl names will be put in 'OPTIONS/prl_exp' field on the main ftbl file.
These names will be written there in a form relative to the main ftbl.
To shorten the writings, it is possible to indicate only one of two .miso/.linp files.
The other one will be guessed if it has canonical extension. If extension is omitted then .miso and .linp files are searched with these extensions. In this case, several parallel experiments can be given with one --eprl option. So that above example can be shorten to:
  'txt2ftbl --mtf ec.netw,glc6.linp,glc6.miso --eprl glc1,glc4'
""")
    parser.add_argument("--inst", action="store_true", default=False, help=
"""Prepare FTBL for instationary case. File 'netw' is supposed to have 
column 'Time' non empty. Isotopic kinetic data will be written to a TSV 
file with 'ikin' extension. Its name will be the same as in FTBL file, 
and FTBL field 'OPTIONS/file_labcin' will contain 'ikin' file name.
""")
    parser.add_argument("--force", action="store_true", default=False, help=
"""Overwrite an existent result file not produced by this script.
NB. If a result file exists and is actually produced by this script, 
then it is silently overwritten even without this option. The script 
detects if it was the creator of a file by searching for a string "// 
Created by 'txt2ftbl" at the first line of the file. By removing or 
editing this comment, user can protect a file from a silent 
overwriting.
""")
    parser.add_argument("netw", default="", nargs="?", help="""
If 'netw' file is not given in any option (neither --mtf nor --prefix), it 
can be given as the only argument NETW, e.g.
  txt2ftbl ecoli.txt
or
  txt2ftbl --mtf ms_nmr_data.miso,glucose.linp ecoli.txt
If 'netw' file name is given both in any option and as an argument, it 
is the argument value that will take precedence.
""")
    if len(argv) == 0:
        parser.print_usage(sys.stderr)
        return 1
    #print("opts=", vars(opts))
    #print("ord=", ord_args)
    # default values
    mtfsuf={"netw", "linp", "miso", "mflux", "mmet", "tvar", "cnstr", "opt", "vmtf", "ftbl"}
    mtf={} # multiplex tsv files to compile
    prl=[] # parallel experiments, 2 or 3-tuples linp+miso(+ftbl)
    # get arguments
    opts = parser.parse_args(argv)
    force=opts.force
    netw=opts.netw
    case_i=opts.inst
    for o,a in ord_args:
        if o == "mtf":
            # make dict {"miso": <file_path>, "netw": ...}
            # a is like "f.netw,exp1.miso,...."
            for v in a.split(","):
                v=v.strip()
                if not v:
                    continue
                if "=" in v:
                    ty,nm=v.split("=", 1)
                    mtf[ty]=nm
                else:
                    mtf[Path(v).suffix[1:]]=v
        elif o == "prefix":
            # a is like "/some/path/file_stem"
            f=Path(a)
            if f.is_dir():
                d=f
                stem="*"
            else:
                d=f.parent
                stem=f.name
            if not d.is_dir():
                werr("directory '%s' from --prefix does not exist"%str(d))
            nfound=0
            #print("d=", d, "\tstem=", stem)
            for suf in mtfsuf-{"ftbl"}:
                li=list(d.glob(stem+"."+suf))
                #print("suf=", suf, "\tli=", li)
                if len(li) > 1:
                    werr("multiple .%s files found:\n\t'%s'"%(suf, "'\n\t'".join(str(v) for v in li)))
                if li:
                    mtf[suf]=str(li[0])
                    nfound += 1
            if nfound == 0:
                werr("No MTF file found with prefix '%s'"%a)
    if "vmtf" in mtf and "ftbl" in mtf:
        werr("'ftbl' and 'vmtf' cannot be simultaneously present in '--mtf' option")
    if netw:
        mtf["netw"]=netw # overwrite previous setting if any
    if opts.eprl and "vmtf" in mtf:
        werr("Option --eprl cannot be used simultaneously with 'vmtf' entry in --mtf")
    if opts.eprl is None:
        opts.eprl=[]
    if "netw" in mtf:
        if mtf["netw"] == "-":
            mtf["netw"]=sys.stdin
        else:
            mtf["netw"]=Path(mtf["netw"])
    # what kind of output we have?
    if "ftbl" in mtf:
        p=Path(mtf["ftbl"])
        if p.suffix == ".ftbl":
            mtf["ftbl"]=p
        elif mtf["ftbl"] == "-":
            mtf["ftbl"]=sys.stdout
        else:
            mtf["ftbl"]=p.parent/(p.name+".ftbl")
    elif "netw" in mtf and not "vmtf" in mtf:
        if mtf["netw"] == sys.stdin:
            mtf["ftbl"]=sys.stdout
        else:
            mtf["ftbl"]=mtf["netw"].with_suffix(".ftbl")
    elif not "vmtf" in mtf:
        mtf["ftbl"]=sys.stdout
    # get opt in preamble
    if "opt" in mtf:
        dfopt=tsv2df(mtf["opt"])
    # prepare prl
    # if no --eprl, get prl_exp from .opt
    if not opts.eprl and "opt" in mtf:
        wd=Path(mtf["opt"]).parent
        fli=[v.replace(";", ",") for v in dfopt.loc[dfopt["Name"] == "prl_exp", "Value"]]
        #print("fli=", fli)
        fli=[str(wd/s.strip()) for v in fli for s in v.split(",") if s.strip()]
        #print("fli2=", fli)
    else:
        fli=opts.eprl
    #import pdb; pdb.set_trace()
    
    for t in fli: # prepare prl list
        for v in t.split(","):
            d={}
            v=v.strip()
            if not v:
                continue
            if "=" in v:
                ty,nm=v.split("=", 1)
            else:
                ty=Path(v).suffix[1:]
                nm=v
            if ty not in ("linp", "miso", "opt", "ftbl"):
                nm=try_ext(v, ["miso"])
                if not nm.is_file():
                    werr("option --eprl expects in argument a list of linp, miso and optionally auxiliary ftbl files instead got type '%s' in '%s"%(ty, t))
                ty="miso"
            d[ty]=nm
            if "linp" not in d:
                p=Path(d.get("miso", "")).with_suffix(".linp")
                if p.is_file():
                    d["linp"]=p
                else:
                    werr("'linp' file was not found for --eprl option '%s'"%t)
            if "miso" not in d:
                p=Path(d.get("linp", "")).with_suffix(".miso")
                if p.is_file():
                    d["miso"]=p
                else:
                    werr("'miso' file was not found for --eprl option '%s'"%t)
            if "opt" not in d:
                p=Path(d.get("miso", "")).with_suffix(".opt")
                if p.is_file():
                    d["opt"]=p
            if "ftbl" not in d:
                d["ftbl"]=str(Path(d["miso"]).with_suffix(".ftbl"))
            d["ftbl"]=Path(d["ftbl"]).with_suffix(".ftbl")
            prl.append(d)
    #print("prl=", prl)

    cmd=f"{me} "+' '.join(v.replace(' ', r'\ ') for v in argv)
    scre=f"// Created by '{cmd}'"
    
    if "vmtf" in mtf:
        vdf=tsv2df(mtf["vmtf"])
        dftbl=Path(mtf["vmtf"]).parent
    else:
        vdf=pa.DataFrame({"ftbl": [mtf["ftbl"]], "iline":["NA"]})
        dftbl=None
        del(mtf["ftbl"])
    if "ftbl" not in vdf:
        werr("'ftbl' column must be present in '%s'"%mtf["vmtf"])
    # check case_i in opt
    if not case_i and "opt" in mtf:
        tmp=dfopt.loc[dfopt["Name"] == "file_labcin", "Value"]
        if len(tmp) == 1:
            warn(f"instationary mode is activated as 'file_labcin' is not empty in '{mtf['opt']}'")
            case_i=True
    # run through all ftbls
    if dftbl is not None:
        # add dftbl to all fields in vmtf
        vdf[vdf != ""]=[dftbl/v for v in vdf[vdf != ""].values.flatten() if v == v]
        
    for ftbl,ligr in vdf.groupby(["ftbl"]).groups.items():
        il=vdf.loc[ligr, "iline"].to_numpy() # strings
        # sanity check
        if len(ligr) > 1:
            werr("'ftbl' column has repeated values in '%s': %s"%(mtf["vmtf"], ", ".join(il)))
        # prepare running mtf, full mtf for one row
        rmtf=mtf.copy()
        rmtf.update((k,v) for k,v in vdf.iloc[ligr[0], :].to_dict().items() if v)

        # what kind of output we have?
        if ftbl != sys.stdout:
            p=Path(ftbl)
            if p.suffix == ".ftbl":
                ftbl=p
            elif ftbl == "-":
                ftbl=sys.stdout
            else:
                ftbl=p.parent/(p.name+".ftbl")
        # check if we can overwrite
        if not force and ftbl != sys.stdout:
            if ftbl.is_file() and ftbl.stat().st_size > 0:
                with ftbl.open(mode="rb") as fc:
                     if scre[:23] != fc.read(23).decode():
                         werr(f"cannot overwrite '{fc.name}' as not created by this script. Use '--force' to go through.")
        ftbl.parent.mkdir(parents=True, exist_ok=True)
        # compile ftbl dict 'dsec'
        # make prl relative to main ftbl
        if ftbl == sys.stdout:
            wd=Path().resolve()
        else:
            wd=Path(ftbl).resolve().parent
        if case_i:
            #import pdb; pdb.set_trace()
            dsec,dclen,df_kin=compile(rmtf, cmd, case_i)
            p=Path(ftbl).resolve()
            p.parent.mkdir(parents=True, exist_ok=True)
            if p.is_file() and not force:
                with p.open(mode="rb") as fc:
                    if scre[:23] != fc.read(23).decode():
                         werr(f"cannot overwrite '{fc.name}' as not created by this script. Use '--force' to go through.")
            if len(df_kin) > 0:
                fkin=p.with_suffix(".ikin")
                with fkin.open("w") as fc:
                    fc.write(scre+f" at {dtstamp()}\nrow_col")
                    #import pdb; pdb.set_trace()
                    #df_kin.reindex(index=sorted(df_kin.index, key=natural_sort_key))
                    df_kin.loc[sorted(df_kin.index, key=natural_sort_key), :].to_csv(fc, sep="\t")
                #print(str(fkin))
                for i,v in enumerate(dsec["opt"]):
                    if v.startswith("\tfile_labcin\t"):
                        dsec["opt"][i]=v.replace("\tfile_labcin\t", "\t//file_labcin\t")
                dsec["opt"].append("\tfile_labcin\t"+str(fkin.relative_to(p.parent)))
            # prl
            prl_li=[]
            for d in prl:
                dsec_prl,dclen,df_kin_prl=compile(d, cmd, case_i, clen=dclen)
                # output ftbl
                p=Path(d["ftbl"]).resolve()
                p.parent.mkdir(parents=True, exist_ok=True)
                if p.is_file() and not force:
                    with p.open() as fc:
                        if scre[:23] != fc.read(23):
                            werr(f"cannot overwrite '{fc.name}' as not created by this script. Use '--force' to go through.")
                if len(df_kin_prl) > 0:
                    fkin=p.with_suffix(".ikin")
                    with fkin.open("w") as fc:
                        fc.write(scre+f" at {dtstamp()}\nrow_col")
                        df_kin_prl.loc[sorted(df_kin_prl.index, key=natural_sort_key), :].to_csv(fc, sep="\t")
                    for i,v in enumerate(dsec_prl["opt"]):
                        if v.startswith("\tfile_labcin\t"):
                            dsec_prl["opt"][i]=v.replace("\tfile_labcin\t", "\t//file_labcin\t")
                    out=p.open("w")
                    out.write(scre+f" at {dtstamp()}\n")
                    dsec_prl["opt"].append("\tfile_labcin\t"+str(fkin.relative_to(p.parent)))
                dsec2out(dsec_prl, out)
                prl_li.append(str(p.relative_to(wd)))
                out.close()
            if prl_li:
                for i,v in enumerate(dsec["opt"]):
                    if v.startswith("\tprl_exp\t"):
                        dsec["opt"][i]=v.replace("\tprl_exp\t", "\t//prl_exp\t")
                dsec["opt"].append("\tprl_exp\t"+"; ".join(prl_li))
        else:
            dsec,dclen=compile(rmtf, cmd)
            # prl
            prl_li=[]
            for d in prl:
                dsec_prl,dclen=compile(d, cmd, clen=dclen)
                # output ftbl
                p=Path(d["ftbl"]).resolve()
                p.parent.mkdir(parents=True, exist_ok=True)
                if p.is_file() and not force:
                    with p.open() as fc:
                        if scre[:23] != fc.read(23):
                             werr(f"cannot overwrite '{fc.name}' as not created by this script. Use '--force' to go through.")
                out=p.open("w")
                out.write(scre+f" at {dtstamp()}\n")
                dsec2out(dsec_prl, out)
                prl_li.append(str(p.relative_to(wd).with_suffix("")))
                out.close()
            if prl_li:
                for i,v in enumerate(dsec["opt"]):
                    if v.startswith("\tprl_exp\t"):
                        dsec["opt"][i]=v.replace("\tprl_exp\t", "\t//prl_exp\t")
                dsec["opt"].append("\tprl_exp\t"+"; ".join(prl_li))
        # output ftbl
        out=ftbl.open("w") if type(ftbl) == type(Path()) else ftbl
        out.write(scre+f" at {dtstamp()}\n")
        dsec2out(dsec, out)
        out.close()
        if case_i and type(ftbl) == type(Path()):
            from ftbl2labcin import main as renum
            renum([str(ftbl)])
        if res_ftbl is not None:
            res_ftbl.append(out.name)
        out.close()
    return 0
if __name__ == "__main__" or __name__ == "influx_si.cli":
    main()
