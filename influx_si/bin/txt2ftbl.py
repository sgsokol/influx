#!/usr/bin/env python3

"""
read a .txt file from a parameter and translate to NETWORK and FLUXES section of .ftbl file.
The generated FTBL file can be then edited by hand to be added
other sections (MEASUREMENTS and so on)
Comments starting with '###' in txt file separate pathways which are numbered
as well as reactions in them. If no explicite name "reac: " is given at
the begining of the line, ractions in ftbl will be named as
"rX.Y" where X is pathhway number and Y is reaction number in the
pathway.
Symbols "+", "(", ")" and ":" are not allowed in metabolite neither reaction names
Empty lines are ignored.

usage: txt2ftbl.py [-h|--help] mynetwork.txt [> mynetwork.ftbl]

OPTIONS
-h, --help print this message and exit

:param: mynetwork the base of an txt file (mynetwork.txt)

:returns: mynetwork.ftbl -- file of the network definition in FTBL format

Copyright 2015, INRA, France
Author: Serguei Sokol (sokol at insa-toulouse dot fr)
License: Gnu Public License (GPL) v3 http://www.gnu.org/licenses/gpl.html
"""

import re
def txt_parse(fname, re_metab=re.compile(r"(?:(?P<coef>[\d\.]*)\s+)?(?:(?P<metab>[^() \t\r]+)\s*)(?:\(\s*(?P<carb>[^()]*)\s*\))?\s*"),
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

Retrun a list with four items:
- a list of carbon exchange reactions
- a list of non carbon echanging reactions
- a list of equalities net and xch
- a list of flux tuples (reac, rev, imposed_sens)
- a list of two lists: left and right metabolites [(metab, clen),]

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
    comment=""
    
    open_here=False
    if isstr(fname):
        open_here=True
        fc=open(fname, "r")
    else:
    	fc=fname
    m_left={} # metab sources
    m_right={} # metab products
    sto={} # stoechiometric matrix dictionary {flux:{metab: coef}}}
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
        li=l.split(":", 1)
        if len(li)==2:
            nm_reac=li[0].strip()
            reac=li[1].strip()
        else:
            nm_reac="r"+str(ipath)+"."+str(ireac)
            reac=li[0].strip()
        in_out=[side.strip() for side in reac.split("->")]
        if len(in_out) != 2:
            raise Exception("Bad syntax. No reaction detected on the line %d"%iline)
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
        if any(len(t) == 0 for side in lr for c,m,t in side):
            raise Exception("Wrong format for labeling pattern on row %d."%iline)
        if not tr_reac:
            if comment:
                resnotr.append(comment)
                comment=""
            resnotr.append([(nm_reac, rev, imposed_sens)]+lr)
            fluxes.append((nm_reac, rev, imposed_sens, "F"))
            # add a row to stoechiometric matrix
            if nm_reac in sto:
                raise Exception("Reaction '%s' was already met (row: %d)"%(nm_reac, iline))
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
                raise Exception("Reaction '%s' was already met (row: %d)"%(nm_r, iline))
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
            m_left.update((m,clen) for m, clen in es.items() if m not in m_left)
            m_right.update((m,clen) for m, clen in ps.items() if m not in m_right)
            rgr=[] # reaction group (for possible long reactions)
            r=[(nm_r, rev, imposed_sens), [], []] # lace for elementary reaction (no more than 2 metabs on each side)
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
    return [res, resnotr, eqs, ineqs, fluxes, [m_left, m_right], sto]

if __name__ == "__main__" or __name__ == "influx_si.cli":
    import sys
    import os
    import stat
    import getopt
    import re
    import math
    import datetime as dt
    from scipy import linalg
    from numpy import diag

    import influx_si
    from tools_ssg import *

    werr=sys.stderr.write
    #stoe_coeff=Exception("Stoechiometric coefficient is different from 1")

    # get arguments
    me=os.path.basename(sys.argv[0])
    def usage():
        werr(__doc__)
    try:
        opts,args=getopt.getopt(sys.argv[1:], "h", ["help"])
    except getopt.GetoptError as err:
        werr(str(err)+"\n")
        usage()
        sys.exit(1)
    for o,a in opts:
        if o in ("-h", "--help"):
            usage()
            sys.exit()
    if len(args) > 1:
        werr("Expecting exactly one txt file name\n")
        usage()
        exit(1)
    elif len(args) == 0:
        # read on stdin
        path_txt=sys.stdin
        fout=sys.stdout
        base="network"
    else:
        base=args[0]
        if base[-4:]==".txt":
            base=base[:-4]
        path_txt=base+".txt"
        #-->
        
        # what kind of output we have?
        mode=os.fstat(1).st_mode
        fout=sys.stdout if stat.S_ISFIFO(mode) or stat.S_ISREG(mode) else  open(path_txt[:-3]+"ftbl", "w")

        # define where to read and write
        # input file is an argument
        fdir=os.path.dirname(base) or "."
        base=os.path.basename(base)

    # Parse .txt file
    try:
        (netw, notr_netw, eqs, ineqs, fluxes, (m_left, m_right), sto)=txt_parse(path_txt)
    except Exception as inst:
        werr(str(inst)+"\n")
        raise
    # build afl matrix: each row is a balance on an internal metabolite, each column is a flux values
    nb_flux=len(sto)
    nm_flux=sorted(sto.keys())
    m_l=set(m_left.keys())
    m_r=set(m_right.keys())
    m_inp=m_l-m_r
    m_out=m_r-m_l
    nm_met=sorted((m_l|m_r) - m_inp - m_out)
    fl2i=dict((fl, i) for i, fl in enumerate(nm_flux))
    met2i=dict((m, i) for i, m in enumerate(nm_met))
    afl=np.zeros((len(nm_met)+len(eqs[0]), len(nm_flux)))
    for fl, row in sto.items():
        for m, c in row.items():
            if m not in met2i:
                continue
            afl[met2i[m], fl2i[fl]]=c
    # store matrix for R reading
    fmat=open("sto_mat.txt", "w")
    fmat.write("\t".join(["row_col"]+nm_flux)+"\n")
    for ir in range(len(nm_met)):
        fmat.write(nm_met[ir]+"\t"+"\t".join(str(v) for v in afl[ir,:])+"\n")
    fmat.close()
    # add net equations
    for i,eq in enumerate(eqs[0]):
        fls=eq[1].split(" - ")
        afl[len(nm_met)+i, fl2i[fls[0]]]=1.
        afl[len(nm_met)+i, fl2i[fls[1]]]=-1.
    if nb_flux > 1:
        # qr(afl)
        q,r,p=linalg.qr(afl, pivoting=True)
        d=diag(r)
        if d[0] == 0.:
            raise Exception("Stoechiometrix matrix is 0")
        rank=sum(abs(d/d[0]) > 1.e-10)
        # free fluxes are the last in pivots
        ff=set(nm_flux[i] for i in p[rank:])
    else:
        ff=[]
    # write header
    fout.write(
"""PROJECT
	NAME	VERSION	FORMAT	DATE	COMMENT
	%s	1		%s	converted by txt2ftbl.py from '%s'

NETWORK
	FLUX_NAME	EDUCT_1	EDUCT_2	PRODUCT_1	PRODUCT_2
"""%(base, dt.datetime.strftime(dt.datetime.now(), "%Y-%m-%d"), path_txt if path_txt != sys.stdin else "stdin"))

    # write the ftbl content
    for row in netw:
        if isstr(row):
            fout.write("%s"%row)
            continue
        fout.write("\t%s"%row[0][0]) # reac name
        # input metabs
        carbs="\t"
        imetab=0
        for (coef, metab, carb) in row[1]:
            imetab=imetab+1
            fout.write("\t%s%s"%(coef+"*" if coef and coef != "1" else "", metab))
            carbs=carbs+"\t%s"%(("#"+carb) if metab else "")
        if imetab==1:
            # complete by empty metabolite
            fout.write("\t")
            carbs=carbs+"\t"
        # output metabs
        imetab=0
        for (coef, metab, carb) in row[2]:
            imetab=imetab+1
            fout.write("\t%s"%metab)
            carbs=carbs+"\t%s"%(("#"+carb) if metab else "")
        fout.write("\n%s\n"%carbs)
    if notr_netw:
        fout.write(
"""
NOTRACER_NETWORK
	FLUX_NAME	EQUATION
""")
        for row in notr_netw:
            if isstr(row):
                fout.write("%s"%row)
                continue
            fout.write("\t%s\t%s\n"%(row[0][0], " = ".join("+".join((c+"*" if c and c != "1" else "")+m  for c,m,t in side) for side in row[1:3])))
    fout.write("""
EQUALITIES
	NET
		VALUE	FORMULA
""")
    for e in eqs[0]:
        fout.write("\t\t%s\t%s\n"%e)
    fout.write("""	XCH
		VALUE	FORMULA
""")
    for e in eqs[1]:
        fout.write("\t\t%s\t%s\n"%e)
    fout.write("""
FLUXES
	NET
		NAME	FCD	VALUE(F/C)	ED_WEIGHT	LOW(F)	INC(F)	UP(F)
""")
    for f, rev, imposed_sens, fd in fluxes:
        fout.write("\t\t%s\t%s\t0.2E0\n"%(f, "F" if f in ff else "D"))
    fout.write("""	XCH
		NAME	FCD	VALUE(F/C)	ED_WEIGHT	LOW(F)	INC(F)	UP(F)
""")
    for f, rev, imposed_sens, fd in fluxes:
        fout.write("\t\t%s\t%s\n"%(f, "C\t0" if not rev else "%s\t0.01E0"%fd))
    fout.write("""
LABEL_INPUT
	META_NAME	ISOTOPOMER	VALUE
""")
    for m in sorted(m_inp):
        clen=m_left[m]
        if clen == 0:
            continue
        fout.write("\t%s\t#%s\t1\n"%(m, "1"*clen))
    # print footer: empty sections (to complete manually by user)
    fout.write("""
INEQUALITIES
	NET
		VALUE	COMP	FORMULA
""")
    for ine in ineqs[0]:
        fout.write("\t\t%s\t%s\t%s\n"%ine)
    fout.write("""	XCH
		VALUE	COMP	FORMULA
FLUX_MEASUREMENTS
	FLUX_NAME	VALUE	DEVIATION
LABEL_MEASUREMENTS
	META_NAME	CUM_GROUP	VALUE	DEVIATION	CUM_CONSTRAINTS
PEAK_MEASUREMENTS
	META_NAME	PEAK_NO	VALUE_S	VALUE_D-	VALUE_D+	VALUE_DD	VALUE_T	DEVIATION_S	DEVIATION_D-	DEVIATION_D+	DEVIATION_DD/T
MASS_SPECTROMETRY
	META_NAME	FRAGMENT	WEIGHT	VALUE	DEVIATION
OPTIONS
	OPT_NAME	OPT_VALUE
""")
    fout.close()
