#!/usr/bin/python

r"""
read a .txt file from a parameter and translate to NETWORK and FLUXES section of .fbtl file.
The generated FTBL file can be then edited by hand to be added
other sections (MEASUREMENTS and so on)
Comments in txt file separate pathways which are numbered
as well as reactions in them. Ractions in ftbl are then names as
"rX.Y" where X is pathhway number and Y is reaction number in the
pathway.
Symbols "+", "(", ")" and ":" are not allowed in metabolite neither reaction names

usage: txt2ftbl.py [-h|--help] mynetwork.txt [> mynetwork.ftbl]

OPTIONS
-h, --help print this message and exit

:param: mynetwork the base of an txt file (mynetwork.txt)

:returns: mynetwork.ftbl -- file of the network definition in FTBL format

Copyright 2015, INRA, France
Author: Serguei Sokol (sokol at insa-toulouse dot fr)
License: Gnu Public License (GPL) v3 http://www.gnu.org/licenses/gpl.html
"""
def txt_parse(fname):
    """Parse txt file from fname which is in format:
# Glycolysis and OPP pathway
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

reversible reaction are represented by "<->" sign

Retrun a list with four items:
- a list of carbon exchange reactions
- a disctionary with a whole transposed stoechiometric matrix as dict()
- a set of carbon echanging reactions (just their names)
- a set of carbon non echanging reactions.

the first item is a list of:
plain string == just a comments
list == reaction items: input, output: lists of tuples (metab, carb, coeff)
"""
    res=list() # main result
    # transpose of stoechiometric matrix st["reac"]={metab1: coef1, metab2: coef2, ...}
    # input metabs have negative coeffs, output metabs have positive ones.
    st=dict()
    ncmet=set() # collection of non carbon exchanging metabolites (e.g. cofactors)
    cmet=set() # collection of carbon exchanging metabolites
    open_here=False
    if isstr(fname):
        open_here=True
        fc=file(fname, "rb")
    else:
    	fc=fname
    ireac=0
    ipath=0
    iline=0
    for l in fc.readlines():
        iline=iline+1
        l=l.strip()
        if l[0]=="#":
            res.append(l)
            ireac=0
            ipath=ipath+1
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
        st[nm_reac]=dict()
        in_out=reac.split("->")
        rev=False
        if in_out[0][-1:] == "<":
            in_out[0]=in_out[0][:-1].strip()
            rev=True
        if len(in_out)!=2:
            raise Exception("Bad syntax. No reaction detected on the line %d"%iline)
        r=[nm_reac, rev]
        scramble=False
        for i in (0,1):
            mets=list()
            csign=1 if i else -1
            for g in re_metab.findall(in_out[i]):
                coef, m, t=g
                if m=="+":
                    continue
                coef=1 if coef==None or coef == "" else float(coef)
                st[nm_reac][m]=csign*coef+st[reac].get(m, 0.)
                if coef != "1":
                    werr("Warning: stoechiometric coefficients different from 1 are not admitted in FTBL format (row %d)\n"%iline)
                    r=nm_reac+"\t"+l
                    continue # exclude this metabolite from NETWORK section
                    #raise stoe_coeff
                if t:
                    cmet.add(m)
                else:
                    ncmet.add(m)
                # check for scrumbling patterns
                li=t.split("+") if t else []
                if len(li) > 1:
                    # append the same metabolite with two different patterns
                    # the repeated letters in the patterns will require manual curation
                    # in the FTBL
                    if not scramble:
                        res.append("***Warning: Manual curation of the carbon patter needed to avoid repeated letters in different metabolite instances")
                        scramble=True
                    for cl in li:
                        # remove 1/2 from the pattern
                        cl=cl.replace("1/2", "").strip()
                        mets.append((g[0], g[1], cl))
                else:
                    mets.append(g)
            if len(mets)!=1 and len(mets)!=2:
                # warning
                werr("Warning: only 1 or 2 metabolites can appear on each side of reaction in FTBL format (row %d)\n"%iline)
                # send this line just as a comment
                r=nm_reac+"\t"+l
                break
            if type(r)==type([]):
                r.append(mets)
        res.append(r)
    if open_here:
        fc.close()
    return list(res, st, cmet, ncmet)
if __name__ == "__main__":
    import sys
    import os
    import stat
    import getopt
    import re
    import math
    import datetime as dt

    from tools_ssg import *

    werr=sys.stderr.write
    re_metab=re.compile(r"(?:(?P<coef>\d?\.?\d*)\s+)?(?:(?P<metab>[^() \t\r]+)\s*)(?:\(\s*(?P<carb>[^()]*)\s*\))?\s*")
    stoe_coeff=Exception("Stoechiometric coefficient is different from 1")

    # get arguments
    me=os.path.basename(sys.argv[0])
    def usage():
        werr(__doc__)
    try:
        opts,args=getopt.getopt(sys.argv[1:], "h", ["help"])
    except getopt.GetoptError, err:
        werr(str(err)+"\n")
        usage()
        sys.exit(1)
    for o,a in opts:
        if o in ("-h", "--help"):
            usage()
            sys.exit()
    if len(args) != 1:
        werr("Expecting exactly one txt file name\n")
        usage()
        exit(1)
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
        netw, st, cmet, ncmet=txt_parse(path_txt)
    except Exception as inst:
        werr(str(inst)+"\n")
        raise
    # write header
    fout.write(
"""PROJECT
	NAME	VERSION	FORMAT	DATE	COMMENT
	%s	1		%s	converted by txt2ftbl.py from '%s'

NETWORK
	FLUX_NAME	EDUCT_1	EDUCT_2	PRODUCT_1	PRODUCT_2


"""%(base, dt.datetime.strftime(dt.datetime.now(), "%Y-%m-%d"), path_txt))

    # write the ftbl content
    fluxes=list()
    for row in netw:
        if isstr(row):
            fout.write("// %s\n"%row)
            continue
        fout.write("\t%s"%row[0]) # reac name
        # input metabs
        carbs="\t"
        imetab=0
        for (coef, metab, carb) in row[2]:
            imetab=imetab+1
            fout.write("\t%s"%metab)
            carbs=carbs+"\t#%s"%carb
        if imetab==1:
            # complete by empty metabolite
            fout.write("\t")
            carbs=carbs+"\t"
        # output metabs
        imetab=0
        for (coef, metab, carb) in row[3]:
            imetab=imetab+1
            fout.write("\t%s"%metab)
            carbs=carbs+"\t#%s"%carb
        if imetab==1:
            # complete by empty metabolite
            fout.write("\t")
            carbs=carbs+"\t"
        fout.write("\n%s\n\n"%carbs)
        fluxes.append(row[0:2])
    # print footer: other empty sections (to complete manually by user)
    fout.write("""
FLUXES
	NET
		NAME	FCD	VALUE(F/C)	ED_WEIGHT	LOW(F)	INC(F)	UP(F)
""")
    for f, rev in fluxes:
        fout.write("\t\t%s\tD\n"%f)
    fout.write("""
	XCH
		NAME	FCD	VALUE(F/C)	ED_WEIGHT	LOW(F)	INC(F)	UP(F)
""")
    for f, rev in fluxes:
        fout.write("\t\t%s\t%s\n"%(f, "C\t0" if not rev else "F\t0.01"))
    fout.write("""
EQUALITIES
	NET
		VALUE	FORMULA

	XCH
		VALUE	FORMULA


INEQUALITIES
	NET
		VALUE	COMP	FORMULA

	XCH
		VALUE	COMP	FORMULA


LABEL_INPUT
	META_NAME	ISOTOPOMER	VALUE

FLUX_MEASUREMENTS
	FLUX_NAME	VALUE	DEVIATION

LABEL_MEASUREMENTS
	META_NAME	CUM_GROUP	VALUE	DEVIATION	CUM_CONSTRAINTS
											// Example

PEAK_MEASUREMENTS
	META_NAME	PEAK_NO	VALUE_S	VALUE_D-	VALUE_D+	VALUE_DD	VALUE_T	DEVIATION_S	DEVIATION_D-	DEVIATION_D+	DEVIATION_DD/T


MASS_SPECTROMETRY
	META_NAME	FRAGMENT	WEIGHT	VALUE	DEVIATION
OPTIONS
	OPT_NAME	OPT_VALUE
""")
