#!/usr/bin/env python3
"""Transform .ftbl file to human readable matrix A and right hand part b
in full cumomer (default), reduced cumomer (option -r) or emu (option --emu) system A*x=b.
Option '-i' should be used for instationary cases.

usage: ./ftbl2cumoAb.py [-h|--help|-r network.ftbl | --prefix PREFIX | --mtf MTF [> network.sys]
If neither output redirection (>) nor a pipe (|) are not used, the output
file name is silently formed from the input one by replacing the suffix .ftbl with .sys. If exists, the output file is rewritten without warning.
If you wish an output on the screen, use

$ ./ftbl2cumoAb.py network.ftbl | cat

"""
# 2008-10-14 sokol: initial release
# 2012-06-26 sokol: added --emu option
# 2014-03-04 sokol: fixed growth fluxes

import sys, os, stat
import math
import time
import getopt
import numpy as np
from fractions import Fraction as Frac

import influx_si
from tools_ssg import *
import C13_ftbl
import txt2ftbl

#import pdb

me=os.path.basename(sys.argv[0]);
def usage():
    print(__doc__)
def seli(l, i):
    return type(l)(l[ii] for ii in i)
def get_net(r, dfc):
    """Get net flux value for a reaction r from a dictionary dfc w/ keys like "d.n.v1" """
    for s in ("d", "f", "c", "g"):
        res=dfc.get(s+".n."+r, None)
        if res != None:
            break
    return res
def net2type(r, dfc):
    """Return net flux type (one of "d", "f", "c", "g") for a reaction r from a dictionary dfc w/ keys like "d.n.v1". Return None if r is not found."""
    for s in ("d", "f", "c", "g"):
        if s+".n."+r in dfc:
            return s
    return None
# get arguments
try:
    opts,args=getopt.getopt(sys.argv[1:], "hri", ["help", "clownr", "emu", "prefix=", "flux=", "mtf="])
except getopt.GetoptError as err:
    print(str(err))
    usage()
    sys.exit(1)
reduced=False
emu=False
clownr=False
case_i=False
has_flux=False
li_ftbl=[]
mtf_opts=[]
for o,a in opts:
    if o in ("-h", "--help"):
        usage()
        sys.exit(0)
    elif o=="-r":
        reduced=True
    elif o=="--emu":
        emu=True
    elif o=="--clownr":
        clownr=True
    elif o=="-i":
        case_i=True
        mtf_opts += ["--inst"]
    elif o=="--flux":
        has_flux=True
        fv=txt2ftbl.tsv2df(a)
        fv.isetitem(1, fv.iloc[:,1].astype(float))
    elif o=="--prefix" or o=="--mtf":
        mtf_opts += [o, a]
    else:
        assert False, "unhandled option '"+o+"'"

if "--mtf" in mtf_opts or "--prefix" in mtf_opts:
    txt2ftbl.main(mtf_opts, li_ftbl)

if not args and not li_ftbl:
    sys.stderr.write("Error: expecting ftbl file name or --prefix/--mtf options\n")
    usage()
    sys.exit(0)
fullsys=not reduced
C13_ftbl.clownr=clownr
C13_ftbl.ffguess=True
fftbl=args[0] if args else li_ftbl[0]

if fftbl and fftbl[-5:] != ".ftbl":
    fftbl+=".ftbl"
if fftbl and not os.path.exists(fftbl):
    sys.stderr.write(me+": file '"+fftbl+"' does not exist.\n")
    sys.exit(1)

# what kind of output we have?
mode=os.fstat(1).st_mode
f=sys.stdout if stat.S_ISFIFO(mode) or stat.S_ISREG(mode) else  open(fftbl[:-4]+"sys", "w")

# parse ftbl
ftbl=C13_ftbl.ftbl_parse(fftbl)

# analyse network
netan=dict()
C13_ftbl.ftbl_netan(ftbl, netan, emu, fullsys, case_i)
f.write("# This is automatically generated text. Don't edit.\n")
f.write("# Generated by "+sys.argv[0]+" at "+time.ctime()+".\n")

# assign A and b weight by weight
# prepare numpy Afl and bfl for dependent flux solving
Afl=np.array(netan["Afl"])
# prepare dictionary with dependent, free, constraint and growth flux values
# (dependent are added later)
dfc_val=dict(("f.n."+f, v) for (f,v) in netan["flux_free"]["net"].items())
dfc_val.update(("f.x."+f, v) for (f,v) in netan["flux_free"]["xch"].items())
dfc_val.update(("c.n."+f, v) for (f,v) in netan["flux_constr"]["net"].items())
dfc_val.update(("c.x."+f, v) for (f,v) in netan["flux_constr"]["xch"].items())
dfc_val.update(("g.n."+f, v) for (f,v) in netan["flux_vgrowth"]["net"].items())

# update dfc by values from file
if has_flux:
    dset=set(netan["vflux"]["net"])
    for i,row in fv.iterrows():
        nm=row[0]
        t=net2type(nm, dfc_val)
        if t is None:
            if nm in dset:
                t="d"
            else:
                continue # ignore names not defined in ftbl
        fnm=t+".n."+nm
        dfc_val[fnm]=row[1]
bfl=np.array( [(sum(dfc_val.get(f, 1.)*v
    for (f,v) in row.items()) if row else 0.) for row in netan["bfl"]] )
# test for linear dependence of rows in Afl
# get vector of col names in bfl
cnm=sorted(set(f for r in netan["bfl"] for f in r.keys() if f))
ncbfl=len(cnm)
i2cnm=dict(enumerate(cnm))
# make numpy matrix from bfl
mbfl=np.array([[r.get(i2cnm.get(i, None), 0.) for i in range(ncbfl)] for r in netan["bfl"]])
afull=np.concatenate((Afl, -mbfl), axis=1)
qf,rf=np.linalg.qr(afull.T)
rd=np.diag(rf)
ikeep=np.where(np.abs(rd) >= 1.e-10)[0]
#pdb.set_trace()
if len(ikeep) < Afl.shape[0]:
    sys.stderr.write(f"Warning: found {Afl.shape[0] - len(ikeep)} linearly dependent rows in stoechiometric matrix.\nThe following rows will be ignored:\n\t"+"\n\t".join(np.array(netan["vrowAfl"])[np.abs(rd) < 1.e-10])+"\n")
    Afl=Afl[ikeep,:]
    bfl=bfl[ikeep]
    netan["Afl"]=seli(netan["Afl"], ikeep)
    netan["vrowAfl"]=seli(netan["vrowAfl"], ikeep)
    netan["bfl"]=seli(netan["bfl"], ikeep)
# solve Afl*d=bfl
#pdb.set_trace()
# exclude linearly redundant rows in Afull

invAfl=None
if has_flux:
    d_avail=True
else:
    try:
        d=np.linalg.solve(Afl, bfl)
        dfc_val.update(("d.n."+f, d[i]) for (i,f) in enumerate(netan["vflux"]["net"]))
        dfc_val.update(("d.x."+f, d[i+len(netan["vflux"]["net"])]) for (i,f) in enumerate(netan["vflux"]["xch"]))
        invAfl=np.matrix(Afl).I
        d_avail=True
    except Exception as err:
        sys.stderr.write("Error: Afl is singular or is not square\n")
        sys.stderr.write("nrow x ncol = %d x %d\n"%Afl.shape)
        sys.stderr.write("dependent net fluxes="+str(netan["vflux"]["net"])+"\n")
        sys.stderr.write("dependent xch fluxes="+str(netan["vflux"]["xch"])+"\n")
        sys.stderr.write("Afl="+str(netan["Afl"])+"\n")
        sys.stderr.write("bfl="+str(netan["bfl"])+"\n")
        sys.stderr.write(str(err)+"\n")
        d_avail=False

# print flux values
f.write("""
Flux values
%(f)s
""" % {
"f": join("\n", sorted(f+"="+str(v) for (f,v) in dfc_val.items()))
})

# stoichiometric equations
if d_avail:
    f.write("""
    Stoichiometric equations:
    Metab:<tab>sum influxes=sum outfluxes
    """)
    for metab,lr in sorted(netan["sto_m_r"].items()):
        f.write("%s:\t"%metab)
        f.write("%(in)s = %(out)s\n"%{
        "in": " + ".join((str(co)+"*" if co != 1. else "")+r+f"({net2type(r, dfc_val)})" for r,co in lr["right"]) or "<entering flux>",
        "out": " + ".join((str(co)+"*" if co != 1. else "")+r+f"({net2type(r, dfc_val)})" for r,co in lr["left"]) or "<exiting flux>",
        })
        if lr["right"] and lr["left"]:
            #print "left | right=", str(lr["right"])+"|"+str(lr["left"]);##
            #pdb.set_trace()
            fl=[join("", ((" + " if i>0 and math.copysign(1, get_net(r, dfc_val)) > 0 else " - " if math.copysign(1, co*get_net(r, dfc_val)) < 0 else " ")+(str(co)+"*" if co != 1. else "")+str(abs(get_net(r, dfc_val))) for i,(r,co) in enumerate(side))) for side in (lr["right"], lr["left"])]
            try:
               dif=abs(eval(fl[0])-eval(fl[1]))
               verdict="OK" if dif < 1.e-9 else "bad ("+str(dif)+")"
            except:
               verdict="No check status"
            f.write("\t%(in)s =%(out)s\t%(verdict)s\n"%{
                "in": fl[0],
                "out": fl[1],
                "verdict": verdict,
            })
else:
    f.write("\n***Stoichiometric equations cannot be checked as the stoichiometric matrix is not invertible.\n")
#pdb.set_trace()
f.write("""
Full flux equations:
metab:net fluxes\t|exchange fluxes\t=b\n
""")
nb_fnet=len(netan["vflux"]["net"])
for (ir,row) in enumerate(netan["Afl"]):
    f.write("%(metab)s:%(fnet)s\t|%(fxch)s\t=%(b)s\n"%{
        "metab": netan["vrowAfl"][ir],
        "fnet": "\t".join(("" if coef==0 else ssign(coef)+(str(abs(coef)) if abs(coef) !=1. else "")+
            "d.n."+netan["vflux"]["net"][i]) for (i,coef) in enumerate(row)
            if i < nb_fnet),
        "fxch": "\t".join(("" if coef==0 else ssign(coef)+(str(abs(coef)) if abs(coef) !=1. else "")+
            "d.x."+netan["vflux"]["xch"][i-nb_fnet]) for (i,coef) in enumerate(row)
            if i >= nb_fnet),
        "b": join(" + ", ((str(coef) if abs(coef) != 1 else "" if coef == 1 else "-") +("*" if (abs(coef) != 1 and fl) else "")+str(fl)
            for (fl,coef) in netan["bfl"][ir].items()), a="0"),
    })

if not invAfl is None:
    # show formulas for dependent fluxes using inverted Afl
    # free, constrained and constant value name to index translator
    nfn=len(netan["vflux_free"]["net"])
    nfx=len(netan["vflux_free"]["xch"])
    nf=nfn+nfx
    ncn=len(netan["vflux_constr"]["net"])
    ncx=len(netan["vflux_constr"]["xch"])
    nc=ncn+ncx
    ngn=len(netan["vflux_growth"]["net"])
    ng=ngn
    
    # constant part
    fcv2i={"": 0}
    
    # free part
    fcv2i.update(("f.n."+k, v+1) for k,v in netan["vflux_free"]["net2i"].items())
    fcv2i.update(("f.x."+k, v+nfn+1) for k,v in netan["vflux_free"]["xch2i"].items())
    
    # constrained part
    fcv2i.update(("c.n."+k, v+nf+1) for k,v in netan["vflux_constr"]["net2i"].items())
    fcv2i.update(("c.x."+k, v+nf+ncn+1) for k,v in netan["vflux_constr"]["xch2i"].items())

    # growth part
    fcv2i.update(("g.n."+k, v+nf+nc+1) for k,v in netan["vflux_growth"]["net2i"].items())
    
    # inverse: from index to name
    fl,i=list(zip(*iter(fcv2i.items())))
    i2fcv=np.array(fl)
    i2fcv[np.array(i)]=fl
    
    i,j,v=[ np.array(i) for i in zip(*((i,fcv2i[f],v) for (i,row) in enumerate(netan["bfl"]) if row for (f,v) in row.items())) ]
    
    f2bfl=np.zeros((Afl.shape[0], nf+nc+ng+1))
    f2bfl[i,j]=v
    fcv2dep=invAfl*np.matrix(f2bfl)
    f.write("""
Dependent fluxes as functions of free and constrained fluxes:
dep.flux=f(free.flux, constr.flux)
----------------------------------
""")
    ndn=len(netan["vflux"]["net"])
    ndx=len(netan["vflux"]["xch"])
    nd=ndn+ndx
    assert fcv2dep.shape[0] == nd, "Bad dimensions dependent fluxes as function of free, constrained and constant values.\nWe should have %d rows, instead we got %d."%(nd, fcv2dep.shape[0])
    ira=np.arange(fcv2dep.shape[1])
    for (ir,row) in enumerate(fcv2dep.A):
        nz=abs(row) >= 1.e-10
        #pdb.set_trace()
        formula=join("", ( (("+" if v > 0 else "-") if abs(abs(v)-1.) < 1.e-10 and f != "" else ("+" if v > 0 else "") + "%s"% Frac.from_float(v).limit_denominator(10000)) + ("*" if f != "" and abs(abs(v)-1) > 1.e-10 else "") + f for (v,f) in zip(row[nz],i2fcv[nz])))
        formula="0" if not formula else formula if formula[0] != "+" else formula[1:]
        f.write("%(dep)s=%(formula)s\n"%{
            "dep": "d.n."+netan["vflux"]["net"][ir] if ir < ndn else "d.x."+netan["vflux"]["xch"][ir-ndn],
            "formula": formula,
        })

# cumomer or emu balance equations
#pdb.set_trace()
measures={"label": {}, "mass": {}, "peak": {}}
o_meas=list(measures.keys()); # ordered measure types
# calculate measure matrices (mapping cumomers to observations)
f.write("""
Measurements:
""")
for meas in o_meas:
    measures[meas]=eval("C13_ftbl.%s_meas2matrix_vec_dev(netan)"%meas)
    # measure vector
    #aff(meas, measures[meas]);##
    for iexp in range(len(measures[meas])):
        f.write("%s:\n"%netan["exp_names"][iexp])
        f.write(meas+":\n"+join("\n", (row["scale"]+" "+str(i)+": "+
            join(", ", iter(row["coefs"].items())) for (i,row) in
            enumerate(measures[meas][iexp]["mat"]) if row["coefs"] is not None))+"\n")

Ab=C13_ftbl.rcumo_sys(netan, emu)
vcumo=netan.get("vrcumo")
if emu:
    f.write("\n# EMU fragment system\n")
elif reduced:
    f.write("\n# reduced to measurable cumomers system\n")
else:
    Ab=netan["cumo_sys"]
    vcumo=netan["vcumo"]
    f.write("\n# full (not reduced to measurable) cumomer system\n")
#aff("A3", Ab["A"][2]);##
#aff("\nb3", Ab["b"][2]);##
#pdb.set_trace()
for (w,A) in enumerate(Ab["A"]):
    f.write("\n# weight %d\n"%(w+1))
    #for r_cumo in sorted(A.keys()):
    for r_cumo in vcumo[w]:
        # output term
        f.write("%(c)s*(%(f)s)" % ({
                "c": r_cumo,
                "f": join("+", A[r_cumo][r_cumo]),
            }))
        # input terms
        term=join("+",("%(c)s*(%(f)s)" % ({
              "c": c_cumo,
              "f": join("+", A[r_cumo][c_cumo]),
            }) for c_cumo in vcumo[w] if c_cumo in A[r_cumo] and c_cumo != r_cumo))
            #}) for c_cumo in sorted(A[r_cumo].keys()) if c_cumo != r_cumo))
        if term:
            f.write("-("+term+")")
        # rhs
        if r_cumo not in Ab["b"][w]:
            f.write("=0.\n")
            continue
        f.write("="+join("+", ("%(f)s*(%(pc)s)" % ({
                "f": fl,
                "pc": join("+", (join("*", Ab["b"][w][r_cumo][fl][i])
                    for i in Ab["b"][w][r_cumo][fl]))
            }) for fl in Ab["b"][w][r_cumo]))+"\n")

#f.close(); # no need to close stdout
#f.write(str(netan["emu_input"])+"\n")
