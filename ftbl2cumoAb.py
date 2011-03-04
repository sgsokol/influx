#!/usr/bin/env python
"""Transform .ftbl file to human readable matrix A and right hand part b
in full (default) or reduced (option -r) system A*x=b.

usage: ./ftbl2cumoAb.py [-h|--help|-r[--DEBUG] network.ftbl > network_symeq.txt
"""
# 2008-10-14 sokol

import sys, os;
import time;
import getopt;
import numpy as np;

sys.path.append('/home/sokol/dev/python');
from tools_ssg import *;
import C13_ftbl;

def usage():
    print(__doc__);
def get_net(r, dfc):
    """Get net flux value for a reaction r from a dictionary dfc w/ keys like "d.n.v1" """
    for s in ("d", "f", "c"):
       res=dfc.get(s+".n."+r, None);
       if res != None:
           break;
    return res;
# get arguments
try:
    opts,args=getopt.getopt(sys.argv[1:], "hr", ["help", "DEBUG"]);
except getopt.GetoptError, err:
    print str(err);
    usage();
    sys.exit(1);
DEBUG=False;
reduced=False;
for o,a in opts:
    if o in ("-h", "--help"):
        usage();
        sys.exit(0);
    elif o=="--DEBUG":
        DEBUG=True;
    elif o=="-r":
        reduced=True;
    else:
        assert False, "unhandled option";
if not args:
    sys.stderr("Expecting ftbl file name\n");
    usage();
fftbl=args[0] if len(args) else "";
if fftbl and fftbl[-5:] != ".ftbl":
    fftbl+=".ftbl";
if fftbl and not os.path.exists(fftbl):
    sys.stderr.write(me+": file '"+fftbl+"' does not exist.\n");
    sys.exit(1);
fftbl=open(fftbl, "r") if fftbl else sys.stdin;

# parse ftbl from stdin
ftbl=C13_ftbl.ftbl_parse(fftbl);

# analyse network
netan=C13_ftbl.ftbl_netan(ftbl);
f=sys.stdout;
f.write("# This is automatically generated text. Don't edit.\n");
f.write("# Generated by "+sys.argv[0]+" at "+time.ctime()+".\n");
f.write("""
# Copyright Metasys, INSA/INRA UMR 792, Toulouse, France.
# assign A and b weight by weight;
""");

# prepare numpy Afl and bfl for dependent flux solving
Afl=np.array(netan["Afl"]);
# prepare dictionary with dependent, free and constraint fluxe values
# (dependent are added later)
dfc_val=dict(("f.n."+f, v) for (f,v) in netan["flux_free"]["net"].iteritems());
dfc_val.update(("f.x."+f, v) for (f,v) in netan["flux_free"]["xch"].iteritems());
dfc_val.update(("c.n."+f, v) for (f,v) in netan["flux_constr"]["net"].iteritems());
dfc_val.update(("c.x."+f, v) for (f,v) in netan["flux_constr"]["xch"].iteritems());
bfl=np.array( [(sum(dfc_val.get(f, 1.)*v
    for (f,v) in row.iteritems()) if row else 0.) for row in netan["bfl"]] );
# solve Afl*d=bfl;
try:
    d=np.linalg.solve(Afl, bfl);
    dfc_val.update(("d.n."+f, d[i]) for (i,f) in enumerate(netan["vflux"]["net"]));
    dfc_val.update(("d.x."+f, d[i+len(netan["vflux"]["net"])]) for (i,f) in enumerate(netan["vflux"]["xch"]));
except Exception as err:
    sys.stderr.write("dependent net fluxes="+str(netan["vflux"]["net"])+"\n");
    sys.stderr.write("dependent xch fluxes="+str(netan["vflux"]["xch"])+"\n");
    sys.stderr.write("Afl="+str(netan["Afl"])+"\n");
    sys.stderr.write("bfl="+str(netan["bfl"])+"\n");
    sys.stderr.write(str(err)+"\n");

# print flux values
f.write("""
Flux values
%(f)s
""" % {
"f": join("\n", sorted(f+"="+str(v) for (f,v) in dfc_val.iteritems()))
});

# stoichiometric equations
f.write("""
Stoichiometric equations:
Metab:<tab>sum influxes=sum outfluxes
""");
for metab in sorted(netan["sto_m_r"].keys()):
    f.write("%s:\t"%metab);
    lr=netan["sto_m_r"][metab];
    f.write("%(in)s=%(out)s\n"%{
    "in": "+".join(lr["right"]) or "<entering flux>",
    "out": "+".join(lr["left"]) or "<exiting flux>",
    });
    if lr["right"] and lr["left"]:
        #print "left | right=", str(lr["right"])+"|"+str(lr["left"]);##
        infl=join("+", (get_net(r, dfc_val) for r in lr["right"]));
        outfl=join("+", (get_net(r, dfc_val) for r in lr["left"]));
        try:
           verdict="OK" if abs(eval(infl)-eval(outfl)) < 1.e-14 else "bad";
        except:
           verdict="No decision";
        f.write("\t%(in)s=%(out)s\t%(verdict)s\n"%{
            "in": infl,
            "out": outfl,
            "verdict": verdict,
        });

f.write("""
Full flux equations:
net fluxes\t|exchange fluxes\t=b\n
""");
nb_fnet=len(netan["vflux"]["net"]);
for (ir,row) in enumerate(netan["Afl"]):
    f.write("%(fnet)s\t|%(fxch)s\t=%(b)s\n"%{
        "fnet": "\t".join(("" if coef==0 else ssign(coef)+
            "d.n."+netan["vflux"]["net"][i]) for (i,coef) in enumerate(row)
            if i < nb_fnet),
        "fxch": "\t".join(("" if coef==0 else ssign(coef)+
            "d.x."+netan["vflux"]["xch"][i-nb_fnet]) for (i,coef) in enumerate(row)
            if i >= nb_fnet),
        "b": join(" + ", (str(coef)+("*" if (coef and fl) else " ")+str(fl)
            for (fl,coef) in netan["bfl"][ir].iteritems()), a="0"),
    });


# cumomer balance equations
#pdb.set_trace();
measures={"label": {}, "mass": {}, "peak": {}};
o_meas=measures.keys(); # ordered measure types
# calculate measure matrices (mapping cumomers to observations)
f.write("""
Measures:
""");
for meas in o_meas:
    measures[meas]=eval("C13_ftbl.%s_meas2matrix_vec_dev(netan)"%meas);
    # measure vector
    #aff(meas, measures[meas]);##
    f.write(meas+": \n"+join("\n", (row["scale"]+" "+str(i)+": "+
        join(", ", row["coefs"].iteritems()) for (i,row) in
        enumerate(measures[meas]["mat"])))+"\n");

if reduced:
    Ab=C13_ftbl.rcumo_sys(netan);
    f.write("\n# reduced to measurable cumomers system\n");
else:
    Ab=netan["cumo_sys"];
    f.write("\n# full (not reduced to measures) system\n");
#aff("A3", Ab["A"][2]);##
#aff("\nb3", Ab["b"][2]);##
for (w,A) in enumerate(Ab["A"]):
    f.write("\n# weight %d\n"%(w+1));
    for r_cumo in sorted(A.keys()):
        # output term
        f.write("%(c)s*(%(f)s)" % ({
                "c": r_cumo,
                "f": join("+", A[r_cumo][r_cumo]),
            }));
        # input terms
        term=join("+",("%(c)s*(%(f)s)" % ({
              "c": c_cumo,
              "f": join("+", A[r_cumo][c_cumo]),
            }) for c_cumo in sorted(A[r_cumo].keys()) if c_cumo != r_cumo));
        if term:
            f.write("-("+term+")");
        # rhs
        if r_cumo not in Ab["b"][w]:
            f.write("=0.\n");
            continue;
        f.write("="+join("+", ("%(f)s*(%(pc)s)" % ({
                "f": fl,
                "pc": join("+", (join("*", Ab["b"][w][r_cumo][fl][i])
                    for i in Ab["b"][w][r_cumo][fl]))
            }) for fl in Ab["b"][w][r_cumo]))+"\n");

#f.close(); # no need to close stdout
