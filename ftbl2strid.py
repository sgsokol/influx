#!/usr/bin/env python
# analyse a .ftbl file and propose
# a free fluxes choices

# Copyright 2008, INRA.
# 2008-12-04 sokol: initial version

# imports
import sys;
import array;
from rpy import *; # creates r object
import copy;

sys.path.append("/home/sokol/dev/python");
sys.path.append("/home/sokol/insa/sysbio/dev/ftbl2sys");
from tools_ssg import *;
#from C13_ftbl import *;
import C13_ftbl;

#<--skip in interactive session
me=sys.argv[0];
def usage():
    sys.stderr.write("usage: "+me+" organism");
if len(sys.argv) < 2:
    # no argument is given
    usage();
    exit(1);
org=sys.argv[1];
#-->
#org="PPP_sd";
n_ftbl=org+".ftbl";

# parse ftbl file
f_ftbl=open(n_ftbl, "r");
fo=sys.stdout;
ftbl=C13_ftbl.ftbl_parse(f_ftbl);
f_ftbl.close();
netan=C13_ftbl.ftbl_netan(ftbl);

# analyse the rank of flux matrix
# dependent flux counts
nb_fln=len(netan["vflux"]["net"]);
nb_flx=len(netan["vflux"]["xch"]);
nb_fl=nb_fln+nb_flx;
Afl=copy.deepcopy(netan["Afl"]);

rAfl=r.matrix([el for el in valval(Afl)], nrow=len(Afl), byrow=True);
# if flux matrix is underdetermined
# then
# propose choice of free fluxes
qr=r.qr(rAfl);
nb_deff=nb_fl-qr["rank"]; # defficient fluxes
if not nb_deff:
    # there is no defficient fluxes
    exit(0);

# prepare list of candidates for free fluxes
fo.write("""
Flux matrix is under-determined. You have to choose %(nb_deff)d
free fluxes from the following list:\n"""% {
    "nb_deff": nb_deff,
});
# get last pivots for candidates for free fluxes
flnames=([ fl+".net" for fl in netan["vflux"]["net"] ] +
    [ fl+".xch" for fl in netan["vflux"]["xch"] ]);
fldef=[flnames[i-1] for i in qr["pivot"][-nb_deff:]];
fo.write("%s\n"%"\n".join(fldef));
