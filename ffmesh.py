#!/usr/bin/env python
# Generate a mesh of starting values for free fluxes
# and check the feasibility (inequality constraints) at each node
# If asked a khi2 valuea are calculated at feasible nodes.

# 2008-12-10 sokol
# Copyright 2008, INRA

# imports
import sys;
import os;
import time;

sys.path.append('/home/sokol/dev/python');
sys.path.append('/home/sokol/insa/sysbio/dev/ftbl2sys');
from tools_ssg import *;
#from C13_ftbl import *;
import C13_ftbl;

# take arguments

#<--skip in interactive session
me=sys.argv[0];
def usage():
    sys.stderr.write("usage: "+me+" mesh.kvh");
if len(sys.argv) < 2:
    usage();
    exit(1);
# set some python constants
ffmesh=sys.argv[1];
#-->
#ffmesh="ffmesh_simple.kvh"

# read mesh parameters from a kvh file
fp=open(ffmesh, "r");
meshpar=kvh2dict(fp);
fp.close();
#aff("m", meshpar);##
wd=os.path.dirname(ffmesh);

if meshpar.get("ff",{}).get("mesh_parameters")=="start\tend\tn":
    # start-end-n type of mesh
    # write R and fortran code which runs through this mesh
    import ftbl2code;
    ftbl2code.DEBUG=DEBUG;
    #org="ex3";
    #org="PPP_exact";
    #DEBUG=True;
    if DEBUG:
        import pdb;
    n_ftbl=(meshpar["ftbl"]["name"] if os.path.isabs(meshpar["ftbl"]["name"])
        else os.path.join(wd,meshpar["ftbl"]["name"]));
    org=ffmesh[:-4];
    n_opt=org+".R";
    n_fort=org+".f";
    n_kvh=org+"_res.kvh";
    f_ftbl=open(n_ftbl, "r");
    try:
        os.system("chmod u+w '%s'"%n_ftbl);
        os.system("chmod u+w '%s'"%n_fort);
    except:
        pass;
    f=open(n_opt, "w");
    ff=open(n_fort, "w");

    # parse ftbl
    ftbl=C13_ftbl.ftbl_parse(f_ftbl);
    f_ftbl.close();

    # analyse network
    # reload(C13_ftbl);

    netan=C13_ftbl.ftbl_netan(ftbl);

    # write initialization part of R code
    ftbl2code.netan2Rinit(netan, org, f, ff);
    f.write("""# open connection to store results
cnct=file("%(n_kvh)s", "w");
descr=c(date="%(date)s", generator="%(generator)s");
obj2kvh(descr, "feasible free fluxes", cnct);
feas_tot=0;
"""%{
        "n_kvh": n_kvh,
        "date": time.strftime("%d-%m-%Y"),
        "generator": join(" ", sys.argv),
    });

    # nested loops on the mesh
    indent=0;
    for (fl,par) in meshpar["ff"]["net"].iteritems():
        f.write(indent*"   "+
        'for (%(fl)s.net in seq(from=%(start)s, to=%(end)s, len=%(n)s)) {\n'%
        {
            "fl": fl,
            "start": par["start"],
            "end": par["end"],
            "n": par["n"],
        });
        indent+=1;
    for (fl,par) in meshpar["ff"]["xch"].iteritems():
        f.write(indent*"   "+
        'for (%(fl)s.xch in seq(from=%(start)s, to=%(end)s, len=%(n)s)) {\n'%
        {
            "fl": fl,
            "start": par["start"],
            "end": par["end"],
            "n": par["n"],
        });
        indent+=1;
    # R: check feasibility of flux set
    f.write(indent*"   "+"param[1:no_ff]=c(%(fls)s);\n"%{
       "fls": ",".join([fl+".net" for fl in meshpar["ff"]["net"].keys()]+
       [fl+".xch" for fl in meshpar["ff"]["xch"].keys()]),
    });
    f.write(indent*"   "+"names(param)=nm_par;\n");
    f.write(indent*"   "+"if (all(ui%*%param-ci>=0)) {\n");
    f.write((indent+1)*"   "+"# this free flux vector is feasible => save it\n");
    f.write((indent+1)*"   "+"feas_tot=feas_tot+1;\n");
    f.write((indent+1)*"   "+'obj2kvh(param[1:no_ff], paste("ff_", feas_tot, sep=""), cnct);\n');
    f.write(indent*"   "+"}\n");
    # close loops
    for (fl,par) in meshpar["ff"]["net"]:
        indent-=1;
        f.write(indent*"   "+"}\n");
    for (fl,par) in meshpar["ff"]["xch"]:
        indent-=1;
        f.write(indent*"   "+"}\n");
else:
    raise NameError("Unknown mesh_parameters value '%s'"%meshpar.get("ff").get("mesh_parameters"));
    exit(1);
