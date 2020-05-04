#!/usr/bin/env python3

"""
This program gets all fluxes (free, dependent and constrained) as well as free pool values from kvh file
and put them in a ftbl file for free/dependent fluxes and free pools. Values are rounded to 9 digits
(default value but can be changed with -r option).

The first input parameter must point to a valid _res.kvh file or to its prefix.
If the second parameter (ftbl name) is omitted, it is derived from kvh file name.

Names in kvh are like "f.n.FLX" or "f.x.FLX" or "pf:Glc6P"
Names in ftbl are like "FLX" in corresponding NET or XCH section
or like "Glc6P" in METABOLITE_POOLS
If a name in kvh file does not have its equivalent in
ftbl file, it is silently ignored.
Optional comments in the original ftbl after the value of
edited flux or pool are lost.

New content is sent to stdout

usage: ./ff2ftbl.py [-h|--help] [-r 9] f[_res.kvh] [f.ftbl] > new_f.ftbl
or: cat f.kvh | ./ff2ftbl.py - f.ftbl > new_f.ftbl
"""
import sys
import os
import getopt
import re
from math import exp

import influx_si
import kvh
def usage(mes=""):
    sys.stderr.write(os.linesep.join([mes, __doc__]))

try:
    opts,args=getopt.getopt(sys.argv[1:], "hr", ["help"])
except getopt.GetoptError as err:
    #pass
    usage(str(err))
    sys.exit(1)

fullsys=False
nround=9
for o,a in opts:
    if o in ("-h", "--help"):
        usage()
        sys.exit(0)
    elif o == "-r":
        # round option
        nround=int(a)
    else:
        assert False, "unhandled option"
#aff("args", args);##
if len(args) < 1 or len(args) > 2:
    usage("Expecting one or two arguments. Got %d."%len(args))
    exit(1)
fkvh=args[0]
if fkvh == "-":
    if len(args) != 2:
        raise Exception("If kvh is read from stdin, FTBL name must be given in the second parameter")
    # read from standart input
    fkvh=sys.stdin
else:
    if os.path.isfile(fkvh):
        if fkvh[-8:] == "_res.kvh":
            fbase=fkvh[:-8]
        else:
            fbase=None
    else:
        fbase=fkvh
        fkvh+="_res.kvh"
if (len(args) == 2):
    ftbl=args[1]
elif not fbase is None:
    ftbl=fbase+".ftbl"
else:
    raise Exception("When _res.kvh and .ftbl have different basecpart of name, both file names must be provided.")

# get free and dependent fluxes from kvh
ff=kvh.kvh_get_matrix(fkvh, ["linear stats", "net-xch01 fluxes (sorted by name)"])
# get constrained fluxes from kvh
try:
    ff+=kvh.kvh_get_matrix(fkvh, ["constrained net-xch01 fluxes"])
except:
    pass;

# convert strings to floats and round it to nround digits (default 9)
ff=dict((row[0], round(float(row[1]), nround)) for row in ff[1:])

#print("ff=", ff)
# read ftbl in a list of lines
with open(ftbl, "r") as f:
    lftbl=f.readlines()
# detect flux definition section in ftbl
ifl=[i for (i,l) in enumerate(lftbl) if re.match("^FLUXES\w*(//.*)*", l)][0]
inet=ifl+[i for (i,l) in enumerate(lftbl[ifl:]) if re.match("\tNET\w*(//.*)*", l)][0]
ixch=inet+[i for (i,l) in enumerate(lftbl[inet:]) if re.match("\tXCH\w*(//.*)*", l)][0]
iend=ixch+[i for (i,l) in enumerate(lftbl[ixch:]) if re.match("[^ \t\r\n/]+", l)][0]
#print(ifl,inet,ixch)

# detect metabolite definition section in ftbl
ipool=[i for (i,l) in enumerate(lftbl) if re.match("^METABOLITE_POOLS\w*(//.*)*", l)]
ipool=ipool[0] if ipool else None

# get metab_scale if any
#try:
#    metab_scale=[float(eval(l.split("\t")[2])) for l in lftbl if l.startswith("\tmetab_scale\t")][0]
#except:
#    metab_scale=1.

# main part
# update values in ftbl
for (fl, v) in ff.items():
    try:
        (fdc, nx, flu)=fl.split(".", 2)
    except:
        try:
            (fdc, flu)=fl.split(":", 2)
            nx=None
        except:
            continue
    #if fdc != "f" and fdc != "pf":
    #    continue;
    if nx == "n":
        # replace the value of net free or dependent flux
        iflu=[i for (i,l) in enumerate(lftbl[inet:ixch]) if l.startswith("\t\t%s\tF\t"%flu) or l.startswith("\t\t%s\tD\t"%flu)]
        if len(iflu)==1:
            iflu=inet+iflu[0]
            #print("iflu=", iflu)
        else:
            continue;
        lftbl[iflu]=re.sub("\t\t[^\t]+\t(.)\t.*", "\t\t%s\t\\1\t%.15g"%(flu, v), lftbl[iflu])
        #print("edited %d=%g"%(iflu,v))
    elif nx == "x":
        # replace the value of xch flux
        iflu=[i for (i,l) in enumerate(lftbl[ixch:iend]) if l.startswith("\t\t%s\tF\t"%flu) or l.startswith("\t\t%s\tD\t"%flu)]
        if len(iflu)==1:
            iflu=ixch+iflu[0]
            #print("iflu=", iflu)
        else:
            continue;
        #lftbl[iflu]="\t\t%s\tF\t%.15g%s"%(flu, abs(v), os.linesep)
        lftbl[iflu]=re.sub("\t\t[^\t]+\t(.)\t.*", "\t\t%s\t\\1\t%.15g"%(flu, abs(v)), lftbl[iflu])
    elif fdc == "pf":
        # replace the value of pool
        ipfu=[i for (i,l) in enumerate(lftbl[ipool:]) if l.startswith("\t%s\t-"%flu)]
        if len(ipfu)==1:
            ipfu=ipool+ipfu[0]
            #print("ipfu=", ipfu)
        else:
            continue;
        lftbl[ipfu]="\t%s\t%.15g%s"%(flu, -v, os.linesep)

sys.stdout.write("".join(lftbl))
