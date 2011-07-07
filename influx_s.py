#!/usr/bin/env python
"""Optimize free fluxes of a given static metabolic network to fit
13C data defined in a FTBL file."
"""
import sys, os, datetime as dt, subprocess as subp
from optparse import OptionParser

def now_s():
    return(dt.datetime.strftime(dt.datetime.now(), "%Y-%m-%d %H:%M:%S"))

pla=sys.platform
# shared object suffix
so="dll" if pla=="win32" or pla=="cygwin" else "dylib" if pla=="darwin" else "so"
# my own name
me=os.path.realpath(sys.argv[0])
# my exec dir
direx=os.path.dirname(me)
direx="." if not direx else direx

# my version
version=file(os.path.join(direx, "influx_version.txt"), "r").read().strip()

# valid options for python
pyopt=set(("--fullsys", "--DEBUG"))

# create a parser for command line options
parser = OptionParser(usage="usage: %prog [options] /path/to/FTBL_file",
    description=__doc__,
    version="%prog "+version)
parser.add_option("--noopt", action="store_true",
    help="no optimization, just use free fluxes as is, to calculate dependent fluxes, cumomers, stats and so on")
parser.add_option("--noscale", action="store_true",
    help="no scaling factors to optimize => all scaling factors are assumed to be 1")
parser.add_option("--meth", type="choice",
    choices=["BFGS", "Nelder-Mead", "nlsic"],
    help="method for optimization, one of nlsic|BFGS|Nelder-Mead. Default: nlsic")
parser.add_option("--fullsys", action="store_true",
    help="calculate all cumomer set (not just the reduced one necesary to simulate measurements)") 
parser.add_option("--irand", action="store_true",
    help="ignore initial approximation for free fluxes from FTBL file and use random values instead")
parser.add_option("--sens",
    help="sensitivity method: SENS can be 'mc[=N]', mc stands for Monte-Carlo. N is the number of Monte-Carlo simulations. Default for N: 10")
parser.add_option("--cupx", type="float",
    help="upper limit for reverse fluxes. Must be in interval [0, 1]")
parser.add_option("--cupn", type="float",
    help="upper limit for net fluxes")
parser.add_option("--clownr", type="float",
    help="lower limit for not reversible free and dependent fluxes. Zero value (default) means no lower limit")
parser.add_option("--np", type="int",
    help="""Number of parallel process used in Monte-Carlo simulations
    Without this option or for NP=0 all available cores in a given node are used""")
parser.add_option("--ln", action="store_true",
    help="Approximate least norm solution is used for increments during the non-linear iterations when Jacobian is rank deficient")
parser.add_option("--zc", action="store_true",
    help="Apply zero crossing strategy for net fluxes")
parser.add_option("--DEBUG", action="store_true",
    help="developer option")
parser.add_option("--TIMEIT", action="store_true",
    help="developer option")
parser.add_option("--prof", action="store_true",
    help="developer option")

# parse commande line
(opts, args) = parser.parse_args()
#print ("opts=", opts)
#print ("args=", args)
if len(args) != 1:
    parser.print_help()
    parser.error("FTBL_file expected in argument")
ft=args[-1]
if ft[-5:] != ".ftbl":
    ft=ft+".ftbl"
if not os.path.exists(ft):
    parser.error("FTBL file '%s' does not exist."%ft)
f=ft[:-5]
flog=open(f+".log", "w")
ferr=open(f+".err", "w")
flog.write(" ".join('"'+v+'"' for v in sys.argv)+"\n")

# now parse commandArgs from the FTBL file
cmd=""
for line in open(ft, "rb"):
    if line[:13] == "\tcommandArgs\t":
        cmd=line[13:]
        break
(cmd_opts, args) = parser.parse_args(cmd.split())
#print ("cmd_opts=", cmd_opts)
if len(args) != 0:
    parser.print_help()
    parser.error("FTBL_file cannot be given in commandArgs option")

# update cmd_opts with runtime options
cmd_opts._update_loose(dict((k,v) for (k,v) in eval(str(opts)).iteritems() if not v is None))
#print("cmd_opts=", cmd_opts)

lopts=[("--"+k, v if type(v)!=type(True) else None) for (k,v) in eval(str(cmd_opts)).iteritems() if not v is None]
#print("lopts=", lopts)
lopts=[str(v) for t in lopts for v in t if not v is None]
#print("lopts2=", lopts)

# code generation
s="code gen: "+now_s()
flog.write(s+"\n")
flog.flush()
print(s)

try:
    # compile static fortran functions if the shared lib is inexistent or too old
    fcumo=os.path.join(direx, "cumo.")
    if not os.path.exists(fcumo+so) or os.path.getmtime(fcumo+"f") >= os.path.getmtime(fcumo+so):
        p=subp.check_call("R CMD SHLIB --clean ".split()+["cumo.f"], cwd=direx, stdout=flog, stderr=ferr)
        flog.flush()

    # generate the R code
    # extract just python options
    opt4py=list(pyopt.intersection(lopts))+[ft]
    pycmd=["python", os.path.join(direx, "ftbl2optR.py")] + opt4py
    #print(pycmd)
    p=subp.check_call(pycmd, stdout=flog, stderr=ferr)
    flog.flush()

    # execute R code
    s="calcul  : "+now_s()
    flog.write(s+"\n")
    flog.flush()
    print(s)
    rcmd="R --no-save --no-restore --slave --args".split()+lopts
    p=subp.check_call(rcmd, stdin=open(f+".R"), stdout=flog, stderr=ferr)
    flog.flush()
except:
    pass

#end
s="end     : "+now_s()
flog.write(s+"\n")
print(s)

ferr.close()
if os.path.getsize(ferr.name) > 0:
    s="=>Check "+ferr.name
    print(s)
    flog.write(s+"\n")
flog.close()
sys.exit(0)
