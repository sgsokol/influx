#!/usr/bin/env python
"""Optimize free fluxes and optionaly metabolite concentrations of a given static metabolic network defined in an FTBL file to fit 13C data provided in the same FTBL file."
"""
import sys, os, datetime as dt, subprocess as subp
from optparse import OptionParser

def now_s():
    return(dt.datetime.strftime(dt.datetime.now(), "%Y-%m-%d %H:%M:%S"))

pla=sys.platform
# my own name
me=os.path.realpath(sys.argv[0])
# my exec dir
direx=os.path.dirname(me)
direx="." if not direx else direx

# my version
version=file(os.path.join(direx, "influx_version.txt"), "r").read().strip()

# valid options for python
pyopt=set(("--fullsys", "--emu", "--DEBUG"))

# create a parser for command line options
parser = OptionParser(usage="usage: %prog [options] /path/to/FTBL_file",
    description=__doc__,
    version="%prog "+version)
parser.add_option(
"--noopt", action="store_true",
    help="no optimization, just use free parameters as is (after a projection on feasability domain), to calculate dependent fluxes, cumomers, stats and so on")
parser.add_option(
"--noscale", action="store_true",
    help="no scaling factors to optimize => all scaling factors are assumed to be 1")
parser.add_option(
"--meth", type="choice",
    choices=["BFGS", "Nelder-Mead", "ipopt", "nlsic"],
    help="method for optimization, one of nlsic|BFGS|Nelder-Mead. Default: nlsic")
parser.add_option(
"--fullsys", action="store_true",
    help="calculate all cumomer set (not just the reduced one necesary to simulate measurements)")
parser.add_option(
"--emu", action="store_true",
    help="simulate labeling in EMU approach")
parser.add_option(
"--irand", action="store_true",
    help="ignore initial approximation for free parameters (free fluxes and metabolite concentrations) from the FTBL file or from a dedicated file (cf --fseries and --iseries option) and use random values drawn uniformly from [0,1] interval")
parser.add_option(
"--sens",
    help="sensitivity method: SENS can be 'mc[=N]', mc stands for Monte-Carlo. N is an optional number of Monte-Carlo simulations. Default for N: 10")
parser.add_option(
"--cupx", type="float",
    help="upper limit for reverse fluxes. Must be in interval [0, 1]. Default: 0.999")
parser.add_option(
"--cupn", type="float",
    help="upper limit for net fluxes. Default: 1.e3")
parser.add_option(
"--cupp", type="float",
    help="upper limit for metabolite pool. Default: 1.e5"),
parser.add_option(
"--clownr", type="float",
    help="lower limit for not reversible free and dependent fluxes. Zero value (default) means no lower limit")
parser.add_option(
"--cinout", type="float",
    help="lower limit for input/output free and dependent fluxes. Must be non negative. Default: 0")
parser.add_option(
"--clowp",
    help="lower limit for free metabolite pools. Must be positive. Default 1.e-8")
parser.add_option(
"--np", type="int",
    help="""Number of parallel process used in Monte-Carlo simulations. Without this option or for NP=0 all available cores in a given node are used""")
parser.add_option(
"--ln", action="store_true",
    help="Approximate least norm solution is used for increments during the non-linear iterations when Jacobian is rank deficient")
parser.add_option(
"--zc", type="float",
    help="Apply zero crossing strategy with non negative threshold for net fluxes")
parser.add_option(
"--fseries",
       help="File name with free parameter values for multiple starting points. Default: '' (empty, i.e. only one starting point from the FTBL file is used)"),
parser.add_option(
"--iseries",
       help="Indexes of starting points to use. Format: '1:10' -- use only first ten starting points; '1,3' -- use the the first and third starting points; '1:10,15,91:100' -- a mix of both formats is allowed. Default: '' (empty, i.e. all provided starting points are used)")
parser.add_option(
"--seed",
       help="Integer (preferably a prime integer) used for reproducible random number generating. It makes reproducible random starting points (--irand) but also Monte-Carlo simulations for sensitivity analysis (--sens mc=N) if executed in sequential way (--np=1). Default: current system value, i.e. random drawing will be varying at each run.")
parser.add_option(
"--DEBUG", action="store_true",
    help="developer option")
parser.add_option(
"--TIMEIT", action="store_true",
    help="developer option")
parser.add_option(
"--prof", action="store_true",
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
print(" ".join('"'+v+'"' for v in sys.argv))

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
    parser.error("FTBL_file has to be given in the command line")

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
    # generate the R code
    # extract just python options
    opt4py=list(pyopt.intersection(lopts))+[ft]
    pycmd=["python", os.path.join(direx, "ftbl2optR.py")] + opt4py
    #print(pycmd)
    p=subp.check_call(pycmd, stdout=flog, stderr=ferr)
    flog.flush()

    # execute R code
    rcmd="R --vanilla --slave --args".split()+lopts
    flog.write("executing: "+" ".join(rcmd)+" <"+f+".R >"+flog.name+" 2>"+ferr.name+"\n")
    s="calcul  : "+now_s()
    flog.write(s+"\n")
    flog.flush()
    print(s)
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
