#!/usr/bin/env python
"""Optimize free fluxes and optionaly metabolite concentrations of a given static metabolic network defined in an FTBL file to fit 13C data provided in the same FTBL file.
"""
import sys, os, datetime as dt, subprocess as subp, re
from optparse import OptionParser
from threading import Thread # threaded parallel jobs
from multiprocessing import cpu_count
from Queue import Queue # threaded parallel jobs
from glob import glob # wildcard expansion

#from pdb import set_trace

def now_s():
    return(dt.datetime.strftime(dt.datetime.now(), "%Y-%m-%d %H:%M:%S"))

def optional_pval(pval):
    def func(option,opt_str,value,parser):
        if parser.rargs and not parser.rargs[0].startswith('-'):
            try:
                val=float(parser.rargs[0])
                parser.rargs.pop(0)
            except:
                val=pval
        else:
            val=pval
        setattr(parser.values,option.dest,val)
    return func

def qworker():
    while True:
        item=q.get()
        #print("item=", item)
        retcode=launch_job(**item)
        qret.put(retcode)
        q.task_done()

def launch_job(ft, fshort, cmd_opts, nb_ftbl, case_i):
    r"""Launch R code generation and then its execution
"""
    #set_trace()
    f=ft[:-5]
    flog=open(f+".log", "wb")
    ferr=open(f+".err", "wb")
    flog.write(" ".join('"'+v+'"' for v in sys.argv)+"\n")

    # code generation
    s="code gen: "+now_s()+"\n"
    flog.write(s)
    flog.flush()
    sys.stdout.write(fshort+s)
    retcode=0

    try:
        if not os.path.exists(ft):
            sys.stderr.write("Error: FTBL file '%s' does not exist.\n"%ft)
            ferr.write("Error: FTBL file '%s' does not exist.\n"%ft)
            retcode=1
            flog.close()
            ferr.close()
            return(retcode);
        # parse commandArgs from the FTBL file
        cmd=" ".join(re.findall("^\tcommandArgs\t(.*?)(?://.*)?$", open(ft, "rb").read(), re.MULTILINE))
        (ftbl_opts, ftbl_args) = parser.parse_args(cmd.split())
        #print ("cmd_opts=", cmd_opts)
        if len(ftbl_args) != 0:
            ferr.write("Warning: argument(s) '%s' from the field commandArgs of '%s' are ignored.\n"%(" ".join(ftbl_args), ft))

        # update ftbl_opts with cmd_opts with runtime options so rt options take precedence
        ftbl_opts._update_loose(dict((k,v) for (k,v) in eval(str(cmd_opts)).iteritems() if not v is None))
        cmd_opts=eval(str(ftbl_opts))
        cmd_opts=dict((k,v) for k,v in cmd_opts.iteritems() if v is not None)
        #print("cmd_opts=", cmd_opts)

        # generate the R code
        # leave python options as are and put R options as argument to --ropts
        opt4py=list(pyopt.intersection("--"+kc for kc in cmd_opts.keys())) + \
            ["--ropts", '"' + "; ".join(k+"="+("'"+v+"'" \
            if isinstance(v, type("")) else "T" if v is True else "F" \
            if v is False else str(v)) for k,v in cmd_opts.iteritems()) + '"'] + \
            (["--case_i"] if case_i else []) + [ft]
        pycmd=["python", os.path.join(direx, "ftbl2optR.py")] + opt4py
        pycmd_s=" ".join(('' if item and item[0]=='"' else '"')+item+('' if item and item[0]=='"' else '"') for item in pycmd)
        flog.write("executing: "+pycmd_s+"\n")
        r_generated=True
        if os.name=="nt":
            retcode=subp.call(pycmd_s, stdout=flog, stderr=ferr, shell=True)
        else:
            retcode=subp.call(pycmd, stdout=flog, stderr=ferr)
        if retcode:
            r_generated=False
        if os.path.getsize(ferr.name) > 0:
            s="=>Check "+ferr.name+"\n"
            sys.stdout.write(s)
            flog.write(s)
            flog.close()
            ferr.close()
            # stop here because of error(s)
            return(retcode)
        flog.close()
        ferr.close()
        if r_generated and "nocalc" not in cmd_opts:
            qres.put(f+".R") # one or many R files in // will be launched outside
        else:
            # we are done, end up all writings
            flog=open(f+".log", "ab")
            s="end     : "+now_s()+"\n"
            flog.write(s)
            sys.stdout.write(fshort+s)
            if os.path.getsize(ferr.name) > 0:
                s="=>Check "+ferr.name+"\n"
                sys.stdout.write(s)
                flog.write(s)
            flog.close()
            return(retcode)
    except:
        pass
    return(retcode)

# my own name
me=os.path.realpath(sys.argv[0])
# my exec dir
direx=os.path.dirname(me)
direx="." if not direx else direx
me=os.path.basename(sys.argv[0])
if me[:8]=="influx_i":
    case_i=True
elif me[:8]=="influx_s":
    case_i=False
else:
    raise Exception("""Cannot determine wether we are in 's' or 'i' case.
My basename must start with 'influx_s' or 'influx_i' instead of '%s'."""%me[:8])

# my version
version=file(os.path.join(direx, "influx_version.txt"), "r").read().strip()

# valid options for python
pyopt=set(("--fullsys", "--emu", "--clownr"))

# create a parser for command line options
parser = OptionParser(usage="usage: %prog [options] /path/to/FTBL_file1 [FTBL_file2 [...]]",
    description=__doc__,
    version="%prog "+version)
parser.add_option(
"--noopt", action="store_true",
    help="no optimization, just use free parameters as is (after a projection on feasibility domain), to calculate dependent fluxes, cumomers, stats and so on")
parser.add_option(
"--noscale", action="store_true",
    help="no scaling factors to optimize => all scaling factors are assumed to be 1")
parser.add_option(
"--meth", type="choice",
    choices=["BFGS", "Nelder-Mead", "nlsic"],
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
"--clowp",  type="float",
    help="lower limit for free metabolite pools. Must be positive. Default 1.e-8")
parser.add_option(
"--np", type="float",
    help="""When integer >= 1, it is a number of parallel threads (on Unix) or subprocesses (on Windows) used in Monte-Carlo (M-C) simulations or for multiple FTBL inputs. When NP is a float number between 0 and 1, it gives a fraction of available cores (rounded to closest integer) to be used. Without this option or for NP=0, all available cores in a given node are used for M-C simulations.""")
parser.add_option(
"--ln", action="store_true",
    help="Least norm solution is used for increments during the non-linear iterations when Jacobian is rank deficient")
parser.add_option(
"--sln", action="store_true",
    help="Least norm of the solution of linearized problem (and not just of increments) is used when Jacobian is rank deficient")
parser.add_option(
"--tikhreg", action="store_true",
    help="Approximate least norm solution is used for increments during the non-linear iterations when Jacobian is rank deficient")
parser.add_option(
"--zc", type="float",
    help="Apply zero crossing strategy with non negative threshold for net fluxes")
parser.add_option(
"--ffguess", action="store_true",
       help="Don't use free/dependent flux definitions from FTBL file(s). Make an automatic guess."),
parser.add_option(
"--fseries",
       help="File name with free parameter values for multiple starting points. Default: '' (empty, i.e. only one starting point from the FTBL file is used)"),
parser.add_option(
"--iseries",
       help="Indexes of starting points to use. Format: '1:10' -- use only first ten starting points; '1,3' -- use the the first and third starting points; '1:10,15,91:100' -- a mix of both formats is allowed. Default: '' (empty, i.e. all provided starting points are used)")
parser.add_option(
"--seed", type="int",
       help="Integer (preferably a prime integer) used for reproducible random number generating. It makes reproducible random starting points (--irand) but also Monte-Carlo simulations for sensitivity analysis. Default: none, i.e. current system value is used, so random drawing will be varying at each run.")
parser.add_option(
"--excl_outliers", action='callback', callback=optional_pval(0.01), dest="excl_outliers",
       help="This option takes an optional argument, a p-value between 0 and 1 which is used to filter out measurement outliers. The filtering is based on Z statistics calculated on reduced residual distribution. Default: 0.01.")
parser.add_option(
"--nocalc", action="store_true",
       help="generate an R code but not execute it.")
if case_i:
    parser.add_option(
"--time_order", type="int",
       help="Integer (1 (default) or 2) indicating the order of time discretization scheme. Order 2 is more precise but more time consuming.")

parser.add_option(
"--TIMEIT", action="store_true",
    help="developer option")
parser.add_option(
"--prof", action="store_true",
    help="developer option")

# parse commande line
(opts, args) = parser.parse_args()
# expand wildcard (for Windows OS)
lglob=[glob(f) or [f] for f in args]
args=[l for f in lglob for l in f]

#print ("opts=", opts)
#print ("args=", args)
# make args unique
args=set(args)
if len(args) < 1:
    parser.print_help()
    parser.error("At least one FTBL_file expected in argument")
dict_opts=eval(str(opts))
if case_i and dict_opts["time_order"] not in (None, 1, 2):
    parser.error("--time_order can take only values 1 or 2. Instead '%d' was given"%dict_opts["time_order"])

print(" ".join('"'+v+'"' for v in sys.argv))
#print("cpu=", cpu_count())
np=dict_opts.get("np")
avaco=cpu_count()
if np > 0 and np < 1:
    np=int(round(np*avaco))
elif np > 1:
    np=int(round(np))
else:
    np=avaco
#print("np=", np)
q=Queue() # arguments for threads
qres=Queue() # results from threads (R file names)
qret=Queue() # returned code from workers

for i in range(min(np, len(args))):
    t=Thread(target=qworker)
    t.daemon=True
    t.start()

nb_ftbl=len(args)
ftpr=[] # proceeded ftbls
(cmd_opts, cmd_args) = (opts, args)
for ft in args:
    if ft[-5:] != ".ftbl":
        ft=ft+".ftbl"
    ftpr.append(ft)
    f=ft[:-5]
    fshort="" if len(args) == 1 else os.path.basename(f)+": "

    item={"ft": ft, "fshort": fshort, "cmd_opts": cmd_opts, "nb_ftbl": nb_ftbl, "case_i": case_i}
    q.put(item)
if not ftpr:
    sys.exit(1)
q.join()

# get names of generated R files
rfiles=[]
while not qres.empty():
    rfiles.append(qres.get())
# get returned codes
rcodes=[]
while not qret.empty():
    rcodes.append(qret.get())

retcode=max(rcodes)
if len(rfiles) > 1:
    # write file parallel.R and launch it
    # //calcul on cluster
    fpar=open("parallel.R", "wb")
    fpar.write("""
    suppressPackageStartupMessages(library(parallel))
    suppressPackageStartupMessages(library(Matrix, warn=F, verbose=F)); # to economize this loading in every parallel worker
    doit=function(fR) {
       f=substr(fR, 1, nchar(fR)-2)
       nm_flog=sprintf("%%s.log", f)
       nm_ferr=sprintf("%%s.err", f)
       fshort=basename(f)
       now=Sys.time()
       flog=file(nm_flog, "ab")
       cat("calcul  : ", format(now, "%%Y-%%m-%%d %%H:%%M:%%S"), "\\n", sep="", file=flog)
       close(flog)
       source(fR)
       now=Sys.time()
       flog=file(nm_flog, "ab")
       cat("end     : ", format(now, "%%Y-%%m-%%d %%H:%%M:%%S"), "\\n", sep="", file=flog)
       if (file.info(nm_ferr)$size > 0) {
           cat("=>Check ", nm_ferr, "\\n", sep="", file=flog)
       }
       close(flog)
       return(retcode)
    }
    nodes=%d
    if (.Platform$OS.type=="unix") {
       type="FORK"
    } else {
       type="PSOCK"
       nodes=rep("localhost", nodes)
    }
    cl=makeCluster(nodes, type)
    flist=c(%s)
    retcode=max(unlist(parLapply(cl, flist, doit)))
    stopCluster(cl)
    q("no", status=retcode)
"""%(min(np, len(rfiles)), ", ".join('"'+f+'"' for f in rfiles)))
    fpar.close()
    # execute R code on cluster
    rcmd="R --vanilla --slave"
    s="//calcul: "+now_s()+"\n"
    sys.stdout.write(s)
    if os.name=="nt":
        retcode=subp.call(rcmd, stdin=open(fpar.name, "rb"), shell=True)
    else:
        retcode=subp.call(rcmd.split(), stdin=open(fpar.name, "rb"))
    # end up writing
    for fR in rfiles:
        f=fR[:-2]
        nm_ferr=f+".err"
        if os.path.getsize(nm_ferr) > 0:
            s="=>Check "+nm_ferr+"\n"
            sys.stdout.write(s)
    s="//end   : "+now_s()+"\n"
    sys.stdout.write(s)
elif len(rfiles)==1:
    # execute only one R code
    f=rfiles[0][:-2]
    rcmd="R --vanilla --slave"
    flog=open(f+".log", "ab")
    flog.write("executing: "+rcmd+" < "+f+".R\n")
    s="calcul  : "+now_s()+"\n"
    flog.write(s)
    flog.close()
    sys.stdout.write(s)
    if os.name=="nt":
        retcode=subp.call(rcmd, stdin=open(f+".R", "rb"), shell=True)
    else:
        retcode=subp.call(rcmd.split(), stdin=open(f+".R", "rb"))
    s="end     : "+now_s()+"\n"
    flog=open(f+".log", "ab")
    flog.write(s)
    sys.stdout.write(s)
    if os.path.getsize(f+".err") > 0:
        s="=>Check "+f+".err"+"\n"
        sys.stdout.write(s)
        flog.write(s)
    flog.close()
sys.exit(retcode)
