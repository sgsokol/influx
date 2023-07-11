#!/usr/bin/env python3
"""Optimize free fluxes and optionaly metabolite concentrations of a given static metabolic network defined in an FTBL file to fit 13C data provided in the same FTBL file.
"""

#import pdb

import sys, os, datetime as dt, subprocess as subp, re, time, traceback
import argparse
from threading import Thread # threaded parallel jobs
from multiprocessing import cpu_count
from queue import Queue # threaded parallel jobs
from glob import glob # wildcard expansion
from pathlib import Path

import influx_si
import txt2ftbl
import ftbl2optR

class ArgumentParser(argparse.ArgumentParser):
    def error(self, message):
        #self.print_help(sys.stderr)
        #self.exit(2, '%s: error: %s\n' % (self.prog, message))def now_s():
        raise Exception(me+": ArgumentParser: "+message)
def now_s():
    return(dt.datetime.strftime(dt.datetime.now(), "%Y-%m-%d %H:%M:%S"))
def qworker():
    while True:
        item=q.get()
        #print("item=", item)
        try:
            retcode=launch_job(**item)
        except Exception as e:
            traceback.print_exc()
            sys.stderr.write("Error: launch_job: "+"%s\n"%str(e)+"\n")
            retcode=1
        qret.put(retcode)
        q.task_done()
def move2tmp(dirres, fp):
    "Move file fp to <dirres>/tmp dir"
    if fp:
        fp=Path(fp)
    else:
        return
    if dirres == "default":
        fp.rename(fp.parent/(fp.stem+"_res")/"tmp"/fp.name)
    elif dirres:
        fp.rename(Path(dirres)/"tmp"/fp.name)
def launch_job(ft, fshort, cmd_opts, nb_ftbl, case_i):
    r"""Launch R code generation and then its execution
"""
    #pdb.set_trace()
    ft=Path(ft)
    f=ft.with_suffix("")
    d=Path(ft).parent
    if not d.is_dir():
        sys.stderr.write("Error: directory of FTBL file '%s' does not exist.\n"%d)
        return(1)
    if dirres:
        if dirres == "default":
            dres=d/(str(f.name)+"_res")
        else:
            dres=Path(dirres)
        dres.mkdir(exist_ok=True, parents=True)
        (dres/"tmp").mkdir(exist_ok=True)
        fres=dres/(f.name+".res") # .res will be replaced with .log/.err
        try:
            flog=fres.with_suffix(".log").open("w", encoding="utf8")
        except Exception as e:
            sys.stderr.write("%s\n"%str(e))
            return(1);
        try:
            ferr=fres.with_suffix(".err").open("w", encoding="utf8")
        except Exception as e:
            sys.stderr.write("%s\n"%str(e))
            flog.close()
            return(1);
    else:
        dres=""
        fres=None
        flog=sys.stdout
        ferr=sys.stderr
    
    #pdb.set_trace()
    flog.write(" ".join('"'+v+'"' for v in sys.argv)+"\n")

    # code generation
    s="code gen: "+now_s()+"\n"
    flog.write(s)
    flog.flush()
    sys.stdout.write(fshort+s)
    retcode=0

    try:
        if not ft.is_file():
            sys.stderr.write("Error: FTBL file '%s' does not exist.\n"%ft)
            retcode=1
            if dirres:
                ferr.write("Error: FTBL file '%s' does not exist.\n"%ft)
                flog.close()
                ferr.close()
            return(retcode);
        # parse commandArgs from the FTBL file
        #pdb.set_trace()
        inp=ft.open("rb").read()
        for co in ["utf-8", "latin9", "utf-16", "utf-32"]:
            try:
                inp=inp.decode(co)
                break
            except:
                pass
        inp=inp.encode("utf-8").decode("utf-8")
        cmd=" ".join(re.findall("^\tcommandArgs\t(.*?)(?://.*)?$", inp, re.MULTILINE))
        ldict_opts=dict_opts.copy()
        ftbl_opts=parser.parse_args(cmd.split())
        if len(ftbl_opts.ftbl):
            ferr.write("Warning: argument(s) '%s' from the field commandArgs of '%s' are ignored.\n"%(" ".join(ftbl_opts.ftbl), ft))

        # update ftbl_opts with cmd_opts with runtime options so rt options take precedence
        ftbl_opts=vars(ftbl_opts)
        ftbl_opts.update(dict((k,v) for (k,v) in vars(cmd_opts).items() if not v is None and v is not False and k != "ftbl"))
        ldict_opts=dict((k,v) for k,v in ftbl_opts.items() if v is not None and v is not False and k != "ftbl")
        if "meth" in dict_opts and dict_opts["meth"]:
            ldict_opts["meth"]=",".join(ldict_opts["meth"])

        # generate the R code
        # leave python options as they are and put R options as arguments to --ropts
        #pdb.set_trace()
        opt4py=["--"+kc for kc in pyoptnota.intersection(list(ldict_opts.keys()))] + \
            [item for kc in pyopta.intersection(list(ldict_opts.keys())) if dict_opts[kc] is not None for item in ("--"+kc, str(ldict_opts[kc]))] + \
            ["--ropts", '"' + "\n".join(k+"="+("'"+v+"'" \
            if isinstance(v, type("")) else ("TRUE" if v is True else ("FALSE" \
            if v is False else ("c('"+"', '".join(str(vv) for vv in v)+"')" if type(v) == list else str(v))))) for k,v in ldict_opts.items() if k not in notropt) + '"'] + \
            (["--case_i"] if case_i else []) + ["--dirres", str(dres)] + [str(ft)]
        r_generated=True
        #pdb.set_trace()
        retcode=ftbl2optR.main(["ftbl2optR.main"] + opt4py)
        if retcode != 0:
            r_generated=False
            if dirres:
                s="=>Check "+ferr.name+"\n"
                sys.stdout.write(s)
                flog.write(s)
            # stop here because of error(s)
            return(retcode)
        if r_generated:
            qres.put((str(ft.with_suffix(".R")), ldict_opts.get("nocalc", False))) # one or many R files in // will be launched outside
        else:
            # we are done, end up all writings
            s="end     : "+now_s()+"\n"
            if dirres:
                flog.write(s)
            sys.stdout.write(fshort+s)
            if os.path.getsize(ferr.name) > 0:
                s="=>Check "+ferr.name+"\n"
                sys.stdout.write(s)
                if dirres:
                    flog.write(s)
            return(retcode)
    except Exception as e:
        #pdb.set_trace()
        sys.tracebacklimit=ldict_opts.get("tblimit", 0)
        ferr.write(("".join(traceback.format_exc())) + "\t(in "+ft.stem+")\n")
        return(1)
    return(retcode)

# my own name
me=os.path.realpath(sys.argv[0])
# my install dir
dirinst=os.path.dirname(os.path.realpath(influx_si.__file__))
dirbin=os.path.join(dirinst, "bin")

me=os.path.basename(me)
if me[:8]=="influx_i":
    case_i=True
elif me[:8]=="influx_s":
    case_i=False
else:
    raise Exception("""Cannot determine wether we are in 's' or 'i' case.
My basename must start with 'influx_s' or 'influx_i' instead of '%s'."""%me[:8])

# my version
version=open(os.path.join(dirinst, "influx_version.txt"), "r").read().strip()

# valid options for python
pyopta=set(("tblimit",))
pyoptnota=set(("fullsys", "emu", "clownr", "ffguess"))
# non valid options for R
notropt=set(("tblimit", "mtf", "prefix", "force"))

# create a parser for command line options
ord_args=[]
class ordAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        ord_args.append(("--"+self.dest, values))
class pvalAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        # treat excl_outliers
        if values is None or values == "":
            # given without arg
            setattr(namespace, self.dest, 0.01)
        elif values == -1:
            # default value => not given on command line
            pass
        else:
            # some value is given in arg
            try:
                val=float(values)
            except:
                raise Exception("--excl_outliers expects a float number, instead got '%s'"%str(values))
            if val <= 0. or val >= 1.:
                raise Exception("--excl_outliers expects a float number in ]0, 1[ interval, instead got '%s'"%str(values))
            setattr(namespace, self.dest, val)
parser = ArgumentParser(description=__doc__)
if True: # just to make a collapsible code block
    parser.add_argument(
    "--noopt", action="store_true",
        help="no optimization, just use free parameters as is (after a projection on feasibility domain), to calculate dependent fluxes, cumomers, stats and so on")
    parser.add_argument(
    "--noscale", action="store_true",
        help="no scaling factors to optimize => all scaling factors are assumed to be 1")
    parser.add_argument(
    "--meth", action="append",
        choices=["BFGS", "Nelder-Mead", "nlsic", "pso"],
        help="method for optimization, one of 'nlsic|BFGS|Nelder-Mead|pso'. Default: 'nlsic'. Multiple occurrences of this option can appear on command line. In this case, specified minimization methods are applied successively, e.g. '--meth pso --meth nlsic' means that 'pso' will be used first, then 'nlsic' will take over from the point where 'pso' ends. In case of multiple methods, it is recommended to start with non-gradient methods like 'pso' or 'Nelder-Mead' and make them follow by gradient based methods like 'nlsic' or 'BFGS'. If 'pso' or 'Nelder-Mead' are indeed used as the first method, it is not recommended to combine them with '--zc' option.")
    parser.add_argument(
    "--fullsys", action="store_true",
        help="calculate all cumomer set (not just the reduced one necesary to simulate measurements)")
    parser.add_argument(
    "--emu", action="store_true",
        help="simulate labeling in EMU approach")
    parser.add_argument(
    "--irand", action="store_true",
        help="ignore initial approximation for free parameters (free fluxes and metabolite concentrations) from the FTBL file or from a dedicated file (cf --fseries and --iseries option) and use random values drawn uniformly from [0,1] interval")
    parser.add_argument(
    "--sens",
        help="sensitivity method: SENS can be 'mc[=N]', mc stands for Monte-Carlo. N is an optional number of Monte-Carlo simulations. Default for N: 10")
    parser.add_argument(
    "--cupx", type=float,
        help="upper limit for reverse fluxes. Must be in interval [0, 1]. Default: 0.999")
    parser.add_argument(
    "--cupn", type=float,
        help="absolute limit for net fluxes: -cupn <= netflux <= cupn. Must be non negative. Value 0 means no limit. Default: 1.e3")
    parser.add_argument(
    "--cupp", type=float,
        help="upper limit for metabolite pool. Default: 1.e5"),
    parser.add_argument(
    "--clownr", type=float,
        help="lower limit for not reversible free and dependent fluxes. Zero value (default) means no lower limit")
    parser.add_argument(
    "--cinout", type=float,
        help="lower limit for input/output free and dependent fluxes. Must be non negative. Default: 0")
    parser.add_argument(
    "--clowp",  type=float,
        help="lower limit for free metabolite pools. Must be positive. Default 1.e-8")
    parser.add_argument(
    "--np", type=float,
        help="""When integer >= 1, it is a number of parallel subprocesses used in Monte-Carlo (MC) simulations or for multiple FTBL inputs. When NP is a float number between 0 and 1, it gives a fraction of available cores (rounded to closest integer) to be used. Without this option or for NP=0, all available cores in a given node are used for MC simulations.""")
    parser.add_argument(
    "--ln", action="store_true",
        help="Least norm solution is used for increments during the non-linear iterations when Jacobian is rank deficient")
    parser.add_argument(
    "--sln", action="store_true",
        help="Least norm of the solution of linearized problem (and not just of increments) is used when Jacobian is rank deficient")
    parser.add_argument(
    "--tikhreg", action="store_true",
        help="Approximate least norm solution is used for increments during the non-linear iterations when Jacobian is rank deficient")
    parser.add_argument(
    "--lim", action="store_true",
        help="The same as --ln but with a function limSolve::lsei()")
    parser.add_argument(
    "--zc", type=float,
        help="Apply zero crossing strategy with non negative threshold for net fluxes")
    parser.add_argument(
    "--ffguess", action="store_true",
           help="Don't use free/dependent flux definitions from FTBL file(s). Make an automatic guess."),
    #parser.add_argument(
    #"--fdfit", action="store_true",
    #       help="Choose starting values for free fluxes such that they match at best free and dependent values given in FTBL file."),
    parser.add_argument(
    "--fseries",
           help="File name with free parameter values for multiple starting points. Default: '' (empty, i.e. only one starting point from the FTBL file is used)"),
    parser.add_argument(
    "--iseries",
           help="Indexes of starting points to use. Format: '1:10' -- use only first ten starting points; '1,3' -- use the the first and third starting points; '1:10,15,91:100' -- a mix of both formats is allowed. Default: '' (empty, i.e. all provided starting points are used)")
    parser.add_argument(
    "--seed", type=int,
           help="Integer (preferably a prime integer) used for reproducible random number generating. It makes reproducible random starting points (--irand) but also Monte-Carlo simulations for sensitivity analysis. Default: none, i.e. current system value is used, so random drawing will be varying at each run.")
    parser.add_argument(
    "--excl_outliers", nargs="?", action=pvalAction,
           help="This option takes an optional argument, a p-value between 0 and 1 which is used to filter out measurement outliers. The filtering is based on Z statistics calculated on reduced residual distribution. Default: 0.01.")
    parser.add_argument(
    "--nocalc", action="store_true",
           help="generate an R code but not execute it.")
    parser.add_argument(
    "--addnoise", action="store_true",
           help="Add centered gaussian noise to simulated measurements written to _res.kvh file. SD of this noise is taken from FTBL file"),
if case_i:
    parser.add_argument(
"--time_order",
       choices=[None, "1", "2", "1,2"],
       help="Time order for ODE solving (1 (default), 2 or 1,2). Order 2 is more precise but more time consuming. The value '1,2' makes to start solving the ODE with the first order scheme then continues with the order 2.")
if True:
    parser.add_argument('--version', action='version', version=version)
    # add MTF options, they are passed through to txt2ftbl
    parser.add_argument(
    "--mtf", action=ordAction, help="option passed to txt2ftbl. See help there.")
    parser.add_argument(
    "--prefix", action=ordAction, help="option passed to txt2ftbl. See help there.")
    parser.add_argument(
    "--eprl", action=ordAction, help="option passed to txt2ftbl. See help there.")
    parser.add_argument(
    "--force", action=ordAction, help="option passed to txt2ftbl. See help there.")
    parser.add_argument(
    "-o", "--out",
           help="output directory. Default: basename of input file without suffix + '_res'. If empty, no output directory is created. In this case log and error messages are directed to standard output/error channels. Non empty OUT can be used when only one input file or MTF set is given.")

    # install helper actions
    parser.add_argument(
    "--copy_doc", action="store_true",
        help="copy documentation directory in the current directory and exit. If ./doc exists, its content is silently owerriten.")
    parser.add_argument(
    "--copy_test", action="store_true",
        help="copy test directory in the current directory and exit. If ./test exists, its content is silently owerriten.")
    parser.add_argument(
    "--install_rdep", action="store_true",
        help="install R dependencies and exit.")
    parser.add_argument(
    "--TIMEIT", action="store_true",
        help="developer option: measure cpu time or not")
    parser.add_argument(
    "--prof", action="store_true",
        help="developer option: do time profiling or not")
    parser.add_argument(
    "--tblimit", type=int, default=0,
        help="developer option: set trace back limit for python error messages")
    parser.add_argument(
    "ftbl",
        help="FTBL file name to process", nargs="*")

if len(sys.argv) == 1:
    parser.print_usage(sys.stderr)
    sys.exit(1)
# parse command line
#pdb.set_trace()
opts = parser.parse_args()

dict_opts=vars(opts)
#pdb.set_trace()

do_exit=False
if dict_opts["copy_doc"]:
    do_exit=True
    from distutils.dir_util import copy_tree
    dsrc=os.path.join(dirinst, "doc")
    ddst=os.path.realpath(os.path.join(".", "doc"))
    print("Copy '%(dsrc)s' to '%(ddst)s'"%{"dsrc": dsrc, "ddst": ddst})
    copy_tree(dsrc, ddst, verbose=1)
if dict_opts["copy_test"]:
    do_exit=True
    from distutils.dir_util import copy_tree
    dsrc=os.path.join(dirinst, "test")
    ddst=os.path.realpath(os.path.join(".", "test"))
    print("Copy '%(dsrc)s' to '%(ddst)s'"%{"dsrc": dsrc, "ddst": ddst})
    copy_tree(dsrc, ddst, verbose=1)
if dict_opts["install_rdep"]:
    do_exit=True
    if os.name == 'nt':
        p=subp.Popen(["Rterm", "--ess", "--no-save", "--no-restore"], stdin=subp.PIPE, stdout=sys.stdout, stderr=sys.stderr)
        p.stdin.write(("source('"+"/".join(dirinst.split(os.path.sep)+["R", "upd_deps.R"])+"')\n").encode())
        p.stdin.flush()
        import msvcrt
        while True:
            if not p.poll() is None:
                break
            if msvcrt.kbhit():
                try:
                    buf=input()+"\n"
                except:
                    break
                #print("got on stdin buf='"+buf+"'")
                p.stdin.write(buf.encode())
                p.stdin.flush()
            time.sleep(0.1)
    else:
        p=subp.Popen(["R", "--interactive", "--no-save", "--no-restore"], stdin=subp.PIPE, stdout=sys.stdout, stderr=sys.stderr, bufsize=1)
        p.stdin.write(("source('"+os.path.join(dirinst, "R", "upd_deps.R")+"')\n").encode())
        p.stdin.flush()
        import select
        po=select.poll()
        po.register(sys.stdin, select.POLLIN)
        while True:
            if not p.poll() is None:
                break
            if po.poll(100):
                try:
                    buf=input()+"\n"
                except:
                    break
                #print("got buf='"+buf+"'")
                p.stdin.write(buf.encode())
                p.stdin.flush()
if do_exit:
    sys.exit(0)

args=opts.ftbl
# expand wildcard (for Windows OS)
if len(args) > 0:
    lglob=[glob(f) or [f] for f in args]
    args=[l for f in lglob for l in f]
    # make args unique
    args=sorted(set(args))
# treat MTF options if any
li_ftbl=[]
mtf_opts=[]
if ord_args:
    if case_i:
        mtf_opts += ["--inst"]
    # check that multiple --prefix is not mixed with other --mtf or --eprl
    ocount={}
    [ocount.update([(t[0], ocount.get(t[0], 0)+1)]) for t in ord_args]
    if "--force" in ocount:
        mtf_opts += ["--force"]
        del(ocount["--force"])
    if ocount["--prefix"] > 1:
        if len(ocount) > 1:
            raise Exception("If several occurrences of '--prefix' options are present, no other MTF options can be used (got '%s')"%"', '".join(k for k in ocount if k != "--prefix"))
        for t in ord_args:
            tmp_ftbl=[]
            txt2ftbl.main(mtf_opts+list(t), tmp_ftbl)
            args += tmp_ftbl
            li_ftbl += tmp_ftbl
    else:
        mtf_opts += [v for t in ord_args for v in t]
        #print("mtf_opts=", mtf_opts)
        txt2ftbl.main(mtf_opts, li_ftbl)
        args += li_ftbl
set_ftbl=set(li_ftbl)
    
if len(args) < 1:
    parser.print_help(sys.stderr)
    parser.error("At least one FTBL_file or MTF set expected in argument")
if len(args) > 1 and opts.out is not None and opts.out:
    parser.error("When several FTBL_files or MTF sets are given in argument, --out can not be set to non empty value (got '%s')"%opts.out)
if opts.out is None:
    dirres="default"
elif len(opts.out) == 0:
    dirres=""
else:
    dirres=opts.out
#pdb.set_trace()
#print((" ".join('"'+v+'"' for v in sys.argv)))
#print("cpu=", cpu_count())
np=dict_opts.get("np")
avaco=cpu_count()
if np and np > 0 and np < 1:
    np=int(round(np*avaco))
elif np and np >= 1:
    np=int(round(np))
else:
    np=avaco
#print("np=", np)
tblimit=dict_opts.get("tblimit")
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
    fshort="" if len(args) == 1 else Path(ft).stem+": "

    item={"ft": ft, "fshort": fshort, "cmd_opts": cmd_opts, "nb_ftbl": nb_ftbl, "case_i": case_i}
    q.put(item)
if not ftpr:
    sys.exit(1)
#pdb.set_trace()
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
    if not retcode and not dict_opts["nocalc"]:
        fpar=open("parallel.R", "w")
        fpar.write("""
    suppressPackageStartupMessages(library(parallel))
    suppressPackageStartupMessages(library(Rcpp))
    dirx="%(dirx)s"
    dirres="%(dirres)s" # can be "default" or empty
    source(file.path(dirx, "opt_cumo_tools.R"))
    doit=function(fR) {
       f=substr(fR, 1, nchar(fR)-2)
       fshort=basename(f)
       if (dirres == "default") {
          nm_flog=sprintf("%%s_res/%%s.log", f, fshort)
          nm_ferr=sprintf("%%s_res/%%s.err", f, fshort)
          flog=file(nm_flog, "a")
       } else {
          flog=base::stdout()
          ferr=base::stderr()
       }
       now=Sys.time()
       
       cat("calcul  : ", format(now, "%%Y-%%m-%%d %%H:%%M:%%S"), "\\n", sep="", file=flog)
       if (dirres == "default") {
          close(flog)
       }
       source(fR)
       now=Sys.time()
       if (dirres == "default") {
          flog=file(nm_flog, "a")
       } else {
          flog=base::stdout()
       }
       cat("end     : ", format(now, "%%Y-%%m-%%d %%H:%%M:%%S"), "\\n", sep="", file=flog)
       if (dirres == "default") {
          if (file.info(nm_ferr)$size > 0) {
              cat("=>Check ", nm_ferr, "\\n", sep="", file=flog)
          }
          close(flog)
       }
       return(retcode)
    }
    # build dyn lib
    nodes=%(np)d
    type="PSOCK"
    cl=makeCluster(nodes, type)
    clusterExport(cl, "dirres")
    flist=c(%(flist)s)
    retcode=max(abs(unlist(parLapply(cl, flist, doit))))
    stopCluster(cl)
    q("no", status=retcode)
    """%{
        "np": min(np, len([fr for fr,noc in rfiles if not noc])),
        "flist": ", ".join('"'+f.replace(os.path.sep, "/")+'"' for f,noc in rfiles if not noc),
        "dirx": os.path.join(dirinst, "R").replace(os.path.sep, "/"),
        "dirres": dirres
        })
        fpar.close()
        # execute R code on cluster
        s="//calcul:  "+now_s()+"\n"
        sys.stdout.write(s)
        rcmd="R --vanilla --slave"
        retcode=subp.run(rcmd.split(), stdin=open(fpar.name, "r")).returncode
        # end up writing
        s="//end   : "+now_s()+"\n"
        sys.stdout.write(s)
    #pdb.set_trace()
    for fR,noc in rfiles:
        f=fR[:-2]
        fp=Path(fR)
        if dirres == "default":
            fpe=(fp.parent/(fp.stem+"_res")/(fp.stem+".err"))
            flog=(fp.parent/(fp.stem+"_res")/(fp.stem+".log"))
        else: # non empty dirres with multiple inputs cannot be
            fpe=""
            flog=None
        move2tmp(dirres, fp)
        fp=fp.with_suffix(".ftbl")
        if str(fp) in set_ftbl:
            move2tmp(dirres, fp)

        if fpe and fpe.stat().st_size > 0:
            s="=>Check "+str(fpe)+"\n"
            with flog.open("a") as flp: flp.write(s)
            sys.stdout.write(s)
else:
    fp=Path(ftpr[0]).with_suffix(".R")
    if dirres == "default":
        flog=(fp.parent/(fp.stem+"_res")/(fp.stem+".log"))
        fpe=(fp.parent/(fp.stem+"_res")/(fp.stem+".err"))
    elif not dirres:
        flog=None
        fpe=None
    else:
        flog=(Path(dirres)/(fp.stem+".log"))
        flog=(Path(dirres)/(fp.stem+".err"))
    # execute only one R code
    if not retcode and rfiles:
        f=rfiles[0][0][:-2]
        rcmd="R --vanilla --slave"
        if flog:
            with flog.open("a") as flp: flp.write("executing: "+rcmd+" < "+f+".R\n")
        s="calcul  : "+now_s()+"\n"
        sys.stdout.write(s)
        if flog:
            with flog.open("a") as flp: flp.write(s)
        retcode=subp.run(rcmd.split(), stdin=open(f+".R", "r")).returncode
    s="end     : "+now_s()+"\n"
    sys.stdout.write(s)
    if flog:
        with flog.open("a") as flp: flp.write(s)
    # move .R (and .ftbl if --prefix)
    move2tmp(dirres, fp)
    if ftpr[0] in set_ftbl:
        move2tmp(dirres, ftpr[0])

    #pdb.set_trace()
    if fpe and fpe.stat().st_size > 0:
        s="=>Check "+str(fpe)+"\n"
        sys.stdout.write(s)
        if flog:
            with flog.open("a") as flp: flp.write(s)
sys.exit(retcode)
