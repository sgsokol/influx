#!/usr/bin/env python3
"""Test some script "a la" unit tests. Test cases are read from
the input file given as the first and unique argument.
The input file must be a plain tabulated text with at most 4 and at least 3 columns:
 - case name, e.g. "too many free fluxes"
 - correct return code of the shell command to execute (e.g. True of False). If the the returned value does not match the expected result, the case is reported as failed
 - shell command to be executed and whose result must be tested
 - optional additional python commands separated by ";" to execute for verifying the succes of the test. All commands must evaluate to True on succesful test. Otherwise the case will be reported as failed. Additional commands are not executed if the main command returned unexpected code.

The stdout and stderr of executed test go concatenated to case_tests.log and case_tests.err respectively. Both files are overwritten without prompt, so take care to save them elsewhere before execution if you want to preserve them.
Different tests are separated by a string "---<date>, <number>:<name>\n"

Usage: case_tests.py [options] cases_influx_s.tab
"""

import sys, os, subprocess as subp, re, platform, tempfile
from shutil import copyfileobj
from glob import glob
from time import time, asctime
import argparse as ap
import threading
from multiprocessing import cpu_count
from queue import Queue

import random # for debug only

def eval_item(item, fd):
    """eval python code in item and return its results.
    In case of exception, return None
    """
    try:
        res=eval(item)
    except Exception as e:
        mes="%s\n\tin code: %s\n"%(str(e), item)
        sys.stderr.write(mes)
        fd.write(mes)
        res=None
    return res
def plural_s(n):
    return "s" if n > 1 else ""
def setvar(k, v):
    globals()[k]=v
    return(None)
def worker(ith):
    while True:
        item = q.get()
        #print(c("ith=", ith, "; item=", item))
        if item is None:
            q.task_done()
            break
        li_res[ith].append(do_case(ith, *item))
        q.task_done()
def do_case(ith, icase, line):
    if icase not in itest:
        return None
    fd_out=open(os.path.join(tempdir.name, ("out%%0%dd.txt"%ndig)%icase), "w")
    fd_err=open(os.path.join(tempdir.name, ("err%%0%dd.txt"%ndig)%icase), "w")
    ok=None
    t0=time()
    nm_t="row %s:%d"%(fcases, icase)
    fields=line.split("\t")
    if len(fields) == 3:
        fields.append("")
    try: # capture exception to log them in err.txt
        #import pdb; pdb.set_trace()
        if len(fields) != 4:
            raise Exception("%s: Three or four fields separated by tabs are expected.\nInstead %d fields are found (%s: %d).\nThe row was '%s'"%(me, len(fields), fcases, icase, line))
        (nm_t, retcode, cmd, testcmd)=(item.strip() for item in fields[:4])
        if not nm_t:
            nm_t="row %s:%d"%(fcases, icase)
        if not retcode:
            raise Exception("%s: the return code (second column) for tested command must not be empty. Give at least True (no error) of False (error) (%s: %d)"%(me, fcases, icase))
        if not cmd:
            raise Exception("%s: the command to execute (third column) must not be empty (%s: %d)"%(me, fcases, icase))
        testcmd=testcmd.strip()
        # prepare path in cmd if we are in dos
        if ondos:
            cmd=cmd.replace("/", os.path.sep)
        # make tests
        #print cmd
        fd_out.write("---%s, %d:%s\n"%(asctime(), icase, nm_t))
        fd_err.write("---%s, %d:%s\n"%(asctime(), icase, nm_t))
        fd_out.flush()
        fd_err.flush()
        
        #devnull=open(os.devnull, "w")
        #if cmd == "echo yes > tmp_case2.txt":
        #    import pdb; pdb.set_trace()
        if ncore > 1:
            print('%d:th%d:%s running "%s" ...'%(icase, ith, nm_t, cmd))
        else:
            print('%d:%s: running "%s" ...'%(icase, nm_t, cmd))
        #pdb.set_trace()
        if not dry:
            p = subp.run(cmd, shell=True, executable='/bin/bash', stdout=fd_out, stderr=fd_err).returncode
            if ncore > 1:
                print('done %d:th%d:%s (%.2f s)'%(icase, ith, nm_t, time()-t0))
        else:
            p = None
    except Exception as e:
        mes="%s\n\tin case: %s\n"%(str(e), icase)
        sys.stderr.write(mes)
        fd_err.write(mes)
        ok = False
        p = None
    if ok is None:
        retcode_ok = dry or (retcode.title()=="True" and p==0) or (retcode.title()=="False" and p!=0) or (retcode.title() not in ("True", "False") and int(retcode)==p)
    else:
        retcode_ok = False
    ok = ok is None and (dry or retcode_ok) # init ok, can be modified by condition tests
    addmes = ""
    if not ok and p == None:
        addmes = " Exception was raised (cf. test_case.err)"
    elif ok and testcmd and not dry:
        with lock: # globals() can be modified in checks via varset()
            res = [(eval_item(item, fd_err), item) for item in testcmd.split(";") if item.strip()]
        #print res
        nb_fail = sum(not t for (t, item) in res)
        if nb_fail:
            addmes = " Failed condition%s: %s."%(plural_s(nb_fail), "; ".join(item for (t, item) in res if not t))
            ok = False
    elif not retcode_ok:
        addmes = " Unexpected code returned by program: '%s' (expected '%s')"%(p, retcode.title())
    t1=time()
    if ok:
        res = "OK for test %d '%s' (%.2f s)\n"%(icase, nm_t, t1-t0)
    else:
        res = "=> Test %d '%s' failed! (%.2f s)%s\n"%(icase, nm_t, t1-t0, addmes)
    fd_out.write(res)
    fd_out.close()
    fd_err.close()
    return (icase, nm_t, ok, res)
me=os.path.basename(sys.argv[0])
ondos=platform.system() == "Windows"
# create a parser for command line options
parser = ap.ArgumentParser(usage="usage: %(prog)s [options] tabulated_file.txt",
    description=__doc__)
if parser:
    parser.add_argument(
    "--itest",
           help="Indexes of tests to be executed. Format: '1:10' -- use only first ten tests; '1,3' -- use the first and third tests; '1:10,15,91:100' -- a mix of both formats is allowed. '3:' means from the third to the end. ':5' means from the first to the fifth test. Negative values counts from the end of case list. E.g. '-1' indicates the last test case. Default: '' (empty, i.e. all provided tests are passed)")
    parser.add_argument(
    "--nmtest",
           help="Names of tests to be executed. Format: 'name1,name2' (a coma separated list) -- use tests who's names are name1'  and 'name2' (cf. the first column of tab file); 'name1:name2' (a begin:end interval) -- use the tests located in tab file between 'name1' and 'name2';  'name1:name2,name3,name4:name5' -- a mix of both formats is allowed. 'name1:' means from the test 'name1' to the end. ':name2' means from the first to the test 'name2'. In a coma separated list, each entry is tested literally against test names, in case of fail, the entries are tried as regular expressions. E.g. a name 'err.*' will fit all test names started with 'err'. Default: '' (empty, i.e. all provided tests are passed). Options --itest and --nmtest are complementary, i.e. a union of both test collections is passed")
    parser.add_argument(
    "--ncore", type=float, default=0.,
           help="core number to use in parallel testing. Default 0 i.e. all available cores. A fractional number between 0 and 1 indicate a fraction of available cores to use.")
    parser.add_argument(
    "-n", "--dry", action="store_true", default=False,
           help="Dry run: show output as if all tests were OK. None of shell command is excecuted")
    parser.add_argument(
    "fcases", 
           help="The file name of test cases")

# parse commande line
args = parser.parse_args()
for k,v in args.__dict__.items():
    globals()[k]=v
# read cases and store them
tests=[]
for line in open(fcases, "r"):
    # strip whites
    line=line.strip()
    # ignore comments and empty lines
    if not len(line) or line[0] == "#":
        continue
    tests.append(line)
nms=dict((li.split("\t")[0],i+1) for i,li in enumerate(tests))

# now that we know how many tests there are and their names, proceed itest and nmtest
nb_te=len(tests)
if itest:
    ili=itest.split(",")
    itest=[]
    for item in ili:
        item=item.strip()
        match=re.match("(-?\d*):(-?\d*)", item)
        if match:
            try:
                beg=int(match.group(1))
            except ValueError:
                beg=1
            if beg < 0:
                beg=nb_te+beg+1
            try:
                end=int(match.group(2))
            except ValueError:
                end=nb_te
            if end < 0:
                end=nb_te+end+1
            itest+=list(range(int(beg), int(end)+1))
        else:
            i=int(item)
            if i < 0:
                i=nb_te+i+1
            itest.append(i)
if nmtest:
    itest=itest or []
    nmli=nmtest.split(",")
    for item in nmli:
        item=item.strip()
        match=re.match("([^,:]*):([^,:]]*)", item)
        if match:
            beg=match.group(1).strip()
            beg=nms.get(beg, -1) or 1
            if beg < 0:
                raise Exception("Case name '%s' not found in tab file"%match.group(1))
            end=match.group(2).strip()
            end=nms.get(end, -1) or nb_te
            if end < 0:
                raise Exception("Case name '%s' not found in tab file"%match.group(2))
            itest+=list(range(int(beg), int(end)+1))
        else:
            i=nms.get(item, -1)
            # check if item is a regular expression
            ma=[] if i!=-1 else [ic for (nm,ic) in nms.items() if re.match(item, nm)]
            if ma:
                itest+=ma
                continue
            if i < 0:
                raise Exception("Case name '%s' not found in tab file"%item)
            itest.append((i))
if itest:
    itest=sorted(set(itest))
else:
    itest=list(range(1, nb_te+1))
ndig=len(str(itest[-1]))

# prepare tempdir
# stdout of i-th case will go to <tempdir>/case<i>/out.txt, stderr to .../err.txt
tempdir=tempfile.TemporaryDirectory()

# prepare threading framework
ncavail = cpu_count()
if ncore == 0. or ncore > ncavail:
    ncore = ncavail
elif ncore > 0. and ncore < 1.:
    ncore = ncavail*ncore
ncore = min(round(ncore), len(itest))
li_res=[[] for i in range(ncore)] # list of results (icase, nm_t, ok, res)
#import pdb; pdb.set_trace()
lock=threading.Lock()
if ncore > 1:
    #li_f=[Do_case() for i in range(ncore)]
    q=Queue()
    ths=[]
    for i in range(ncore):
        t = threading.Thread(target=worker, args=(i,))
        t.daemon = True # makes the thread interrupt on ctrl-c
        t.start()
        ths.append(t)
    # fill q
    [q.put((i+1,l)) for i,l in enumerate(tests)]

    # run tests
    t00=time()
    q.join() # block until all tasks are done

    # stop workers
    for t in ths:
        q.put(None)
    for t in ths:
        t.join()
else:
    t00=time()
    for (icase, line) in enumerate(tests):
        ith = icase%ncore
        res = do_case(ith, icase+1, line)
        li_res[ith].append(res)
# gather all out's and err's test.out and test.err
for ftype in ["out", "err"]:
    with open("test."+ftype, "wb") as ofd:
        for f in sorted(glob(os.path.join(tempdir.name, "%s*.txt"%ftype))):
            with open(f, "rb") as fd:
                copyfileobj(fd, ofd)
tempdir.cleanup()

nb_ok=sum(res[2] for res_th in li_res for res in res_th if res)
nb_ko=sum(not res[2] for res_th in li_res for res in res_th if res)
li_ko=[res[:2] for res_th in li_res for res in res_th if res and not res[2]]
print("---\nIn total, run %d test%s of which %d failed%s. (%.2f s)"%(nb_ok+nb_ko, plural_s(nb_ko), nb_ko, " ("+", ".join(nm_t for i,nm_t in li_ko)+" ["+",".join(str(i) for i,nm_t in li_ko)+"])" if nb_ko else "", time()-t00))

sys.exit(nb_ko)
