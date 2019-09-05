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

import sys, os, subprocess as subp, re, platform
from time import time, asctime
#from optparse import OptionParser
import argparse as ap
import threading as th
from multiprocessing import cpu_count

def plural_s(n):
    return "s" if n > 1 else ""
def setvar(k, v):
    globals()[k]=v
    return(None)

me=os.path.basename(sys.argv[0])
ondos=platform.system() == "Windows"
# create a parser for command line options
parser = ap.ArgumentParser(usage="usage: %(prog)s [options] tabulated_file.txt",
    description=__doc__)
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
        match=re.match("(\d*):(\d*)", item)
        if match:
            beg=match.group(1) or 1
            if beg < 0:
                beg=nb_te-beg+1
            end=match.group(2) or nb_te
            if end < 0:
                end=nb_te-end+1
            itest+=list(range(int(beg), int(end)+1))
        else:
            itest.append(int(item))
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
itest=sorted(set(itest))
if not itest:
    itest=list(range(1, nb_te+1))

# prepare threading framework
ncavail = cpu_count()
if ncore == 0. or ncore > ncavail:
    ncore = ncavail
elif ncore > 0. and ncore < 1.:
    ncore = ncavail*ncore
ncore = round(ncore)
# run tests
fd_log=open("case_tests.log", "w")
fd_err=open("case_tests.err", "w")
t00=time()
icase=0
nb_ok=0
li_ko=[] # list of failed test: (icase, name)
for line in tests:
    icase+=1
    if icase not in itest:
        continue
    fields=line.split("\t")
    if len(fields) not in (3, 4):
        raise Exception("%s: Three or four fields separated by tabs are expected.\nInstead %d fields are found (%s: %d).\nThe row was '%s'"%(me, len(fields), fcases, icase, line))
    if len(fields) == 3:
        fields.append("")
    (nm_t, retcode, cmd, testcmd)=(item.strip() for item in fields)
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
    fd_log.write("---%s, %d:%s\n"%(asctime(), icase, nm_t))
    fd_err.write("---%s, %d:%s\n"%(asctime(), icase, nm_t))
    fd_log.flush()
    fd_err.flush()
    
    t0=time()
    #devnull=open(os.devnull, "w")
    #import pdb; pdb.set_trace()
    print('%d:%s: running "%s" ...'%(icase, nm_t, cmd))
    if not dry:
        if ondos:
            p=subp.call(cmd, stdout=fd_log, stderr=fd_err, shell=True)
        else:
            p=subp.call(cmd.split(), stdout=fd_log, stderr=fd_err)
    #print p
    ok=dry or (retcode.title()=="True" and p==0) or (retcode.title()=="False" and p!=0) or (retcode.title() not in ("True", "False") and int(retcode)==p)
    addmes=""
    if ok and testcmd and not dry:
        res=[(eval(item), item) for item in testcmd.split(";") if item.strip()]
        #print res
        nb_fail=sum(not t for (t, item) in res)
        if nb_fail:
            addmes=" Failed condition%s: %s."%(plural_s(nb_fail), "; ".join(item for (t, item) in res if not t))
            ok=False
    t1=time()
    if ok:
        print("OK for test %d '%s' (%.2f s)"%(icase, nm_t, t1-t0))
        nb_ok+=1
    else:
        print("=> Test %d '%s' failed! (%.2f s)%s"%(icase, nm_t, t1-t0, addmes))
        li_ko.append((icase, nm_t))

fd_log.close()
fd_err.close()

nb_ko=len(li_ko)
print("---\nIn total, %d/%d test%s failed%s. (%.2f s)"%(nb_ko, nb_ok+nb_ko, plural_s(nb_ko), " ("+", ".join(nm_t for i,nm_t in li_ko)+" ["+",".join(str(i) for i,nm_t in li_ko)+"])" if nb_ko else "", time()-t00))

sys.exit(nb_ko)
