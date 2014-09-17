#!/usr/bin/env python
"""Test some script "a la" unit tests. Test cases are read from
the input file given as the first and unique argument.
The input file must be a plain tabulated text with at most 4 and at least 3 columns:
 - case name, e.g. "too many free fluxes"
 - correct return code of the next shell command (e.g. True of False). If the correct result is not returned, the case is reported as failed
 - shell command to be executed and whose result must be tested
 - optional additional python commands separated by ";" to execute for verifying succes of the test. All commands must evaluate to True on succesful test. Otherwise the case will be reported as failed. Additional commands are not executed if the main command returned unexpected code.

The stdout and stderr of executed test go concatenated to case_tests.log and case_tests.err respectively. Both files are overwritten without prompt, so take care to save them elsewhere before execution if you want to preserve them.
Different tests are separated by "---date, number:name\n"
"""

import sys, os, datetime as dt, subprocess as subp, re
from time import time, asctime
from optparse import OptionParser

def plural_s(n):
    return "s" if n > 1 else ""
def setvar(k, v):
    globals()[k]=v
    return(None)

me=os.path.basename(sys.argv[0])
# create a parser for command line options
parser = OptionParser(usage="usage: %prog [options] tabulated_file.txt",
    description=__doc__,
    version="%prog 1.0")
parser.add_option(
"--itest",
       help="Indexes of tests to be executed. Format: '1:10' -- use only first ten tests; '1,3' -- use the the first and third tests; '1:10,15,91:100' -- a mix of both formats is allowed. '3:' means from the third to the end. ':5' means from the first to the fifth test. Default: '' (empty, i.e. all provided tests are passed)")

# parse commande line
(opts, args) = parser.parse_args()
if len(args) != 1:
    parser.print_help()
    parser.error("The file name of test cases must be given as the first and unique argument.")
fcases=args[0]

# read cases and store them
tests=[]
for line in open(fcases, "rb"):
    # strip whites
    line=line.strip()
    # ignore comments and empty lines
    if not len(line) or line[0] == "#":
        continue
    tests.append(line)

# now that we know how many tests there are, proceed itest
nb_te=len(tests)
itest=opts.itest
if itest:
    ili=itest.split(",")
    itest=[]
    for item in ili:
        item=item.strip()
        match=re.match("(\d*):(\d*)", item)
        if match:
            beg=match.group(1) or 1
            end=match.group(2) or nb_te
            itest+=range(int(beg), int(end)+1)
        else:
            itest.append(int(item))
    itest=set(itest)
else:
    itest=set(range(1, nb_te+1))

# run tests
fd_log=open("case_tests.log", "wb")
fd_err=open("case_tests.err", "wb")
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
    
    # make tests
    #print cmd
    fd_log.write("---%s, %d:%s\n"%(asctime(), icase, nm_t))
    fd_err.write("---%s, %d:%s\n"%(asctime(), icase, nm_t))
    fd_log.flush()
    fd_err.flush()
    
    t0=time()
    #devnull=open(os.devnull, "w")
    print 'Running "%s" ...'%cmd
    if os.name=="nt":
        p=subp.call(cmd, stdout=fd_log, stderr=fd_err, shell=True)
    else:
        p=subp.call(cmd.split(), stdout=fd_log, stderr=fd_err)
    #print p
    ok=(retcode.title()=="True" and p==0) or (retcode.title()=="False" and p!=0) or (retcode.title() not in ("True", "False") and int(retcode)==p)
    addmes=""
    if ok and testcmd:
        res=[(eval(item), item) for item in testcmd.split(";") if item.strip()]
        #print res
        nb_fail=sum(not t for (t, item) in res)
        if nb_fail:
            addmes=" Failed condition%s %s: %s."%(plural_s(nb_fail), "are" if nb_fail > 1 else "is", "; ".join(item for (t, item) in res if not t))
            ok=False
    t1=time()
    if ok:
        print "OK for test %d '%s' (%.2f s)"%(icase, nm_t, t1-t0)
        nb_ok+=1
    else:
        print "=> Test %d '%s' has failed! (%.2f s)%s"%(icase, nm_t, t1-t0, addmes)
        li_ko.append((icase, nm_t))

fd_log.close()
fd_err.close()

nb_ko=len(li_ko)
print "---\nIn total, %d/%d test%s failed%s. (%.2f s)"%(nb_ko, nb_ok+nb_ko, plural_s(nb_ko), " ("+", ".join(nm_t for i,nm_t in li_ko)+" ["+",".join(str(i) for i,nm_t in li_ko)+"])" if nb_ko else "", time()-t00)

sys.exit(nb_ko)
