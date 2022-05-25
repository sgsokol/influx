#!/usr/bin/env python3

"""
Make correspond measurement row id in labcin file to
those in ftbl file, both files are supposed to be produced by txt2ftbl.py
Labcin file is read from ftbl/OPTIONS/file_labcin.
FTBL is left unchanged while labcin is rewritten in-place without warning.

Copyright 2022, INRAE, INSA, CNRS
Author: Serguei Sokol (sokol at insa-toulouse dot fr)
License: Gnu Public License (GPL) v2 http://www.gnu.org/licenses/gpl.html
"""

import sys, os
from glob import glob # wildcard expansion
from pathlib import Path
import argparse
import influx_si
from C13_ftbl import ftbl_parse, ftbl_netan
from txt2ftbl import tsv2df, try_ext

me=os.path.basename(sys.argv[0] or "ftbl2labcin")
def werr(mes):
    raise Exception(f"{me}: "+str(mes))
def warn(mes):
    sys.stderr.write(f"Warning! {me}: "+str(mes)+"\n")

def ftbl_id(ftbl, d, netan, iprl=0):
    "make row id in labcin equal to those in ftbl"
    #import pdb; pdb.set_trace()
    # read labcin
    flabcin = [it["OPT_VALUE"] for it in d.get("OPTIONS", []) if it["OPT_NAME"] == "file_labcin"]
    if not flabcin:
        return # nothing to do, silent return
    if len(flabcin) > 1:
        werr("the field 'file_labcin' is not unique in OPTIONS of '%s'"%ftbl.name)
    flabcin=ftbl.parent/flabcin[0]
    # read flabcin
    try:
        # default comment='#'
        df_cin = tsv2df(flabcin, append_iline=None)
    except:
        # comment='//'
        df_cin = tsv2df(flabcin, comment="//", append_iline=None)
    # build dict miso rows => ftbl rows
    m2f=dict((co[-1].strip(), str(i+1))
        for i,row in enumerate(ftbl.read_text().split("\n"))
        for li in (row.split("//"),) if len(li) > 1
        for co in (li[-1].split(":"),) if len(co) > 1)
    # get id of measurements in flabcin
    mid = df_cin.iloc[:, 0].to_numpy().astype(str)
    # get id in ftbl
    fid = [row["id"]
        for mtype in ('label_meas', 'peak_meas', 'mass_meas')
        for dit in netan[mtype][iprl].values()
        for d2 in dit.values()
        for row in d2.values()
    ]
    # partial id without last field
    fidp=[":".join(li[:-1]) for v in fid for li in (v.split(":"),)]
    sid=set(fid) # for fast literal presence test
    sidp=set(fidp) # for fast partial presence test
    # produce new id if needed => nid
    nid=fid.copy()
    for iid,rid in enumerate(mid):
        if rid in sid:
            continue; # full match, nothing to do
        # get partial match without last field
        li = rid.split(":")
        if not ":".join(li[:-1]) in sidp:
            # partial match not found, leave it as is
            continue
        # from labcin get miso row number and check for equality
        if li[-1] in m2f:
            li[-1] = m2f[li[-1]]
            nid[iid] = ":".join(li)
    df_cin = df_cin.assign(row_col=nid)
    df_cin.to_csv(flabcin, sep="\t", index=False, quoting=False)

def main(argv=sys.argv):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("ftbl",
        help="FTBL file name to process", nargs="+")
    # parse command line
    opts = parser.parse_args(argv)
    args = opts.ftbl
    # expand wildcard (for Windows OS)
    if len(args) > 0:
        lglob = [glob(f) or [f] for f in args]
        args = [l for f in lglob for l in f]
        # make args unique
        args = sorted(set(args))
    for f in args:
        ftbl=try_ext(f, ["ftbl"])
        try:
            d = ftbl_parse(str(ftbl))
        except Exception as e:
            warn(e)
            continue
        netan = dict()
        ftbl_netan(d, netan, case_i=True)
        ftbl_id(ftbl, d, netan)
        if "OPTIONS" in d and "prl_exp" in netan["opt"]:
            iprl=1 # the main ftbl has 0 index
            for pftbl in netan["opt"]["prl_exp"].split(";"):
                pftbl = pftbl.strip()
                if not pftbl:
                    continue
                fpftbl = try_ext(ftbl.parent/pftbl, ["ftbl"])
                dftbl = ftbl_parse(str(fpftbl))
                ftbl_id(fpftbl, dftbl, netan, iprl)
                iprl += 1
    return 0
if __name__ == "__main__":
    main()
