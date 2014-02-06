#!/usr/bin/env python
"""Parse ftbl file from stdin or from first parameter
and write netan in kvh format on stdout
usage: ftbl2netan.py [network[.ftbl]] > netan.txt
"""
if __name__ == "__main__":
    import sys, os, getopt
    import tools_ssg
    import C13_ftbl
    def usage():
        print(__doc__)

    try:
        opts,args=getopt.getopt(sys.argv[1:], "h", ["help", "emu", "clownr", "fullsys", "DEBUG"])
    except getopt.GetoptError, err:
        print str(err)
        usage()
        sys.exit(1)
    DEBUG=False
    emu=False
    clownr=False
    fullsys=True
    for o,a in opts:
        if o in ("-h", "--help"):
            usage()
            sys.exit(0)
        elif o=="--emu":
            emu=True
        elif o=="--DEBUG":
            DEBUG=True
        elif o=="--clownr":
            clownr=True
        elif o=="--fullsys":
            fullsys=True
        else:
            assert False, "unhandled option"
    C13_ftbl.clownr=clownr
    fftbl=args[0] if len(args) else ""
    if fftbl and fftbl[-5:] != ".ftbl":
        fftbl+=".ftbl"
    if fftbl and not os.path.exists(fftbl):
        sys.stderr.write(me+": file '"+fftbl+"' does not exist.\n")
        sys.exit(1)
    fftbl=open(fftbl, "r") if fftbl else sys.stdin

    ftbl=C13_ftbl.ftbl_parse(fftbl)
    netan=C13_ftbl.ftbl_netan(ftbl, emu, fullsys)
    tools_ssg.dict2kvh(netan)
    # calculate measure matrices
    if "measures" not in netan:
        measures=dict()
        for meas in ("label", "mass", "peak"):
            measures[meas]=eval("C13_ftbl.%s_meas2matrix_vec_dev(netan)"%meas)
        tools_ssg.dict2kvh({"measures": measures})
    if emu:
        tools_ssg.dict2kvh({"rcumo_sys": C13_ftbl.rcumo_sys(netan, emu)})
    else:
        tools_ssg.dict2kvh({"rcumo_sys": C13_ftbl.rcumo_sys(netan)})
    tools_ssg.dict2kvh({"vrcumo": netan["vrcumo"]})
    tools_ssg.dict2kvh({"rcumo2i0": netan["rcumo2i0"]})
    tools_ssg.dict2kvh({"rcumo_input": netan["rcumo_input"]})
    if emu:
        tools_ssg.dict2kvh({"emu_input": netan["emu_input"]})
        tools_ssg.dict2kvh({"vemu_input": netan["vemu_input"]})
        tools_ssg.dict2kvh({"vemu": netan["vemu"]})
        tools_ssg.dict2kvh({"emu2i0": netan["emu2i0"]})
    sys.exit(0)
