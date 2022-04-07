#!/usr/bin/env python3
"""Parse ftbl file from stdin or from first parameter
and write netan in kvh format on stdout
usage: ftbl2netan.py network[.ftbl] [-h] [-i] [--emu] [--clownr] [--fullsys]  [> network.netan]
"""
if __name__ == "__main__" or __name__ == "influx_si.cli":
    import sys, os, getopt, stat
    sys.tracebacklimit=None
    import influx_si
    import tools_ssg
    import C13_ftbl
    def usage():
        print(__doc__)

    me=os.path.basename(sys.argv[0])
    try:
        opts,args=getopt.getopt(sys.argv[1:], "hi", ["help", "emu", "clownr", "fullsys"])
    except getopt.GetoptError as err:
        print(str(err))
        usage()
        sys.exit(1)
    emu=False
    clownr=False
    fullsys=False
    case_i=False
    for o,a in opts:
        if o in ("-h", "--help"):
            usage()
            sys.exit(0)
        if o=="-i":
            case_i=True
        elif o=="--emu":
            emu=True
        elif o=="--clownr":
            clownr=True
        elif o=="--fullsys":
            fullsys=True
        else:
            assert False, "unhandled option"
    if not args:
        sys.stderr("Expecting ftbl file name\n")
        usage()
        
    C13_ftbl.clownr=clownr
    C13_ftbl.case_i=case_i
    fftbl=args[0]
    if fftbl and fftbl[-5:] != ".ftbl":
        fftbl+=".ftbl"
    if not os.path.exists(fftbl):
        sys.stderr.write(me+": file '"+fftbl+"' does not exist.\n")
        sys.exit(1)

    # what kind of output we have?
    mode=os.fstat(1).st_mode
    f=sys.stdout if stat.S_ISFIFO(mode) or stat.S_ISREG(mode) else  open(fftbl[:-4]+"netan", "w")
    
    ftbl=C13_ftbl.ftbl_parse(fftbl)
    netan=dict()
    #import pdb; pdb.set_trace()
    try:
        C13_ftbl.ftbl_netan(ftbl, netan, emu, fullsys, case_i)
        
    except Exception as e:
        sys.stderr.write("ftbl2netan: Exception\n"+str(e)+"\n")
        tools_ssg.dict2kvh(netan, f)
        raise e
        #sys.exit(1)
    #import pdb; pdb.set_trace()
    tools_ssg.dict2kvh(netan, f)
    # calculate measure matrices
    if "measures" not in netan:
        measures=dict()
        for meas in ("label", "mass", "peak"):
            measures[meas]=eval("C13_ftbl.%s_meas2matrix_vec_dev(netan)"%meas)
        tools_ssg.dict2kvh({"measures": measures}, f)
    if emu:
        tools_ssg.dict2kvh({"rcumo_sys": C13_ftbl.rcumo_sys(netan, emu)}, f)
    else:
        tools_ssg.dict2kvh({"rcumo_sys": C13_ftbl.rcumo_sys(netan)}, f)
    tools_ssg.dict2kvh({"vrcumo": netan.get("vrcumo", [])}, f)
    tools_ssg.dict2kvh({"rcumo2i0": netan.get("rcumo2i0", [])}, f)
    tools_ssg.dict2kvh({"rcumo_input": netan["rcumo_input"]}, f)
    if emu:
        tools_ssg.dict2kvh({"emu_input": netan["emu_input"]}, f)
        tools_ssg.dict2kvh({"vemu_input": netan.get("vemu_input", [])}, f)
        tools_ssg.dict2kvh({"vemu": netan.get("vemu", [])}, f)
        tools_ssg.dict2kvh({"emu2i0": netan.get("emu2i0", {})}, f)
    sys.exit(0)
