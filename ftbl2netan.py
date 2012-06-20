#!/usr/bin/env python
"""Parse ftbl file from stdin or from first parameter
and write netan in kvh format on stdout
usage: ftbl2netan.py [network[.ftbl]] > netan.txt
"""
if __name__ == "__main__":
    import sys, os, getopt;
    import tools_ssg;
    import C13_ftbl;
    def usage():
        print(__doc__);

    try:
        opts,args=getopt.getopt(sys.argv[1:], "h", ["help", "DEBUG"]);
    except getopt.GetoptError, err:
        print str(err);
        usage();
        sys.exit(1);
    DEBUG=False;
    for o,a in opts:
        if o in ("-h", "--help"):
            usage();
            sys.exit(0);
        elif o=="--DEBUG":
            DEBUG=True;
        else:
            assert False, "unhandled option";
    fftbl=args[0] if len(args) else "";
    if fftbl and fftbl[-5:] != ".ftbl":
        fftbl+=".ftbl";
    if fftbl and not os.path.exists(fftbl):
        sys.stderr.write(me+": file '"+fftbl+"' does not exist.\n");
        sys.exit(1);
    fftbl=open(fftbl, "r") if fftbl else sys.stdin;

    ftbl=C13_ftbl.ftbl_parse(fftbl);
    netan=C13_ftbl.ftbl_netan(ftbl);
    tools_ssg.dict2kvh(netan);
    # calculate measure matrices
    if "measures" not in netan:
        measures=dict();
        for meas in ("label", "mass", "peak"):
            measures[meas]=eval("C13_ftbl.%s_meas2matrix_vec_dev(netan)"%meas);
        tools_ssg.dict2kvh({"measures": measures});
    tools_ssg.dict2kvh({"rcumo_sys": C13_ftbl.rcumo_sys(netan)});
    sys.exit(0);
