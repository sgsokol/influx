#!/usr/bin/env python3
"""Parse ftbl file from first argument
and write the resulting dictionary in kvh format on kvh file
usage: ftbl2kvh.py network[.ftbl] [> network.kvh]
"""
if __name__ == "__main__" or __name__ == "influx_si.cli":
    import sys, os, getopt, stat;
    import influx_si;
    import tools_ssg;
    import C13_ftbl;
    import txt2ftbl;
    
    #import pdb;
    
    def usage():
        print(__doc__);

    me=os.path.basename(sys.argv[0])
    try:
        opts,args=getopt.getopt(sys.argv[1:], "hi", ["help", "prefix="]);
    except getopt.GetoptError as err:
        print(str(err));
        usage();
        sys.exit(1);
    case_i=False
    li_ftbl=[]
    for o,a in opts:
        if o in ("-h", "--help"):
            usage();
            sys.exit(0);
        elif o == "-i":
            case_i=True
        elif o=="--prefix":
            mtf_opts=[]
            if case_i:
                mtf_opts += ["--inst"]
            txt2ftbl.main(mtf_opts+["--prefix", a], li_ftbl)
        else:
            assert False, "unhandled option '"+o+"'";
    fftbl=args[0] if len(args) else li_ftbl[0] if li_ftbl else "";
    if fftbl and fftbl[-5:] != ".ftbl":
        fftbl+=".ftbl";
    if fftbl and not os.path.exists(fftbl):
        sys.stderr.write(me+": file '"+fftbl+"' does not exist.\n");
        sys.exit(1);

    # what kind of output we have?
    #pdb.set_trace()
    mode=os.fstat(1).st_mode
    f=sys.stdout if stat.S_ISFIFO(mode) or stat.S_ISREG(mode) else  open(fftbl[:-4]+"kvh", "w")

    ftbl=C13_ftbl.ftbl_parse(fftbl);
    tools_ssg.dict2kvh(ftbl, f);
    sys.exit(0);
