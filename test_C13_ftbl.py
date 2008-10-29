# 2008-01-23 sokol
# python session for C13_ftbl.py test
#
import pdb;
import sys;
#from itertools import *;

sys.path.append('/home/sokol/dev/python');
sys.path.append('/home/sokol/insa/sysbio/dev/ftbl2sys');
from tools_ssg import *;
#from C13_ftbl import *;
import C13_ftbl;

#fname='ex5.ftbl';
#fname='PPP_exact.ftbl';
#fname='example2.ftbl';
fname='antoniewicz-2006/simple.ftbl';

#pdb.set_trace();

#reload(C13_ftbl);
fin=open(fname, "r");
ftbl=C13_ftbl.ftbl_parse(fin);
fin.close();

#aff("ftbl['PROJECT']", ftbl['PROJECT']);
#aff("ftbl['NETWORK']", ftbl['NETWORK']);
#sec='INEQUALITIES';
#sec='MASS_SPECTROMETRY'
sec='PEAK_MEASUREMENTS'
aff("ftbl['"+sec+"']", ftbl[sec]);

#pdb.set_trace();
#reload(C13_ftbl);
netan=C13_ftbl.ftbl_netan(ftbl);
d=C13_ftbl.aglom_loop1(netan["cumo_sys"]["A"][0]);
aff("loop", d["loop"]);
k=C13_ftbl.lowtri(d["na"]);
fp=open("na.pbm", "w");
C13_ftbl.mat2pbm(d["na"], k, fp);
fp.close();
exit(0);

#e="flux_inequal";
#e="label_meas";
#e="peak_meas";
e="Afl";
#aff("matx_meas_label", C13_ftbl.label_meas2matrix_vec_dev(netan));
#e="mass_meas";
#aff("matx_meas_mass", C13_ftbl.mass_meas2matrix_vec_dev(netan));
#aff("matx_meas_peak", C13_ftbl.peak_meas2matrix_vec_dev(netan));
aff("netan["+e+"]", netan[e]);
aff("netan", netan);

# print stoichometric matrix
metabs=sorted(netan['subs']|netan['prods']);
# print header row
f=open("stocheom_matrix.txt","w");
f.write("flux\t" + "\t".join(metabs)+"\n");
for fl in sorted(netan['sto_r_m'].keys()):
    coefs=netan['sto_r_m'][fl];
    f.write("%s" % fl);
    for m in metabs:
        f.write("\t%s" % coefs.get(m,0));
    f.write("\n");

f.close();

# print fwd-rev flux matrix
fwd_rev=[fl+".fwd" for fl in netan['reac']];
fwd_rev+=[fl+".rev" for fl in set(netan['reac'])-set(netan['notrev'])];
fwd_rev.sort();
aff('fwd_rev', fwd_rev);

f=open('fwd_rev_matrix.txt','w');
f.write('flux\t' + '\t'.join(fwd_rev)+'\tb'+'\n');

# stoechiometric part of matrix and b in A*f=b
for (metab,coefs) in netan['flux_m_r'].iteritems():
    if (metab in netan['input'] or metab in netan['output']):
        continue;
    f.write("%s" % metab);
    for fl in fwd_rev:
        f.write("\t%s" % coefs.get(fl,0));
    f.write("\t0.\n");

# fixed uptakes
#upts=[reac for reac in netan['formula'] if netan['formula'][reac]['left']&netan['input']];
# for each measured uptake give its value in b.
# fixme: uptake should go in cost function to minimize, not in b
for (reac,valdev) in netan['flux_measured'].iteritems():
    f.write("%s" % reac+'.uptake');
    for fl in fwd_rev:
        f.write("\t%s" % ('1' if (reac+'.fwd')==fl else '0'));
    f.write("\t%s\n"%valdev['value']);

# xch01 and net constraints (not implemented yet, waiting first real case
# to ask a meaning to constraint xch01 to some value #0 and #1)
# go back from xch01 to xch then to .fwd, .rev
# xch=beta*xch01/(1-xch01) (beta seems to be 1)
# .fwd=xch-min(-net,0)
# .rev=xch-min(vnet,0)

# free net fluxes part
#.fwd-.rev=net
# -1 in b part is meaning less. It is just to show where free values go.
for reac in netan['flux_free']['net']:
    if (reac in netan['flux_measured']):
        continue;
    f.write("%s" % reac+'.net.free');
    for fl in fwd_rev:
        f.write("\t%s" % ('1' if (reac+'.fwd')==fl else ('-1' if (reac+'.rev')==fl else '0')));
    f.write("\t-1.\n");

# free xch fluxes part
# workaround: xch free flux is just replaced by .fwd
# -1 in b part is meaning less. It is just to show where free values go.
for reac in netan['flux_free']['xch']:
    f.write("%s" % reac+'.fwd.free');
    for fl in fwd_rev:
        f.write("\t%s" % ('1' if (reac+'.fwd')==fl else '0'));
##        print i, reac, fl, (reac+'.fwd')==fl;
##        print 's=', ('1' if (reac+'.fwd')==fl else '0');
    f.write("\t-1.\n");

f.close();

# cumomer balance equations
pdb.set_trace();
for w in xrange(1,netan['Cmax']+1):
    f=open('cumo_w'+str(w)+'_system.txt', 'w');
##    cumos=netan['cumo_sys']['A'][w-1].keys(); # weight 1 equations have all metabolites
    # order cumos along pathways
    # starts are cumomers having input cumomers in rhs
    starts=[cumo for cumo,fluxes in netan['cumo_sys']['b'][w-1].iteritems() \
        if any((cur_cumo.split(':')[0] in netan['input']) \
        for cur_cumo in valval(fluxes.values()))];
    # complete starts by all others cumomers
    starts+=[c for c in netan['cumo_sys']['A'][w-1] if not c in starts];
    cumo_paths=C13_ftbl.cumo_path(starts, netan['cumo_sys']['A'][w-1], set());
    # order
    cumos=[cumo for cumo in valval(cumo_paths)];
    #d=[c for c in netan['cumo_sys']['A'][w-1] if not c in cumos]
    if len(cumos) != len(netan['cumo_sys']['A'][w-1]):
        raise "wrongCumomerNumber";
    #metab_paths=netan['metab_paths']; # ordered by pathways
    # first row is column names
    f.write("cumomers\t%s\tb\n" % '\t'.join(cumos));
    for cumo in cumos:
        # first column is cumomer name
        f.write("%s\t" % cumo);
        # matrix A
        for cur_cumo in cumos:
            term=netan['cumo_sys']['A'][w-1][cumo].get(cur_cumo,'0');
            f.write("%s\t" % term);
        # last column is the rhs b
        term='';
        if cumo in netan['cumo_sys']['b'][w-1]:
            term='+'.join(flux+'*('+'*'.join(cumo_contibs)+')' \
                for flux,cumo_contibs in netan['cumo_sys']['b'][w-1][cumo].iteritems());
        if not term:
           term='0';
        f.write("%s\n" % term);
    
    f.close();

