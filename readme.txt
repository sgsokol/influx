# 2008-08-22 sokol
# gradient based optimization

# 2008-08-21 sokol
# debugging residuals
# ir2isc had an error.
./ftbl2optR.py PPP_exact DEBUG
# set initial scale value to sum(v*m)/sum(m*m)

# c13flux on exact ex5
date; ~/insa/sysbio/soft/13CFlux/FBINS.20010920/CooolEvoAlpha -nc 5 -vr 0.25 ex5.ftbl\
   > ex5_c13.txt; date
# ok

# fragment mask in mass section counted carbons in wrong sens
date; ./ftbl2optR.py PPP_exact && R --no-save --slave < PPP_exact_opt.R > PPP_exact_opt.log; date
starting cost value:
[1] 0.2328660
# in c13_flux: 0.232869
# OK => good initial cost value

# try not exact initial approximation
date; ./ftbl2optR.py PPP_s && R --no-save --slave < PPP_s_opt.R > PPP_s_opt.log; date
# convergence is slow

# 2008-08-20 sokol
# recup ppp_exact.ftbl
# first converge from example ppp.ftbl
date; ~/insa/sysbio/soft/13CFlux/FBINS.20010920/CooolEvoAlpha -nc 5 -vr 0.25  PPP_s.ftbl\
   > PPP_s_c13.txt; date
# search for last "inserted in 0"
# go up to
Residuum:	0.232634	with following values
upt__NET:	1.02
emp1__NET:	0.509827
ppp2__XCH:	0.807418
ppp3__XCH:	0.772251
ppp4__XCH:	0.798751
ppp5__XCH:	0.206189
ppp6__XCH:	0.202575
date; ~/insa/sysbio/soft/13CFlux/FBINS.20010920/CooolEvoAlpha -nc 5 -vr 0.25  PPP_s.ftbl\
   > PPP_s_c13.txt; date
Residuum:	0.232871	with following values
upt__NET:	1.02
emp1__NET:	0.509824
ppp2__XCH:	0.80725
ppp3__XCH:	0.718871
ppp4__XCH:	0.798779
ppp5__XCH:	0.212754
ppp6__XCH:	0.202523

# insert this flux values in .ftbl
# see if fluxes move
date; ~/insa/sysbio/soft/13CFlux/FBINS.20010920/CooolEvoAlpha -nc 5 -vr 0.25  PPP_exact.ftbl \
   > PPP_exact.txt; date
mer aoû 20 14:30:12 CEST 2008
mer aoû 20 14:31:14 CEST 2008
# no inserted in 0

# restart C13 flux => another solution (stochastic algorithm)
date; ~/insa/sysbio/soft/13CFlux/FBINS.20010920/CooolEvoAlpha -nc 5 -vr 0.25  PPP_exact.ftbl \
   > PPP_exact.txt; date
mer aoû 20 14:44:35 CEST 2008
mer aoû 20 14:49:14 CEST 2008
# no inserted in 0 => this is solution (no better approximation)

# now our script
date; ./ftbl2optR.py PPP_exact && R --no-save --slave < PPP_exact_opt.R > PPP_exact_opt.log; date


#
# 2008-08-18 sokol
# added initial approx for flux & cumos output
date; ./ftbl2optR.py ex5 && R --no-save --slave < ex5_opt.R > ex5_opt.log; date

# 2008-07-31 sokol
cd ftbl
date; ~/insa/sysbio/soft/13CFlux/FBINS.20010920/CooolEvoAlpha -nc 5 -vr 0.25 PPP_s.ftbl\
   > PPP_s_res.txt; date
date; ./ftbl2optR.py PPP_s && R --no-save --slave < PPP_s_opt.R > PPP_s_opt.log
cd ..
date; ./ftbl2optR.py ex4 && R --no-save --slave < ex4_opt.R > ex4_opt.log; date
date; ~/insa/sysbio/soft/13CFlux/FBINS.20010920/CooolEvoAlpha -nc 5 -vr 0.25 ex4.ftbl\
   > ex4_res.txt; date; date
date; ./ftbl2optR.py ex5 && R --no-save --slave < ex5_opt.R > ex5_opt.log; date

# 2008-07-30 sokol
# trying PPP_s.ftbl
~/insa/sysbio/soft/13CFlux/FBINS.20010920/Ftbl2Flx PPP_s.ftbl
# no error messages

~/insa/sysbio/soft/13CFlux/FBINS.20010920/CumoNet PPP_s.ftbl

# PPP_s diverges
# rewritten cumo_sys without A.B vitual metabolites => still diverges
# added debug out in R code => output flux is reversed (net < 0, xch=0 => fwd=0,rev>0)
# added standart ineqality net >= 0 when xch constrained to zero

# 2008-07-29 sokol
# some tests with C13Flux in dedicated dir
mkdir ftbl
cd ftbl
cp ../ex4.ftbl .
~/insa/sysbio/soft/13CFlux/FBINS.20010920/CooolEvoAlpha -nc 5 -vr 0.25  ex4.ftbl\
   > ex4_c13.txt
# no solution found even with the exact initial approximation
# exact solution is found when initial approx is slightly different
# no unique solution without inequality v2.net >= 0

# converges to negative flux v2.net=-0.5
# is a true solution
# add a constraint v2.net >= 0
./ftbl2optR.py ex4 && R --no-save --slave < ex4_opt.R

# starting from v2.net=1 => divergence
fwd flux vector:
      v1       v2       v3       v4
     1.0 206577.5      0.0      1.0
rev flux vector:
      v1       v2       v3       v4
     0.0      0.0 206576.5      0.0
# v3.net<0 => add constraint v3.net >= 0
# => converge (after correcting a bug in constraints)

# 2008-07-28 sokol
# added scale constraints to ui, ci
# passed to Nelder-Mead min method
./ftbl2optR.py PPP_s
# test generated R code
R --no-save < PPP_s_opt.R > PPP_s_opt.log

cd /home/sokol/insa/sysbio/soft/13CFlux/Example/essai
../../FBINS.20010920/CumoNet PPP.ftbl
# see matrix 
ncl PPP.ftbl.fdb 

# simple network with two branches 0.25 and 0.75
./ftbl2optR.py ex4
# test generated R code
R --no-save < ex4_opt.R > ex4_opt.log
# convergence point depends on starting point => non unique solution
# scaling free params are may be too much

# test with C13_flux
~/insa/sysbio/soft/13CFlux/FBINS.20010920/CooolEvoAlpha -nc 5 -vr 0.25  ex4.ftbl \
   > ex4_c13.txt

# 2008-07-27 sokol
# added inequalities constraints ui, ci
./ftbl2optR.py ex3
# test generated R code
R --save < ex3_opt.R > ex3_opt.log

# 2008-07-24 sokol
# rewriten Afl, bfl in net-xch terms
# added vfluxes to netan
./ftbl2optR.py ex3

# 2008-07-18 sokol
# test ftbl2optR.py
./ftbl2optR.py ex3

# 2008-07-11 sokol
# start ftbl2optR.py script writing R code for optimization problem

# 2008-06-12 sokol
# add flux inequality analysis in C13_ftbl.py
# add label measures analysis in C13_ftbl.py

# 2008-06-10 sokol
# R profiling
# in R :
Rprof(paste(pref, ".Rprof", sep=""));
Rprof(NULL);

# in shell (profile)
date; R --no-save --silent --args --prof --pref PPP < cumo_solv.R > PPP_solv.log; date;
R CMD Rprof PPP.Rprof > PPP_matrid_prof.txt

# no profile matrid
date; R --no-save --silent --args --pref PPP < cumo_solv.R > PPP_solv.log; date;
mar jun 10 15:33:08 CEST 2008
mar jun 10 15:33:14 CEST 2008
6 s

# dense matrix inversion
mar jun 10 15:18:52 CEST 2008
mar jun 10 15:18:58 CEST 2008
6 s

# 10 times execution
# matrid
mar jun 10 15:45:17 CEST 2008
mar jun 10 15:46:05 CEST 2008
48 s

# dense
mar jun 10 15:47:10 CEST 2008
mar jun 10 15:47:56 CEST 2008
46 s
# no visible gain, wait and see for minimization problem

# profile dense
date; R --no-save --silent --args --prof --pref PPP < cumo_solv.R > PPP_solv.log; date;
R CMD Rprof PPP.Rprof > PPP_dense_prof.txt

# 2008-06-09 sokol
# smw.solve() est reécrit comme une qr.solve() pour la classe
# matridm qui étend la classe matrid par ajout de m colonnes
date; R --no-save --args --pref PPP < cumo_solv.R > PPP_solv.log; date;
lun jun  9 18:28:49 CEST 2008
lun jun  9 18:28:56 CEST 2008
# 7 s with matrid solver, not better than dense matrix

# 2008-06-06 sokol
# fin de debugage R smw.solve() pour l'inversion de la matrice presque
# tridiagonale

# 2008-05-28 sokol
# test de package R tridiag_qr
mkdir ~/R/lib
R CMD build matrid
R CMD INSTALL --with-package-versions -l ~/R/lib matrid_1.0.tar.gz

# 2008-04-24 sokol
# temps de calcul pour example1
date; R --no-save --args --pref example1 < cumo_solv.R > example1_solv.log; date;
jeu avr 24 10:22:27 CEST 2008
jeu avr 24 10:22:34 CEST 2008
# 7 s

# temps de calcul pour PPP
# generate R code
chmod 700 PPP_sym.R;
./ftbl2symAb.py < PPP.ftbl > PPP_sym.R;
chmod 440 PPP_sym.R

# generate flux matrix
chmod 700 PPP_fl_matrix.txt;
./ftbl2fl_matrix.py < PPP.ftbl > PPP_fl_matrix.txt;
chmod 440 PPP_fl_matrix.txt;

# execute R
date; R --no-save --args --pref PPP < cumo_solv.R > PPP_solv.log; date;
jeu avr 24 15:40:14 CEST 2008
jeu avr 24 15:40:43 CEST 2008
# 29 s dense qr solve

ven avr 25 13:57:05 CEST 2008
ven avr 25 13:57:27 CEST 2008
# 22 s without matrix print

ven avr 25 17:21:10 CEST 2008
ven avr 25 17:21:17 CEST 2008
# 7 s without Matrix package neither prints

# ex3 with 2 to 2 reaction
chmod 700 ex3_fl_matrix.txt;
./ftbl2fl_matrix.py < ex3.ftbl > ex3_fl_matrix.txt;
chmod 440 ex3_fl_matrix.txt;
# generate R code
chmod 700 ex3_sym.R;
./ftbl2symAb.py < ex3.ftbl > ex3_sym.R;
chmod 440 ex3_sym.R

R --no-save --args --pref ex3 < cumo_solv.R > ex3_solv.log;


# 2008-04-23 sokol
# in cumo_solv.R passage to sparse QR
R --no-save --args --pref example1 < cumo_solv.R > example1_solv.log

# 2008-04-22 sokol
# first solve of cumomer system

# prepare R code
chmod 700 example1_sym.R;
./ftbl2symAb.py < example1.ftbl > example1_sym.R;
chmod 440 example1_sym.R

# use generated R code
R --no-save --args --pref example1 < cumo_solv.R > example1_solv.log

# get FAKUB soft for solving almost tridiagonal systems
wget http://www.netlib.org/toms/470 -O fakub.txt

sed -e '1,3 {s/C/#/}; 1 {i\
#!/bin/sh
}' fakub.txt > \
   fakub.sh

chmod 755 fakub.sh
./fakub.sh

