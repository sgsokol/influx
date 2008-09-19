# 2008-09-09 sokol
# full fragments
date; ./isotrace.py GLC#111111 < PPP.ftbl > glc_full.txt; date
mar sep  9 11:24:19 CEST 2008
mar sep  9 12:00:17 CEST 2008
# 355558 paths


# 2008-09-08 sokol
# fragment trace
./isotrace.py A#11 < ex5.ftbl > tmp.txt
./isotrace.py GLC#000001 < PPP.ftbl > glc6.txt
date; ./isotrace.py GLC#000111 < PPP.ftbl > glc6.txt; date
lun sep  8 17:33:04 CEST 2008
lun sep  8 17:33:04 CEST 2008
date; ./isotrace.py GLC#001111 < PPP.ftbl > glc6.txt; date
lun sep  8 17:33:26 CEST 2008
lun sep  8 17:33:35 CEST 2008
# 10k+ paths

date; ./isotrace.py GLC#010111 < PPP.ftbl > glc6.txt; date
lun sep  8 17:36:38 CEST 2008
lun sep  8 17:36:47 CEST 2008
# 10k+ paths

date; ./isotrace.py GLC#100111 < PPP.ftbl > glc6.txt; date
lun sep  8 17:37:52 CEST 2008
lun sep  8 17:38:00 CEST 2008
# 9659 paths

date; ./isotrace.py GLC#100000 < PPP.ftbl > glc1.txt; date
lun sep  8 17:41:43 CEST 2008
lun sep  8 17:41:49 CEST 2008
# 7647 paths

date; ./isotrace.py GLC#110000 < PPP.ftbl > glc1.txt; date

# 2008-09-05 sokol
# python profiling
date; python -m cProfile -o prof.txt ./isotrace.py GLC#100000 < PPP.ftbl > glc1.txt; date; ./profiler_stats.py prof.txt > pstats.txt
# 10k paths in 12 s

# empX reversible
date; ./isotrace.py GLC#100000 < PPP_s.ftbl > glc1.txt; date;

# 2008-09-04 sokol
# comparer Cumonet and isotrace.py
~/insa/sysbio/soft/13CFlux_2005/FBINS.20050329/CumoNet PPP_glc1.ftbl\
   > PPP_glc1.txt;
# plenty of Ery4P

# 2008-09-03 sokol
# modeling label propagation
./isotrace.py A#01 < ex5.ftbl > tmp.txt
./isotrace.py GLC#000001 < PPP.ftbl > glc6.txt
date; ./isotrace.py GLC#100000 < PPP.ftbl > glc1.txt; date

# 2008-09-01 sokol
# perf test on am1 (rémi) metabolic network
mkdir am1
cd am1
cp ~/insa/sysbio/data/remi-2008-01/Network_AM1_SerGlyox_v3.1.ftbl am1.ftbl

date; ~/insa/sysbio/soft/13CFlux_2005/FBINS.20050329/CooolEvoAlpha -nc 5 -vr 0.25 am1.ftbl\
   > am1.txt; date
preparing matrizes A,b : s:0 mus:2
dimStage of stage 0: 17
dimStage of stage 1: 54
dimStage of stage 2: 77
dimStage of stage 3: 63
dimStage of stage 4: 30
dimStage of stage 5: 8
dimStage of stage 6: 1
solving all stages (using prepared ones): s:0 mus:2388
Residuum:	0.321923	with following values
meoh_upt__NET:	3.51091
co2_upt__NET:	4.34556
bf_D2pg__NET:	1.05458
bf_smalate__NET:	0.834649
bf_accoa__NET:	0.399167
MCOAL__XCH:	0.99
SDH__XCH:	0.00343696
cumGroup	value	error	deviation	weight_in_SqS	CumConstraint
gly		0.322	-0.00628276	0.034	0.0341463	$00
gly		0.019	0.0117945	0.031	0.144755	$10
gly		0.605	0.000606201	0.014	0.0018749	$01
gly		0.054	-0.0112709	0.03	0.141147	$11
fluxMeas	value	error	deviation	weight_in_SqS
Group Scales:	Value:
gly:Groupe1	0.994847
inserted in 0 with value 0.321923

# ftbl2opt
date; ../ftbl2optR.py am1 &&
   R CMD SHLIB am1_opt.f &&
   R --no-save --silent --args --prof < am1_opt.R > am1_opt.log &&
   R CMD Rprof am1.Rprof > am1_prof.txt; date;
Erreur dans drop(.Call("La_dgesv", a, as.matrix(b), tol, PACKAGE = "base")) :
  sous-programme Lapack dgesv : le système est exactement singulier
Calls: constrOptim ... trisparse_solv -> solve -> solve.default -> drop -> .Call

# bfgs->nelder-mead
achieved minimum:
[1] 0.3808918
> print(param);
          bf_D2pg          bf_accoa        bf_smalate           co2_upt 
        1.2538893         0.1188700         0.6695074         4.0113320 
         meoh_upt             MCOAL               SDH label;gly;Groupe1 
        3.3414046         0.9552031         0.1451210         1.0104724 
# solution is clearly different from 13c_flux
# our minimum is worse

# nelder-mead -> simulated annealing


# 2008-08-29 sokol
# debugged dense2trid()
> print(cost);
[1] 728.1227
# OK

# now, see time of "smw" method
> res=constrOptim(param, cumo_cost, grad=cumo_grad,
+    ui, ci, mu = 1e-04, control=list(trace=1, maxit=100),
+    method="BFGS", outer.iterations=10, outer.eps=1e-05,
+    no_f, no_w, no_cumos, invAfl, p2bfl, bp, fc,
+    imeas, measmat, measvec, measinvvar, ir2isc, fmn, invfmnvar, ifmn);
initial  value 702.981670 
iter  10 value 114.188376
iter  20 value 35.981999
iter  30 value 35.852185
iter  40 value 35.375709

# while "dense" convergence is different
> res=constrOptim(param, cumo_cost, grad=cumo_grad,
+    ui, ci, mu = 1e-04, control=list(trace=1, maxit=100),
+    method="BFGS", outer.iterations=10, outer.eps=1e-05,
+    no_f, no_w, no_cumos, invAfl, p2bfl, bp, fc,
+    imeas, measmat, measvec, measinvvar, ir2isc, fmn, invfmnvar, ifmn);
initial  value 702.981670 
iter  10 value 114.204253
iter  20 value 35.979473
iter  30 value 35.549511
iter  40 value 34.220132
iter  50 value -1.048248
# pb seems to be error accumulation
# => trid qr is not sufficiently stable
# bad.

# 2008-08-28 sokol
# lin method: dense -> smw
> print(cost);
[1] 439.2654
# bad. It is different from 728.1227 (dense)

# 2008-08-27 sokol
# initialize for matrid passes by fortran

# 2008-08-26 sokol
# profiling
# lin method: dense
date; ./ftbl2optR.py PPP_s && R --no-save --silent --args --prof < PPP_s_opt.R > PPP_s_opt.log; date;
mar aoû 26 10:57:55 CEST 2008
mar aoû 26 11:00:03 CEST 2008
R CMD Rprof PPP_s.Rprof > PPP_s_prof.txt
   %       self        %       total
 self     seconds    total    seconds    name
 38.71     40.58     62.90     65.94     "trisparse_solv"
 15.07     15.80     21.04     22.06     "fwrv2Acumo"
  9.23      9.68     11.88     12.46     "fwrv_x2bcumo"

# ~65% of time is spent in vector and matrix constructions
# bad

# fortran subroutine fwrv2Acumo() to generate the cumomer matrix
date; ./ftbl2optR.py PPP_s &&
   R CMD SHLIB PPP_s_opt.f &&
   R --no-save --silent --args --prof < PPP_s_opt.R > PPP_s_opt.log &&
   R CMD Rprof PPP_s.Rprof > PPP_s_prof.txt; date;
mar aoû 26 14:02:05 CEST 2008
mar aoû 26 14:03:17 CEST 2008
   %       self        %       total
 26.40     15.12     34.36     19.68     "fwrv_x2bcumo"
 18.37     10.52     18.40     10.54     ".Call"
 10.06      5.76     10.37      5.94     "matrix"
  5.97      3.42      5.97      3.42     "*"
# fortran subroutine fwrv2Abcumo() to generate the cumomer matrix AND b
mar aoû 26 18:57:00 CEST 2008
mar aoû 26 18:57:53 CEST 2008
   %       self        %       total
 self     seconds    total    seconds    name
 28.82     11.40     28.87     11.42     ".Call"
 14.16      5.60     14.31      5.66     "matrix"
  9.20      3.64     49.49     19.58     "solve.default"


# 2008-08-25 sokol
# time reduction
# outer-inner iterations = 10, 100

# DEBUG=0
date; ./ftbl2optR.py PPP_s && R --no-save --slave < PPP_s_opt.R > PPP_s_opt.log; date
lun aoû 25 09:56:44 CEST 2008
lun aoû 25 09:58:49 CEST 2008
# time reduction: 1h+ -> 2min

# grad=NULL (R internal differentiation)
Erreur dans dR(theta, theta.old, ...) :
  impossible de trouver la fonction "grad"
Calls: constrOptim -> optim -> <Anonymous> -> gr -> dR
Exécution arrêtée

# method: BFGS -> "L-BFGS-B"
date; ./ftbl2optR.py PPP_s && R --no-save --slave < PPP_s_opt.R > PPP_s_opt.log; date
Erreur dans optim(theta.old, fun, gradient, control = control, method = method,  :
  L-BFGS-B nécessite des valeurs finies de 'fn'
Calls: constrOptim -> optim
Exécution arrêtée
# DEBUG=1 didn't reveal where fn is infinite

# method: BFGS, trace=1 -> 0, DEBUG=0
date; ./ftbl2optR.py PPP_s && R --no-save --slave < PPP_s_opt.R > PPP_s_opt.log; date
lun aoû 25 10:16:43 CEST 2008
lun aoû 25 10:18:46 CEST 2008
# no time reduction due to trace=0

# tridiagsolve: dense -> smw
lun aoû 25 19:04:45 CEST 2008
lun aoû 25 19:15:22 CEST 2008

# 2008-08-22 sokol
# gradient based optimization
# added gradient function to opt_cumo_tols.R
# optim method is BFGS
date; ./ftbl2optR.py PPP_s && R --no-save --slave < PPP_s_opt.R > PPP_s_opt.log; date
ven aoû 22 10:35:51 CEST 2008
ven aoû 22 12:11:04 CEST 2008

# outer iterations 100->10
ven aoû 22 14:10:26 CEST 2008
ven aoû 22 15:04:40 CEST 2008

# inner iterations 100->50
ven aoû 22 16:37:20 CEST 2008
Exécution arrêtée
ven aoû 22 18:18:38 CEST 2008
# still not converged

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

