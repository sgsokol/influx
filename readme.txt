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

