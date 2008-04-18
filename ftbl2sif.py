#!/usr/bin/python
# read a .ftbl file from stdin and translate it to .sif file for Cytoscape
# on stdout.
# usage: ftbl2sif.py < fnetw.ftbl > fnetw.sif
# 
# 2008-01-23 sokol


import sys;
import re;
from C13_ftbl import ftbl_parse;

sys.path.append('/home/sokol/dev/python')
from tools_ssg import *;

ftbl=ftbl_parse(sys.stdin);
#aff("ftbl", ftbl);
for nflux in ftbl['NETWORK']:
	flux=ftbl['NETWORK'][nflux];
	# skip carbon section
	if not len(flux.get('FLUX_NAME')): continue
	path=re.sub(r'[0-9]+$', '', flux['FLUX_NAME']);
	if flux.get('EDUCT_1') and flux.get('PRODUCT_1'):
		print flux['EDUCT_1'], path, flux.get('PRODUCT_1');
	if flux.get('EDUCT_2') and flux.get('PRODUCT_1'):
		print flux['EDUCT_2'], path, flux.get('PRODUCT_1');
	if flux.get('EDUCT_1') and flux.get('PRODUCT_2'):
		print flux['EDUCT_1'], path, flux.get('PRODUCT_2');
	if flux.get('EDUCT_2') and flux.get('PRODUCT_2'):
		print flux['EDUCT_2'], path, flux.get('PRODUCT_2');
