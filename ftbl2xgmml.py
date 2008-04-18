#!/usr/bin/python
# read a .ftbl file from stdin and translate it to .xgmml file on stdout.
# 
# usage: ftbl2xggml.py < fnetw.ftbl > fnetw.xggml
# 
# 2008-01-24 sokol


import sys;
import re;
from C13_ftbl import ftbl_parse;

sys.path.append('/home/sokol/dev/python')
from tools_ssg import *;

ftbl=ftbl_parse(sys.stdin);
#aff("ftbl", ftbl);
# print header
print '<?xml version="1.0"?>
<!DOCTYPE graph PUBLIC "-//John Punin//DTD graph description//EN" "http://www.cs.rpi.edu/~puninj/XGMML/xgmml.dtd">
<graph
   directed="1"
   id="2"
   Vendor="ftbl2xgmml.py (C) sokol at insa-toulouse dot fr"
   Layout="circular"
   label="Project '+str(ftbl['PROJECT']['NAME'])+', v. '+\
   str(ftbl['PROJECT']['VERSION'])+'">';

for nflux in ftbl['NETWORK']:
	flux=ftbl['NETWORK'][nflux];
	# skip carbon section
	if not len(flux.get('FLUX_NAME')): continue
	path=re.sub(r'[0-9]+$', '', flux['FLUX_NAME']);
	if flux.get('EDUCT_1') : print outformat % (
			flux['EDUCT_1'],
			flux['FLUX_NAME'],
			path,
			"TRUE",
			"metabolite",
			"reaction",
			"");
	if flux.get('EDUCT_2'): print outformat % (
			flux['EDUCT_2'],
			flux['FLUX_NAME'],
			path,
			"TRUE",
			"metabolite",
			"reaction",
			"");
	if flux.get('PRODUCT_1'): print outformat % (
			flux['FLUX_NAME'],
			flux['PRODUCT_1'],
			path,
			"TRUE",
			"reaction",
			"metabolite",
			"");
	if flux.get('PRODUCT_2'): print outformat % (
			flux['FLUX_NAME'],
			flux['PRODUCT_2'],
			path,
			"TRUE",
			"reaction",
			"metabolite",
			"");
# print footer
print "</graph>"
