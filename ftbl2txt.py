#!/usr/bin/python
# read a .ftbl file from stdin and translate it to .txt file for Cytoscape
# on stdout.
# txt format for cytoscape is like:
# source target interaction boolean_attribute string_attribute floating_point_attribute
# 
# usage: ftbl2txt.py < fnetw.ftbl > fnetw.txt
# 
# 2008-01-23 sokol


import sys;
import re;
from C13_ftbl import ftbl_parse;

sys.path.append('/home/sokol/dev/python')
from tools_ssg import *;

ftbl=ftbl_parse(sys.stdin);
#aff("ftbl", ftbl);
# print header
outformat="%s\t%s\t%s\t%s\t%s\t%s\t%s";
print outformat % (
	"source",
	"target",
	"interaction",
	"boolean_attribute",
	"source_attribute",
	"target_attribute",
	"floating_point_attribute");

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
