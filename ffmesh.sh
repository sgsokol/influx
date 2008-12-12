#!/bin/sh
date;
direx=/home/sokol/insa/sysbio/dev/ftbl2sys;
me=$(basename $0)
DEBUG=""
[ "$me" = "ffmeshd.sh" ] && DEBUG="DEBUG"
fkvh="$1"
base="${fkvh%.kvh}"
$direx/ffmesh.py "$1" $DEBUG &&
   R CMD SHLIB "$base".f &&
   date && R --no-save --silent < "$base".R \
   > "$base".log 2> "$base".err;
date;
