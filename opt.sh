date;
direx=/home/sokol/insa/sysbio/dev/ftbl2sys;
me=$(basename $0)
DEBUG=""
[ "$me" = "optd.sh" ] && DEBUG="DEBUG"

$direx/ftbl2optR.py $1 $DEBUG &&
   R CMD SHLIB $1.f &&
   date && R --no-save --silent --args --meth $2 < $1.R \
   > $1.log 2> $1.err;
date;
