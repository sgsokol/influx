# 2019-09-04 sokol@insa-toulouse.fr

# adapt tests to new dir structure
# replace symlinks on ftbl's by a real copy
find -type l \( -name '*'.ftbl -o -name '*'.txt -o -name '*'.R \) | while read f; do
   ref=$(readlink "$f")
   fref=$(readlink -e "$f")
   if [ "x$fref" == "x" ]; then
      # the link points to nothing => add "../" in front of it
      if [ -e "../$ref" -a "x$(readlink -e ../$ref)" != "x" ]; then
         echo "cp '../$ref' '$f'" | tee >> ftbl_origin.txt
         rm -f "$f" && cp "../$ref" "$f" 
      else
         echo 1>&2 "could not unref '../$ref' for '$f'"
      fi
   else
      echo "cp '$fref' '$f'" | tee >> ftbl_origin.txt
      rm -f "$f" && cp "$fref" "$f"
   fi
done

# 2019-09-10 sokol@insa-toulouse.fr
# make hard links for e_coli.ftbl (and other repeated files) for concurrent running

# in python3
import csv, os
fd=csv.reader(open("cases_influx_si-v4.4.4.tab"), delimiter="\t")
for l in fd: #print(l[2] if len(l)>2 else "")
   if len(l) < 3:
      continue
   fn=l[2].split()[-1]
   if not fn.endswith(".ftbl"):
      fn += ".ftbl"
   if not os.path.exists(fn):
      continue
   # hard link fn by adding 1,2,.... before ftbl
   i = 1
   while True:
      ft = fn[:-5]+"_%d_.ftbl"%i
      if not os.path.exists(ft):
         break
      i += 1
   # we have fond free name slot
   print("linking '%s' (source) -> '%s' (dest)"%(fn,ft))
   os.link(fn, ft)
