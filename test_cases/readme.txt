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
