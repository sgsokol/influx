# 2011-03-06 sokol
# function to manipulate kvh format
obj2kvh=function(obj, objname, conct=stdout(), indent=0) {
   # write in kvh file an R object: scala, vector, matrix and list.
   # other objects are written as an output string of toString() function
   # conct is a connection opened for writing;
   # indent is tab offset for object name;
   #
   # Exemple of use:
   # m=matrix(1:6,2,3);
   # fcn=file("m.kvh", "w");
   # obj2kvh(m, "m", fcn);
   # close(fcn);
   # To add a content to existent file use "a" as open mode
   # fcn=file("m.kvh", "a");
   cls=class(obj);
   indent=max(indent,0);
   open_here=F
   if (class(conct)[1]=="character") {
      # open a file for writing if non existent yet otherwise raise an error
      if (file.exists(conct)) {
         stop(sprintf("Cannot write to existent file '%s'. To overwrite it, open a connection in 'w' or 'wb' mode.", conct))
      }
      conct=file(conct, "wb")
      open_here=T
   }
   if (is.vector(obj) && !is.list(obj) && (length(obj)==1 && is.null(names(obj)))) {
      # unnamed scalar
      cat(rep("\t", indent, sep=""), sep="", file=conct);
      cat(c(esc_kvh_k(objname), "\t", obj, "\n"), sep="", file=conct);
   } else if (cls=="matrix" || substring(cls, nchar(cls)-5)=="Matrix") {
      # matrix name
      cat(rep("\t", indent, sep=""), sep="", file=conct);
      cat(c(esc_kvh_k(objname), "\n"), sep="", file=conct);
      # place for row names
      cat(rep("\t", indent+1, sep=""), sep="", file=conct);
      cat("row_col\t", sep="", file=conct);
      # column names if any else column numbers
      if (dim(obj)[[2]] > 0) {
         if (length(dimnames(obj)[[2]])==dim(obj)[[2]]) {
            cat(dimnames(obj)[[2]], sep="\t", file=conct);
         } else {
            cat(paste("V", 1:dim(obj)[[2]], sep=""), sep="\t", file=conct);
         }
      }
      cat("\n", file=conct);
      # row name followed by the row values
      if (dim(obj)[[1]] > 0) {
         if (length(dimnames(obj)[[1]])==dim(obj)[[1]]) {
            rownms=dimnames(obj)[[1]];
         } else {
            rownms=1:dim(obj)[[1]];
         }
         for (i in 1:dim(obj)[[1]]) {
            cat(rep("\t", indent+1), sep="", file=conct);
            cat(c(rownms[i], obj[i,]), sep="\t", file=conct);
            cat("\n", file=conct);
         }
      }
   } else if (is.vector(obj) && !is.list(obj) && NROW(obj) >= 1) {
      # vector
      cat(rep("\t", indent, sep=""), sep="", file=conct);
      cat(esc_kvh_k(objname), file=conct);
      cat("\n", file=conct);
      # row name followed by the row values
      if (length(names(obj))==NROW(obj)) {
         rownms=names(obj);
      } else {
         rownms=1:NROW(obj);
      }
      for (i in 1:NROW(obj)) {
         cat(rep("\t", indent+1), sep="", file=conct);
         cat(c(rownms[i], esc_kvh_v(obj[i])), sep="\t", file=conct);
         cat("\n", file=conct);
      }
   } else if (is.list(obj)) {
      # list => recursive call if the list is not empty
      cat(rep("\t", indent, sep=""), sep="", file=conct);
      cat(esc_kvh_k(objname), file=conct);
      if (length(obj) > 0) {
         cat("\n", file=conct);
         # row name followed by the row values
         if (length(names(obj))==length(obj)) {
            rownms=names(obj);
         } else {
            rownms=1:length(obj);
         }
         for (i in 1:length(obj)) {
            obj2kvh(obj[[i]], rownms[i], conct, indent+1);
         }
      } else {
         cat("\t\n", file=conct);
      }
   } else if (length(obj)==1 && (class(obj) == "numeric" ||
         class(obj) == "character")) {
      # scalar value
      cat(rep("\t", indent, sep=""), sep="", file=conct);
      cat(c(esc_kvh_k(objname), esc_kvh_v(format(obj[[1]]))), sep="\t", file=conct);
      cat("\n", file=conct);
   } else {
      # unlnown type, write its string value
      cat(rep("\t", indent, sep=""), sep="", file=conct);
      cat(c(esc_kvh_k(objname), esc_kvh_v(format(obj))), sep="\t", file=conct);
      cat("\n", file=conct);
   }
   if (open_here) {
      close(conct)
   }
   return(NULL)
}
esc_kvh_k=function(s) {
   return(gsub("([\t\\\n])", "\\\\\\1", s));
}
esc_kvh_v=function(s) {
   return(gsub("([\\\n])", "\\\\\\1", s));
}
kvh2list=function(fp, lev=0, indent=0) {
#   Read a kvh file from fp stream descriptor
#   and organize its content in a named nested list
#   list(k1=v1, k2=list(k2.1=v2.1))
#   If at some level the first key is row_col
#   then a matrix with names row and columns is tried
#   to be read.
#   If fp is a string, it is used in open() operator
#   for binary reading
# NB: repeted key names in the same hierarchical scope will be
# silently overwritten. So only the last key=value pair will be returned
# NB: all content is of character type.


   # check the stream
   open_here=FALSE;
   if (is.character(fp)) {
      #cat(paste("open ", fp), "\n");
      oldloc=Sys.getlocale("LC_CTYPE");
      Sys.setlocale("LC_CTYPE", "en_US")
      fp=file(fp, "rb");
      seek(fp, 0);
      open_here=TRUE;
   }
   # error control
   if (lev < 0 || indent < 0) {
      stop(sprintf("lev=%d, indent=%d both must be positive", lev, indent));
   }
   if (lev < indent) {
      stop(sprintf("lev=%d, indent=%d, lev must be greater or equal to indent", lev, indent));
   }
#    don't know an equivalent of tell() in R :(
#    if lev > fp.tell():
#        raise NameError("lev=%d, file position=%d, lev must be less or equal to file position"%(lev[0], fp.tell()));
#    if indent[0] > fp.tell():
#        raise NameError("indent=%d, file position=%d, indent must be less or equal to file position"%(indent[0], fp.tell()));
   # algorithm:
   # advance to requested indent (=level)
   # if not sucsessful return an empty list
   # read a key 
   # if sep==\t read value
   # elif sep=\n
   #     recursive call
   #     if no result at the level+1 put empty value
   # else put empty value
   tlist=list();
   key="";
   val="";
   while (TRUE) {
      # make current position to point to the begining of a key
      # so go through an appropriate tab number for the current level
      while (indent < lev) {
         # set file pointer to the requested indentation
         char=readChar(fp, 1);
         if (length(char) && char!="\t") {
            if (length(char) && nchar(char)==1) {
                seek(fp,-1,origin="current");
            }
            break;
         } else if (length(char) == 0) {
            break;
         }
         indent=indent+1;
      }
      if (indent < lev) {
         # we could not read till the requested level
         # so the current level is finished;
         if (open_here) {
            Sys.setlocale("LC_CTYPE", oldloc);
            close(fp);
         }
         return(list(kvh=tlist, indent=indent));
      }
      ksv=kvh_read_keyval(fp);
      key=ksv$key;
      sep=ksv$sep;
      val=ksv$val
      if (sep=="\t") {
         #tlist[[key]]=kvh_read_val(fp);
         tlist[[key]]=val
         indent=0;
      } else if (sep=="\n") {
         lev=lev+1;
         indent=0;
         nextlev=kvh2list(fp, lev, indent);
         nextlist=nextlev$kvh;
         indent=nextlev$indent;
         if (length(nextlist) && names(nextlist)[1] == "row_col") {
            #cat("row_col list names=", names(nextlist), sep="\n");
            cnames=unlist(strsplit(nextlist[["row_col"]], "\t"));
            nc=length(cnames);
            all_col_eq=TRUE;
            m=matrix("", 0, nc);
            for (s in nextlist[-1]) {
               r=unlist(strsplit(s, "\t"));
               if (length(r) != nc) {
                  all_col_eq=FALSE;
                  break;
               }
               m=rbind(m, r);
            }
            if (all_col_eq) {
               dimnames(m)=list(names(nextlist)[-1], cnames);
               nextlist=m;
            }
         }
         lev=lev-1;
         if (length(nextlist)==0) {
             # no value and no deeper level
             tlist[[key]]="";
         } else {
             tlist[[key]]=nextlist;
         }
      } else {
         # we are at the end of file
         if (indent || nchar(key)) {
            tlist[[key]]="";
         }
         if (open_here) {
            Sys.setlocale("LC_CTYPE", oldloc);
            close(fp);
         }
         if (lev) {
            indent=0;
            return(list(kvh=tlist, indent=indent));
         } else {
            return(tlist);
         }
      }
   }
}

kvh_read_key=function(fp) {
   #  Read a string from the current position till the first unescaped
   #  \t, \n or the end of the stream fp.
   #  Return a list (key=key, sep=sep). sep="" at the end of the stream
   key="";
   while (TRUE) {
      char=readChar(fp, 1);
      if (!length(char) || !nchar(char)) {
          return(list(key=key, sep=""));
      } else if (char=="\\") {
         # try to read next char if any
         nextchar=readChar(fp, 1);
         if (!length(nextchar) || !nchar(nextchar)) {
            # end of file
            return(list(key=key,sep=""));
         } else {
            # just add the escaped char
            key=paste(key, nextchar, sep="");
         }
      } else if (char=="\t" || char=="\n") {
         return(list(key=key, sep=char));
      } else {
         # just add a plain char
         key=paste(key, char, sep="");
      }
   }
}
kvh_read_val=function(fp) {
   #  Read a string from current position till the first unescaped
   #  \n or the end of file.
   #  Return the read string.
   val="";
   while (TRUE) {
      char=readChar(fp, 1);
      if (!length(char) || !nchar(char) || char=="\n") {
         return(val);
      } else if (char=="\\") {
          # try to read next char if any
          nextchar=readChar(fp, 1);
          if (!length(nextchar) || !nchar(nextchar)) {
             # end of file
             return(val);
          } else {
             # just add escaped char
             val=paste(val, nextchar, sep="");
          }
      } else {
         # just add a plain char
         val=paste(val, char, sep="");
      }
   }
}
kvh_read_keyval=function(fp) {
   #  Read a string from the current position till the first unescaped
   #  \n or the end of the stream fp.
   #  Return a list (key=key, sep=sep, val=val). sep="" at the end of the stream
   totrow=""
   while (TRUE) { # read till the end of line is unescaped \n or fp end
      row=readLines(fp, 1, warn=F)
      if (length(row)==0 && nchar(totrow)==0) {
         # we were at the end
         return(list(key="", val="", sep=""))
      }
      totrow=paste(totrow, row, sep="", collapse="")
      # see if we need to read further
      rend=regexpr("\\\\*$", row)
#print(paste("row=", row, "rend=", rend))
      if (attr(rend, "match.length") %% 2 == 0) {
         # the row end is not escaped => it is completely read
         break
      } else {
         # append the next row
         next
      }
   }
   # split the totrow in key, sep, val
   rsep=gregexpr("[\\]*\t", totrow)
   isep=attr(rsep[[1]], "match.length")
   isep=which(isep != -1 && isep %% 2 == 1)
   if (length(isep)) {
      # there is a separator tab
      endkey=rsep[[1]][isep[1]]+attr(rsep[[1]], "match.length")[isep[1]]-2
      key=substring(totrow, 1, endkey)
      sep="\t"
      val=substring(totrow, endkey+2)
   } else {
      # no tab separator
      key=totrow
      sep="\n"
      val=NULL
   }
   # strip backslashes
   keyval=as.list(gsub("\\\\(.)", "\\1", c(key=key, val=val)))
   keyval$sep=sep
   return(keyval)
}
kvh_get_matrix=function(f, v) {
   # return a plain matrix stored in the field v of a kvh file f.
   # v can be a vector of subfield keys.
   
   # read the file in a vector
   cont=readLines(f)
   ncont=length(cont)
   # get start line number by grep successively all fields in v
   # and the indent
   indent=0
   nstart=1
   for (it in v) {
      #res=system(sprintf("tail -n +%d '%s' 2>&1| grep -m 1 -n '\\t*%s$'", nstart+1, f, it), intern=T, ignore.stderr=T)
      nstart=nstart+grep(sprintf("\\t*\\Q%s\\E$", it), cont[nstart:ncont])[1]
      indent=which(nchar(strsplit(cont[nstart], "\t", fix=T)[[1]])>0)[1]-1
   }
   # get end number of the matrix in the kvh
   #res=system(sprintf("tail -n +%d '%s' 2>&1| grep -m 1 -n -E ^'%s\\w'", nstart+1, f, join("", rep("\t", indent))), intern=T)
   nend=grep(sprintf("^\\t{0,%d}[^\t]", indent-1), cont[nstart:ncont])
   if (length(nend)) {
      nend=nstart+nend[1]-2
   } else {
      nend=ncont
   }
   d=matrix(unlist(strsplit(cont[nstart:nend], "\t")), byrow=T, nrow=nend-nstart+1)[,-(1:indent), drop=F]
   rownames(d)=d[,1]
   if (rownames(d)[1]=="row_col") {
      colnames(d)=d[1,]
      d=d[-1,-1,drop=F]
   } else {
      d=d[,-1,drop=F]
   }
   dn=dimnames(d)
   wop=options()$warn
   options(warn=-1)
   num<-as.numeric(d)
   options(warn=wop)
   if (!any(is.na(num))) {
      d=matrix(num, nrow=nrow(d))
      dimnames(d)=dn
   }
   return(d)
}
