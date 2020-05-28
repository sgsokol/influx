import sys
import codecs
from re import match

#import pdb

def isstr(s):
    return isinstance(s, str)
class Obj():
    def __init__(**kwargs):
        self.__dict__.update(kwrags)
def kvh2tlist(fp, lev=[0], indent=[0], strip=False):
    """
    Read a kvh file from fp stream descriptor
    and organize its content in list of tuples [(k1,v1), (k2,[(k2.1, v2.1)])]
    If fp is a string, it is used in open() operator
    """
    # check the stream
    open_here=False;
    if isstr(fp):
        fp=codecs.open(fp, "r", encoding="utf-8-sig");
        fp.seek(0);
        open_here=True;
    # error control
    if lev[0] < 0 or indent[0] < 0:
        raise NameError("lev=%d, indent=%d both must be positive"%(lev[0], indent[0]));
    if lev[0] < indent[0]:
        raise NameError("lev=%d, indent=%d, lev must be greater or equal to indent"%(lev[0], indent[0]));
    if fp != sys.stdin and lev[0] > fp.tell():
        raise NameError("lev=%d, file position=%d, lev must be less or equal to file position"%(lev[0], fp.tell()));
    if fp != sys.stdin and indent[0] > fp.tell():
        raise NameError("indent=%d, file position=%d, indent must be less or equal to file position"%(indent[0], fp.tell()));
    # algorithm:
    # advance to requested indent (=level)
    # if not sucsessful return an empty list
    # read a key 
    # if sep==\t read value
    # elif sep==\n
    #     recursive call with increased indentation
    #     if no result at the level+1 put empty value
    # else put empty value
    tlist=[];
    key="";
    val="";
    while True:
        # current position is supposed to point to the begining of a key
        # so go through an appropriate tab number for the current level
        while indent[0] < lev[0]:
            char=fp.read(1);
            if char!="\t":
                if char!="":
                    fp.seek(-1,1);
                break;
            indent[0]+=1;
        if indent[0] < lev[0]:
            # we could not read till the requested level
            # so the current level is finished;
            if open_here:
                fp.close();
            return tlist;
        (key,sep)=kvh_read_key(fp, strip);
        if sep=="\t":
            tlist.append((key, kvh_read_val(fp, strip)));
            indent[0]=0;
        elif sep=="\n":
            lev[0]+=1;
            indent[0]=0;
            nextlist=kvh2tlist(fp, lev, indent, strip);
            lev[0]-=1;
            if len(nextlist)==0:
                # no value and no deeper level
                tlist.append((key, ""));
            else:
                tlist.append((key, nextlist));
        else:
            # we are at the end of file
            if indent[0] or key:
                tlist.append((key, ""));
            indent[0]=0;
            lev[0]=0;
            if open_here:
                fp.close();
            return tlist;

def kvh_read_key(fp, strip=False):
    """Read a string from the current position till the first unescaped \t, \n or the end of stream fp.

:returns: tuple (key, sep), sep=None at the end of the stream

"""
    #pdb.set_trace();##
    key="";
    while True:
        char=fp.read(1);
        if char=="\\":
            # try to read next char if any
            nextchar=fp.read(1);
            if nextchar=="":
                # end of file
                return (key.strip() if strip else key, None);
            else:
                # just add escaped char
                key+=nextchar;
        elif char=="\t" or char=="\n":
            return (key.strip() if strip else key, char);
        elif char=="":
            return (key.strip() if strip else key, None);
        else:
            # just add a plain char
            key+=char;
def kvh_read_val(fp, strip=False):
    """
    Read a string from current position till the first unescaped
    \n or the end of file.
    Return the read string."""
    val="";
    while True:
        char=fp.read(1);
        if char=="\\":
            # try to read next char if any
            nextchar=fp.read(1);
            if nextchar=="":
                # end of file
                return val.strip() if strip else val;
            else:
                # just add escaped char
                val+=nextchar;
        elif char=="\n" or char=="":
            return val.strip() if strip else val;
        else:
            # just add a plain char
            val+=char;
def kvh_tlist2dict(tlist):
    """
    Translate a tlist structure read from a kvh file to
    a hierarchical dictionnary. Repeated keys at the same level
    of a dictionnary are silently overwritten"""
    return dict((k,(v if isstr(v) else kvh_tlist2dict(v))) for (k,v) in tlist);
def kvh_tlist2obj(tlist):
    """
    Translate a tlist structure read from a kvh file to
    a hierarchical dictionnary. Repeated keys at the same level
    of a dictionnary are silently overwritten"""
    return Obj(**dict((k,(v if isstr(v) else kvh_tlist2obj(v))) for (k,v) in tlist));
def kvh2dict(fp, strip=False):
    r"""
    Read a kvh file from fp pointer then translate its tlist
    structure to a returned hierarchical dictionnary.
    Repeated keys at the same level of a dictionnary are
    silently overwritten"""
    return kvh_tlist2dict(kvh2tlist(fp, strip=strip));
def kvh2obj(fp, strip=False):
    r"""
    Read a kvh file from fp pointer then translate its tlist
    structure to a returned object hierarchy.
    Repeated fields at the same level of an object are
    silently overwritten"""
    return kvh_tlist2obj(kvh2tlist(fp, strip=strip));
def dict2kvh(d, fp=sys.stdout, indent=0):
    r"""dict2kvh(d, fp=sys.stdout, indent=0)
    Write a nested dictionary on the stream fp (stdout by default).
    """
    open_here=False;
    if isstr(fp):
        open_here=True;
        fp=open(fp, "w");
    for (k,v) in d.items():
        fp.write("%s%s" % ("\t"*indent, escape(str(k), "\t\\\n")));
        if type(v) == type({}) or type({}) in type(v).__bases__:
            # recursive call with incremented indentation
            fp.write("\n");
            dict2kvh(v, fp, indent+1);
        elif "__dict__" in dir(v) and v.__dict__:
            # recursive call with incremented indentation
            fp.write("\n");
            dict2kvh(v.__dict__, fp, indent+1);
        else:
            fp.write("\t%s\n" % escape(str(v), "\\\n"));
    if open_here:
        fp.close();
def tlist2kvh(d, fp=sys.stdout, indent=0):
    r"""tlist2kvh(d, fp=sys.stdout, indent=0)
    Write a (hierarchichal) list of 2-tuples on the stream fp (stdout by default).
    """
    open_here=False;
    if isstr(fp):
        open_here=True;
        fp=open(fp, "w");
    for (k,v) in d:
        fp.write("%s%s" % ("\t"*indent, escape(str(k), "\t\\\n")));
        if type(v) == type([]):
            # recursive call with incremented indentation
            fp.write("\n");
            tlist2kvh(v, fp, indent+1);
        else:
            fp.write("\t%s\n" % escape(str(v), "\\\n"));
    if open_here:
        fp.close();
def kvh_getv_by_k(kvt, kl):
    r"""kvh_getv_by_k(kvt, kl)->None|String|kvh tlist
    get value from kvt (kvh tlist) according to the key hierarchy
    defined in the list of keys kl. Return None if no key is found
    """
    for (k,v) in kvt:
        if k==kl[0]:
            # found
            if len(kl) == 1:
                return(v);
            elif len(kl) > 1:
                # recursive call
                return(kvh_getv_by_k(v, kl[1:]));
def escape(s, spch="|&;<>()$`\\\"' \t\n*?[#~=%", ech="\\"):
    r"""escape(s, spch="|&;<>()$`\\\"' \t\n*?[#~=%", ech="\\")
escape special characters in s. The special characters are listed in spch.
Escaping is done by putting an ech string before them.
Default spch and ech corresponds to quoting Shell arguments
in accordance with
http://www.opengroup.org/onlinepubs/009695399/utilities/xcu_chap02.html
Example: os.system("ls %s" % escape(file_name_with_all_meta_chars_but_newline));
.. note:

1. Escaped <newline> is removed by a shell if not put in a single-quotted string (' ')

2. A single-quote character even escaped cannot appear in a single-quotted string
"""
    return "".join((ech+c if c in spch else c) for c in s);

def kvh_get_matrix(fp, keys):
    """Get matrix or vector whose key suite is in a list keys from a kvh file given in fp
    (file pointer of file name). For big kvh files, this function can be much faster
    than kvh2tlist()+kvh_getv_by_k()
    Return a matrix which is a list of lists (rows). The first item in each row is the
    row name. In case of matrix (i.e. "row_col" is present in kvh file), the very first row contain column names."""
    
    if isinstance(fp, str):
        with open(fp, "r") as fp:
            cont=fp.readlines()
    else:
        cont=fp.readlines()
    # keys can be a list of subfield keys.
    ncont=len(cont)
    # get start line number by grep successively all fields in v
    # and the indent
    indent=0
    nstart=0
    #pdb.set_trace()
    for k in keys:
        for i in range(nstart, ncont):
            if match(r"^\t{%d,}%s\r?$"%(indent,escape(k)), cont[i]):
                nstart=i+1
                o=match(r"^\t*", cont[i])
                indent=len(o.group(0))
                break
        else:
            raise NameError("The key '%s' was not found in kvh file '%s'."%(k, fp.name))
    # get end number of the matrix row in the kvh
    nend=ncont
    for i in range(nstart, ncont):
        if match(r"^\t{0,%d}[^\t]"%indent, cont[i]):
            nend=i
            break
    # get matrix
    d=[s[(indent+1):].rstrip("\r\n").split("\t") for s in cont[nstart:nend]]
    return d
