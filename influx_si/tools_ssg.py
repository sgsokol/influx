import sys
import operator
import string
import numpy as np
from math import *
from kvh import *
letters=string.ascii_lowercase

def aff(name, obj, ident=0, f=sys.stdout):
    r"""
    print formatted object:
     name=obj
    """
    saveout=sys.stdout
    sys.stdout=f
    to=type(obj)
    if to == type({}):
        # dictionary case
        print('|'*ident+'+'+str(name)+' #'+str(len(obj)))
        for k in obj:
            aff('{'+str(k)+'}', obj[k], ident+1, f)
    elif to == type(()):
        # tuple case
        print('|'*ident+'+'+str(name)+' #'+str(len(obj)))
        for (i,k) in enumerate(obj):
            aff('('+str(i)+')', obj[i], ident+1, f)
    elif to == type([]):
        # list case
        print('|'*ident+'+'+str(name)+' #'+str(len(obj)))
        for (i,k) in enumerate(obj):
            aff('['+str(i)+']', obj[i], ident+1, f)
    else:
        print('%s%s: %s' % ('|'*ident+'+', name, obj))
    sys.stdout=saveout

def asort(d):
    r"""sorts a dictionnary by value preserving key=value association the result is a list of tuples (key,value)"""
    return sorted(list(a.items()), key=operator.itemgetter(1))

def iterbit(i, size=0):
    r"""iterator on bits in integer starting from 0-position. The iterator stops at highest non-zero bit"""
    i=int(i)
    moveb=1
    b_no=0
    while (moveb <= i and size==0) or (b_no < size):
        yield (1 if (moveb&i) else 0)
        moveb<<=1
        b_no+=1

def iternumbit(i, size=0):
    r"""iterator on bits and its number in integer starting from 0-position. The iterator yields tuples (n,bit). If optional size is zero then it stops at highest non-zero bit. If not, it will stop at bit number size-1."""
    i=int(i)
    moveb=1
    b_no=0
    while (moveb <= i and size==0) or (b_no < size):
        yield (b_no, 1 if (moveb&i) else 0)
        moveb<<=1
        b_no+=1

def sumbit(i):
    r""":returns: sum of bits in an integer"""
    return sum(iterbit(i))

def strbit32(i):
    r""":returns: a string of 0-1s (in chunk of 4) in an 32 bit integer"""
    i=int(i)
    moveb=1<<31
    res=''
    for b in range(32):
        res+=str(int((moveb&i) > 0))
        res+=' ' if not (b+1)%4 and b < 31 else ''
        moveb>>=1
    return res
def strbit(i,size=0):
    r""":returns: the lowest part of integer as string binary representation"""
    return ''.join(('1' if b else '0') for b in iterbit(i,size) )[::-1]
def rstrbit(i,size=0):
    r""":returns: the integer as reversed string binary representation. The lowest bit is on the left side"""
    return ''.join(('1' if b else '0') for b in iterbit(i,size) )

def setbit32(i, nb):
    r"""set a bit number nb (0 based) in an integer i"""
    if (nb < 0 or nb > 31):
       return i; # do nothing to i
    i|=(1<<nb)

def valval(o, keepNone=True):
    r""":returns: an iterator over values of values, i.e. collapsing values of fisrt two nested lists in one list, for example."""
    if isstr(o):
        yield o
    else:
        for i1 in o:
            if isstr(i1):
                yield i1
            else:
                try:
                    for i2 in i1:
                        if keepNone:
                            yield i2
                        elif i2 is not None:
                            yield i2
                            
                except:
                    if keepNone:
                        yield i1
                    elif i1 is not None:
                            yield i1

def join(c,l,p='',s='',a=''):
    r"""join the items of the list (or iterator) l separated by c. Each item is prefixed with p and suffixed with s. If the join result is empty for any reason, an alternative a is returned. p, s and a are optional"""
    i=0
    return c.join(p+str(i)+s for i in l) or a
def joint(c,l,p='',s='',a=''):
    r"""join "true" items of the list (or iterator) l separated by c. Each item is prefixed with p and suffixed with s. If the join result is empty for any reason, an alternative a is returned. p, s and a are optional"""
    i=0
    return c.join(p+str(i)+s for i in l if i) or a

def list2count(l, incr=1):
    r"""count values in a (short) list l incrementing the counter by optional incr.
    
    :returns: a dictionary {item:count}"""
    dico=dict((x,l.count(x)*incr) for x in set(l))
    return dico

def strbit2int(s):
    r"""translate a string of 0's and 1's interpreted as bits to an integer all characters different from 0,1 are silently ignored."""
    res=0
    movb=1
    for c in s[::-1]:
        if c=='1':
            res+=movb
        elif c!='0':
            continue
        movb<<=1
    return res

def setcharbit(s,ch,i):
    r"""set character ch in a string s everywhere a corresponding bit of i is set"""
    res=''.join((ch if b else s[-b_no-1]) for (b_no,b) in iternumbit(i))
    #print 'in:', s, ch, strbit(i)
    #print 'r', res
    res=s[:len(s)-len(res)]+res[::-1]
    #print 'final', res
    return res

def expandbit(i,pos):
    r"""copy bits set to 1 in i to the result position given in the list pos. length of pos must be greater or equal to bitlength of i"""
    #print i, pos;##
    return sum(b<<pos[b_no] for (b_no,b) in iternumbit(i))

def isstr(s):
    r""":returns: True if the argument is a string"""
    return isinstance(s,str)

def trd(l, d, p="", s="", a=""):
    r"""translate items in an iterable l by a dictionary d, prefixing translated items by optional p and suffixing them by optional s. If an item is not found in the dictionnary alternative string a is used. If a==None, the item is left unchanged. No prefix or suffix are applied in both case.

    :returns: iterator"""
    return (p+str(d[i])+s if i in d else i if a==None else a for i in l)
def icumo2iiso(icumo, size):
    r""":returns: iterator on isotopomers composing a given icumo. size is carbon number"""
    nzero=size-sumbit(icumo)
    zpos=[n for (n,b) in iternumbit(icumo, size) if not b]
    return (expandbit(iz,zpos)+icumo for iz in range(1<<nzero))
def reverse(it):
    r"""reverse order of an iterable"""
    for i in range(len(it),0,-1):
        yield it[i-1]
def ssign(i, sp="+", sm="-"):
    r""":returns: a string of i sign: sp (i>=0) or sm (i<0). """
    return sp if float(i) >= 0 else sm
def arr2pbm(A, fp):
    r"""Write an image map of non-zero entries of matrix A to file pointer fp. Matrix A is an array"""
    fp.write("P1\n%d %d\n" % A.shape[::-1])
    for row in A:
        p=0
        for c in row:
            fp.write("1 " if c else "0 ")
            p+=2
            if p >= 69:
                fp.write("\n")
                p=0
        fp.write("\n")
        p=0
def cumsum(l, tot=0):
   r""":Returns: an iterable of the length len(l)+1 with cumulated sum of items in l. First element in cumsum is equal to initial value of tot. Result depends on the meaning of "+" operator for l items and of tot type.
   
.. doctest:

>>> list(cumsum("abc",tot=""))
['', 'a', 'ab', 'abc']

>>> list(cumsum(xrange(1,5)))
[0, 1, 3, 6, 10]
   """
   for i in l:
       yield tot
       tot+=i
   yield tot
def wxlay2py(kvt, parent=[None]):
    r"""
    :returns: a string with python code generating wxWindow widget layout described in kvh tlist sturcture
    """
    res=""
    # get the kvh tuples
    for (k,v) in kvt:
        if k[:3]=="wx." or k[:3]=="wx_"and type(v)==type([]):
            ## produce the code for this widget
            # call class init
            varname=kvh_getv_by_k(v, ["varname"])
            res+=(varname+"=") if varname else ""
            param=kvh_getv_by_k(v, ["param"])
            if param:
                param=param.replace(".parent", str(parent[-1]))
                if len(parent) > 1:
                    param=param.replace(".top", str(parent[1]))
            res+=k+"("+(param or "")+");\n"
            # recursivly create children
            res+=wxlay2py(v, parent+[varname])
            # call methods if any and varname is set
            cl=kvh_getv_by_k(v, ["callmeth"])
            if varname and cl:
                for (kc,vc) in cl:
                    if vc:
                        # if key and value store the result in key
                        vc=vc.replace(".parent", str(parent[-1]))
                        vc=vc.replace(".self", varname)
                        kc=kc.replace(".self", varname)
                        if len(parent) > 1:
                            vc=vc.replace(".top", str(parent[1]))
                        res+=kc+"="+(varname+"."
                            if vc[:2] != "wx" else "")+vc+";\n"
                    elif kc:
                        # just call the key content
                        if kc[:1] != "#":
                            kc=kc.replace(".parent", str(parent[-1]))
                            kc=kc.replace(".self", varname)
                            if len(parent) > 1:
                                kc=kc.replace(".top", str(parent[1]))
                            res+=(varname+"." if vc[:2] != "wx"
                                else "")+kc+";\n"
        else:
            # we don't know what it is, just silently ignore
            pass
    return(res)

def ulong(i):
    r"""ulong(i) -> workarounded ulong"""
    return(-((i^0xffffffff)+1))

def read_table(fp, sep=None, skip=0, header=False):
    r"""read_table(f) -> dict(mat, col_names) read a plain text file f in a numpy mat. If some columns are not numerical, they are replaced by np.nan. If header=True, number of column names in the first row after skip must be the same as the number of values in each following row.
    """
    open_here=False
    if isstr(fp):
        fp=open(fp, "rb")
        fp.seek(0)
        open_here=True
    i=0
    while i < skip:
        fp.readline()
        i+=1
        continue
    line=fp.readline()
    fields=line.strip().split(sep)
    # count fields, create array
    col_names=fields if header else None
    # proceed normal (i.e. data) line
    data=np.asmatrix(np.ndfromtxt(fp, delimiter=sep))
    if open_here:
        fp.close()
    return {"mat": data, "col_names": col_names}
