import sys
s=sys.argv[0]
s=(s[-4:] if s[-4:].lower == ".exe" else s)+".py"
res=exec(open(s).read())
def cli():
   return(res);
