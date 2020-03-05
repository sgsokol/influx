import sys
import re
def cli():
   return(res);
s=sys.argv[0]
#print("s before='"+s+"'")
s=re.sub(r'(-script\.pyw?|\.exe)?$', '', s)
s += ('.py' if s[-3:] != '.py' else '')
#print("s after='"+s+"'")
res=exec(open(s).read())
