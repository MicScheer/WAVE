
# Example 11 MUST be setup first!!

Example = '11a'


ReplList = [ \
[" XSHIFT=0.0"," XSHIFT=0.056"], \
[" KBUNDUMAG=2"," KBUNDUMAG=1"], \
[" IPERIODG=0"," IPERIODG=-1"], \
[" XPERWMX=0."," XPERWMX=1.12"], \
[" PERIODG=0.8"," PERIODG=0.056"], \
[" PEROFFG=-0.4"," PEROFFG=-0.028"], \
]


Fexample = 'wave.in'
FexampleN = 'pex11a.in'

Fex = open(Fexample,'r')
Fout = open(FexampleN,'w')

import re

for line in Fex:
  for i in range(len(ReplList)):
    rep = ReplList[i]
    lino = re.sub(rep[0],rep[1],line)
    if lino != line:
      print("\n----\n" + line)
      print(lino)
      break
  #endfor i in len(range(ReplList))
  Fout.write(lino)
#endfor line in Fex

Fout.close()
Fex.close()

