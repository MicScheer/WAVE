
ReplList = [ \
[" CODE='WAVE.EXAMPLE'"," CODE='Example 3: The Concept of Source Points'"], \
[" MYINUM=5000"," MYINUM=100000"], \
[" KELLIP=1"," KELLIP=0"], \
[" KHALBASY=0"," KHALBASY=1"], \
[" B0HALBASY=2."," B0HALBASY=6."], \
[" PKHALBASY=1"," PKHALBASY=0"], \
[" ZLHALBASY=0.2"," ZLHALBASY=0.7"], \
[" AHWPOL=21."," AHWPOL=1."], \
[" IPIN=1"," IPIN=0"], \
[" IFOLD=1"," IFOLD=0"], \
[" IEFOLD=1"," IEFOLD=3"], \
[" IBRILL=1"," IBRILL=0"], \
[" FREQLOW=105."," FREQLOW=50."], \
[" FREQHIG=117."," FREQHIG=25000."], \
[" NINTFREQ=13"," NINTFREQ=2001"], \
]

# Note: The line containing PKHALBASY is necessary to reset this variable
# which is overwritten by the line setting KHALBASY.

Nexample = 3


import platform, re
System = platform.system().upper()

if System == 'WINDOWS': ps = "\\"
else: ps = "/"

Fexample = '..' + ps + 'input' + ps + 'wave.example'
FexampleN = '..' + ps + 'input' + ps + 'wave_example_' + str(Nexample) + '.in'

Fex = open(Fexample,'r')
Fout = open(FexampleN,'w')

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

print(FexampleN + " is ready now. Please copy it to wave.in to use it")

