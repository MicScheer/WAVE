
ReplList = [ \
[" CODE='WAVE.EXAMPLE'"," CODE='Planar Undulator'"], \
[" KHALBASY=0"," KHALBASY=1"], \
[" PKHALBASY=1"," PKHALBASY=0"], \
[" KELLIP=1"," KELLIP=0"], \
[" ZLHALBASY=0.2"," ZLHALBASY=0.050"], \
[" AHWPOL=21."," AHWPOL=99."], \
[" NHHALBASY=0"," NHHALBASY=3"], \
[" HHALBASY=1000."," HHALBASY=500."], \
[" PINW=0.003"," PINW=0.002"], \
[" PINH=0.003"," PINH=0.002"], \
[" FREQLOW=105."," FREQLOW=495."], \
[" FREQHIG=117."," FREQHIG=505."], \
]
# Note: The line containing PKHALBASY is necessary to reset this variable
# which is overwritten by the line setting KHALBASY.

Nexample = 1


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

