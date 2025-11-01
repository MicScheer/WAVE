
ReplList = [ \
[" CODE='WAVE.EXAMPLE'"," CODE='Example 2: Tracking through Magnetic Fields'"], \
[" BYGOFF=0.0"," BYGOFF=0.00005"], \
[" KELLIP=1"," KELLIP=0"], \
[" KHALBASY=0"," KHALBASY=1"], \
[" B0HALBASY=2."," B0HALBASY=1."], \
[" PKHALBASY=1"," PKHALBASY=0"], \
[" AHWPOL=21."," AHWPOL=99."], \
[" ISPEC=1"," ISPEC=0"], \
]

# Note: The line containing PKHALBASY is necessary to reset this variable
# which is overwritten by the line setting KHALBASY.

Nexample = 2


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

