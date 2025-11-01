
Example = '12b'
FileBmap="wave_example12a_bmap.dat"

ReplList = [ \
[" CODE='WAVE.EXAMPLE'"," CODE='Example " + Example + ": Acc. Phys.: TM from B-Map'"], \

[" IUNDULATOR=1"," IUNDULATOR=0"], \
[" IENELOSS=1"," IENELOSS=0"], \
[" MYINUM=5000"," MYINUM=1000"], \

[" IBSUPER=1"," IBSUPER=0"], \
[" IMAGSPLN=-9999"," IMAGSPLN=0"], \
[" KELLIP=1"," KELLIP=0"], \
[" IRFILB0=0"," IRFILB0=-6"], \

[" ISPEC=1"," ISPEC=0"], \
[" ISPECINT=1"," ISPECINT=0"], \
[" IFOLD=1"," IFOLD=0"], \
[" IEFOLD=0"," IEFOLD=0"], \

[" DELTAZ= 0.0"," DELTAZ=0.0001"], \
[" DELTAZP=0.0"," DELTAZP=0.00005"], \
[" DELTAY= 0.0"," DELTAY=0.0001"], \
[" DELTAYP=0.0"," DELTAYP=0.00005"], \
[" DLAPER=5.8"," DLAPER=3.8"], \

["bmap.ntup",FileBmap], \
]

import platform, re
System = platform.system().upper()

if System == 'WINDOWS': ps = "\\"
else: ps = "/"

Fexample = '..' + ps + 'input' + ps + 'wave.example'
FexampleN = '..' + ps + 'input' + ps + 'wave_example_' + Example + '.in'

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

