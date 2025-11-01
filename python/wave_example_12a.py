
Example = '12a'

FullGap = 11.4
HalfGap = FullGap / 2.

FileBmap="wave_example12a_bmap.dat"

import os
try: os.remove(FileBmap)
except: pass

ReplList = [ \
[" CODE='WAVE.EXAMPLE'"," CODE='Ex. 12a: Acc. Phys.: Transfer-mat., Multipoles, B-Map'"], \

[" IUNDULATOR=1"," IUNDULATOR=0"], \

[" MYINUM=5000"," MYINUM=1000"], \
[" IENELOSS=1"," IENELOSS=0"], \
[" XSTART=9999."," XSTART=-1.9"], \
[" XSTOP=9999."," XSTOP=1.9"], \

[" IMAGSPLN=-9999"," IMAGSPLN=0"], \
[" IBSUPER=1"," IBSUPER=0"], \
[" KELLIP=1"," KELLIP=0"], \
[" KBREC=0"," KBREC=1"], \
[" IRECU=0"," IRECU=1"], \

[" URECBC\(1\)=1.22"," URECBC(1)=1.3"], \
[" URECGAP\(1\)=11.6"," URECGAP(1)=" + str(HalfGap)], \
[" URSHIFT\(1\)=10."," URSHIFT(1)=0.0"], \
[" NURANMOD=1"," NURANMOD=0"], \

[" ISPEC=1"," ISPEC=0"], \
[" ISPECINT=1"," ISPECINT=0"], \
[" IFOLD=1"," IFOLD=0"], \
[" IEFOLD=0"," IEFOLD=0"], \

[" DELTAZ= 0.0"," DELTAZ=0.0001"], \
[" DELTAZP=0.0"," DELTAZP=0.00005"], \
[" DELTAY= 0.0"," DELTAY=0.0001"], \
[" DELTAYP=0.0"," DELTAYP=0.00005"], \
[" DLAPER=5.8"," DLAPER=3.8"], \

[" IWBTAB=0"," IWBTAB=1"], \
[" IWBMAP=0"," IWBMAP=6"], \

[" YMAPMN=-0.007"," YMAPMN=-" +  str((HalfGap-0.5)/1000.)], \
[" YMAPMX=\+0.007"," YMAPMX=+" + str((HalfGap-0.5)/1000.)], \
[" NMAPY=3"," NMAPY=11"], \
[" NMAPZ=31"," NMAPZ=61"], \

["bmap.dat",FileBmap], \
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
    #endif
  #endfor i in len(range(ReplList))
  Fout.write(lino)
#endfor line in Fex

Fout.close()
Fex.close()

print(FexampleN + " is ready now. Please copy it to wave.in to use it")

