
ReplList = [ \
[" CODE='WAVE.EXAMPLE'"," CODE='Example 8: Magnetic Field Errors'"], \
[" KELLIP=1"," KELLIP=0"], \
[" KBREC=0"," KBREC=1"], \
[" IRECU=0"," IRECU=1"], \
[" KRECPER\(1\)=30"," KRECPER(1)=100"], \
[" URECLX\(1\)=14."," URECLX(1)=7.5"], \
[" URECLZ\(1\)=40."," URECLZ(1)=80."], \
[" URECGAP\(1\)=11.6"," URECGAP(1)=2.5"], \
[" IUHELI\(1\)=1"," IUHELI(1)=0"], \
[" NURANMOD=1"," NURANMOD=0"], \
[" IPIN=1"," IPIN=0"], \
[" IFOLD=1"," IFOLD=0"], \
[" IEFOLD=1"," IEFOLD=0"], \
[" FREQLOW=105."," FREQLOW=1100."], \
[" FREQHIG=117."," FREQHIG=1110."], \
[" NINTFREQ=13"," NINTFREQ=201"], \
]

Nexample = 8


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

