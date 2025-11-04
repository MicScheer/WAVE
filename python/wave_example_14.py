
ReplList = [ \
[" CODE='WAVE.EXAMPLE'"," CODE='Example 14: Effective Undulator Source'"], \
[" MYINUM=5000"," MYINUM=2000"], \
[" KELLIP=1"," KELLIP=0"], \
[" KHALBASY=0"," KHALBASY=1"], \
[" IFREQ2P=3"," IFREQ2P=-1"], \
[" IFOLD=1"," IFOLD=0"], \
[" IPHASE=0"," IPHASE=1"], \
[" IEFOLD=1"," IEFOLD=0"], \
[" ZLHALBASY=0.2"," ZLHALBASY=0.056"], \
[" AHWPOL=21."," AHWPOL=161."], \
[" NHHALBASY=0"," NHHALBASY=1"], \
[" HHALBASY=1000."," HHALBASY=300."], \
[" MPINZ=21"," MPINZ=151"], \
[" MPINY=21"," MPINY=151"], \
[" PINW=0.003"," PINW=0.005"], \
[" PINH=0.003"," PINH=0.005"], \
[" FREQLOW=105."," FREQLOW=295."], \
[" FREQHIG=117."," FREQHIG=305."], \
[" NPHASEZ=31"," NPHASEZ=101"], \
[" NPHASEY=31"," NPHASEY=101"], \
[" IWIGNER=0"," IWIGNER=-1"], \
["NWIGTHETAY=31","NWIGTHETAY=1"], \
["NWIGEFOLD=15","NWIGEFOLD=0"], \
]

Nexample = 14


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

