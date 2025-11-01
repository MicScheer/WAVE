
ReplList = [ \
[" CODE='WAVE.EXAMPLE'"," CODE='Example 11: Hybrid-Undulator with UNDUMAG'"], \
[" IUNDULATOR=1"," IUNDULATOR=0"], \
[" KELLIP=1"," KELLIP=0"], \
[" KBUNDUMAG=0"," KBUNDUMAG=4"], \
[" IMAGSPLN=-9999"," IMAGSPLN=0"], \
[" ISPEC=1"," ISPEC=0"], \
]

Nexample = 13


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

