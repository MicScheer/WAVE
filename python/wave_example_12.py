
ReplList = [ \
[" CODE='WAVE.EXAMPLE'"," CODE='Example 12: Accelerator Physics'"], \
[" IUNDULATOR=1"," IUNDULATOR=0"], \
[" KELLIP=1"," KELLIP=0"], \
[" KBREC=0"," KBREC=1"], \
[" NURANMOD=1"," NURANMOD=0"], \
[" ISPEC=1"," ISPEC=0"], \
[" ISPECINT=1"," ISPECINT=0"], \
[" IFOLD=1"," IFOLD=0"], \
[" IEFOLD=0"," IEFOLD=0"], \
[" IENELOSS=1"," IENELOSS=0"], \
[" IMAGSPLN=-9999"," IMAGSPLN=0"], \
[" DELTAZ=0.0"," DELTAZ=0.001"], \
[" DELTAZP=0.0"," DELTAZP=0.001"], \
[" DELTAY=0.0"," DELTAY=0.001"], \
[" DELTAYP=0.0"," DELTAYP=0.001"], \
[" DLAPER=5.8"," DLAPER=3.8"], \
[" IOPTIC=0"," IOPTIC=1"], \
[" IGENFUN=0"," IGENFUN=1"], \
]

Example = '12a'


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

