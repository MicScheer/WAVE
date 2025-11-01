
ReplList = [ \
[" CODE='WAVE.EXAMPLE'"," CODE='Global Magnetic Field Options '"], \
[" DMYENERGY=1.722"," DMYENERGY=17."], \
[" IENELOSS=1"," IENELOSS=0"], \
[" XSTART=9999."," XSTART=-1.5"], \
[" XSTOP=9999."," XSTOP=1.5"], \
[" IHTRSMP=1"," IHTRSMP=10"], \
[" KHALBASY=0"," KHALBASY=1"], \
[" ISPECINT=1"," ISPECINT=0"], \
[" IPIN=1"," IPIN=0"], \
[" IBRILL=1"," IBRILL=0"], \
[" IFOLD=1"," IFOLD=0"], \
[" IEFOLD=1"," IEFOLD=0"], \
[" IPERIODG=0"," IPERIODG=1"], \
[" PERIODG=0.8"," PERIODG=3."], \
[" PEROFFG=-0.4"," PEROFFG=-1.5"], \
[" B0ELLIPH=0.353553390593"," B0ELLIPH=0.0"], \
[" B0ELLIPV=0.353553390593"," B0ELLIPV=0.45872917839670901d0"], \
[" XLELLIP=0.056"," XLELLIP=0.05"], \
[" PERELLIP=31."," PERELLIP=50."], \
[" PKHALBASY=1"," PKHALBASY=0"], \
[" B0HALBASY=2."," B0HALBASY=4.0501305277737230D-02"], \
[" AHWPOL=21."," AHWPOL=1."], \
[" XCENHAL=0.0"," XCENHAL=1.4"], \
[" PINCEN\(1\)=10."," PINCEN(1)=250."], \
[" FREQLOW=105."," FREQLOW=149900."], \
[" FREQHIG=117."," FREQHIG=150100."], \
[" NINTFREQ=13"," NINTFREQ=201"], \
]

# Note: The line containing PKHALBASY is necessary to reset this variable
# which is overwritten by the line setting KHALBASY.

Nexample = 7


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

