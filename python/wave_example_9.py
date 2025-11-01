
ReplList = [ \
[" CODE='WAVE.EXAMPLE'"," CODE='Example 9: Real Beams, Part II'"], \
[" IUNDULATOR=1"," IUNDULATOR=0"], \
[" NLPOI=-9999"," NLPOI=0"], \
[" DMYENERGY=1.722"," DMYENERGY=17."], \
[" IENELOSS=-1"," IENELOSS=1"], \
[" XSTART=9999."," XSTART=0."], \
[" XSTOP=9999."," XSTOP=180."], \
[" IHTRSMP=1"," IHTRSMP=10"], \
[" KMAGSEQ=0"," KMAGSEQ=1"], \
[" IMAGSPLN=-9999"," IMAGSPLN=0"], \
[" KELLIP=1"," KELLIP=0"], \
[" KMAGSPLN=-9999"," KMAGSPLN=0"], \
[" WGWINFC=9999."," WGWINFC=100."], \
[" ISPECINT=1"," ISPECINT=0"], \
[" IPBRILL=0"," IPBRILL=1"], \
[" IBRILL=1"," IBRILL=0"], \
[" PINW=0.003"," PINW=0.00007"], \
[" PINH=0.003"," PINH=0.00007"], \
[" IFOLD=1"," IFOLD=0"], \
[" IEFOLD=1"," IEFOLD=0"], \
[" IPERIODG=0"," IPERIODG=1"], \
[" PERIODG=0.8"," PERIODG=6."], \
[" PEROFFG=-0.4"," PEROFFG=0.0"], \
[" BETFUN=8.71"," BETFUN=14.3"], \
[" BETAH=9999."," BETAH=20.2"], \
[" BETAPH=9999."," BETAPH=0.0"], \
[" EPS0H=7.70e-9"," EPS0H=4.0e-11"], \
[" BETFUNV=4.36"," BETFUNV=20.2"], \
[" BETAV=9999."," BETAV=14.3"], \
[" BETAPV=9999."," BETAPV=0.0"], \
[" EPS0V=15.4e-11"," EPS0V=4.0e-11"], \
[" BSIGZ\(1\)=275.e-6"," BSIGZ(1)=28.4e-6"], \
[" BSIGZP\(1\)=28.1e-6"," BSIGZP(1)=1.41e-6"], \
[" BSIGY\(1\)=22.5e-6"," BSIGY(1)=23.9e-6"], \
[" BSIGYP\(1\)=6.80e-6"," BSIGYP(1)=1.67e-6"], \
[" ESPREAD=0.001"," ESPREAD=0.0001"], \
[" PINCEN\(1\)=10."," PINCEN(1)=250."], \
[" FREQLOW=105."," FREQLOW=84990."], \
[" FREQHIG=117."," FREQHIG=85010."], \
[" NINTFREQ=13"," NINTFREQ=21"], \
[" IUBUNCH=0"," IUBUNCH=1"], \
]

Nexample = 9


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

