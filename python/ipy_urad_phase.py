# +PATCH,//WAVE/PYTHON
# +DECK,msh_ipylogon  ,T=PYTHON.

print("\n--- urad_phase ipylogon.py ---\n")
#plt.style.use('seaborn-dark')

import sys,os
pwd = os.getcwd()

import m_hbook as m
from m_hbook import *
from numpy import *

#reakpoint()
args=sys.argv; nargs = len(args)

if args[1] == "last":
  try:
    Farg = open("ipylogon.arg","r")
    argl = Farg.readlines()
    Farg.close()
    args = []
    for arg in argl: args.append(arg.strip())
    #print(args[1:])
  except: pass
#endif

nargs = len(args)

if not args[1] == "last" and nargs > 1:
  Farg = open("ipylogon.arg","w")
  for arg in args: Farg.write(arg + "\n")
  Farg.close()
#endif

seed(0)

NL = '\n'
BL = ' '

C0 = 1+0j
CI = 0+1j

sely0 = 'y==0'
selz0 = 'z==0'
s0sumfld=1.0
s0sumfldh=1.0
s0sumfdp=1.0
s0sumfdph=1.0

ezrmaxfld = 1.0
ezimaxfld = 1.0
s0rmaxfld = 1.0
ezrmaxfdp = 1.0
ezimaxfdp = 1.0
s0rmaxfdp = 1.0
    
ezrmaxfldh = 1.0
ezimaxfldh = 1.0
s0rmaxfldh = 1.0
ezrmaxfdph = 1.0
ezimaxfdph = 1.0
s0rmaxfdph = 1.0
    
ntuples = 0
histos = 1

UradPhase = 0

Tanaka = 1

if Tanaka:
  
  ntau = ncread("ntau","x:y:z:t:betx:bety:betz:tau",'test_tanaka.tau')
  ninfo(ntau)
  
  nfld = ncread("nfld","egam:z:y:erx:ery:erz:eix:eiy:eiz:fd",'test_tanaka.fld')
  ninfo(nfld)
  
  nrad = ncread("nrad","ie:x:y:z:betx:bety:betz:rx:ry:rz:tau:derz:deiz:erz:eiz:",'test_tanaka.rad')
  ninfo(nrad)
  
  nwfd = ncread("nwfd","egam:fd:erz:eiz",'~/wav/stage/wave_rad.dat')
  
  zone(2,1)
  npll(nfld,"egam:fd")
  
  sel="egam<200"
  
  fdmax1 = nfld.query(sel).fd.max()
  ezmax1 = max(nfld.query(sel).erz.max(),nfld.query(sel).eiz.max())
  wfdmax1= nwfd.query(sel).fd.max()
  wezmax1 = max(nwfd.query(sel).erz.max(),nwfd.query(sel).eiz.max())
  
  fdmax2 = nfld.fd.max()
  ezmax2 = max(nfld.erz.max(),nfld.eiz.max())
  wfdmax2= nwfd.fd.max()
  wezmax2 = max(nwfd.erz.max(),nwfd.eiz.max())
  
  #print(fdmax2,wfdmax2,wfdmax2/fdmax2)
  print(ezmax1,g5(wezmax1/ezmax1))
  print(ezmax2,g5(wezmax2/ezmax2))
  #print(ezmax2/ezmax1/(wezmax2/wezmax1))
 
  nex()
  npl(nfld,"egam:erz/"+str(sqrt(fdmax2)))
  xstat(0.7)
  #npllrs(nwfd,"egam:erz/"+str(sqrt(wfdmax2)))
  npllgs(nwfd,"egam:-eiz/"+str(sqrt(wfdmax2)))
#endif Tanaka

  wc()
  npl(nrad,"x:derz","ie==2")

if UradPhase:
  
  set_console_title("Plot urad_phase")
  #reakpoint()
  
  fnam = open("urad_phase.nam","r")
  flines = fnam.readlines()
  fnam.close()
  
  fpin = open("urad_phase.pin",'r')
  
  pin = fpin.readline().strip().split()
  pincen = fpin.readline().strip().split()
  words = fpin.readline().strip().split()
  modepin = int(words[0])
  ifold = int(words[1])
  ifixphase = int(words[2])
  ifieldprop = int(words[3])
  nelec = int(words[4])
  ihbunch=int(words[5])
  words = fpin.readline().strip().split()
  betah = float(words[0])
  emith = float(words[1])
  betav = float(words[2])
  emitv = float(words[3])
  espread = float(words[4])
  words = fpin.readline().strip().split()
  npinzprop = int(words[0])
  npinyprop = int(words[1])
  pinxprop = float(words[2])
  pinwprop = float(words[3])
  pinhprop = float(words[4])
  
  fpin.close()
  
  npinz = int(pin[0])
  npiny = int(pin[1])
  pinw = float(pin[2])
  pinh = float(pin[3])
  
  pinx = float(pincen[0])
  piny = float(pincen[1])
  pinz = float(pincen[2])
  
  ymin = piny - pinh/2.
  ymax = piny + pinh/2.
  
  zmin = pinz - pinw/2.
  zmax = pinz + pinw/2.
  
  if npiny > 1:
    dy = (ymax-ymin)/(npiny-1)
  else:
    dy = pinh
  #endif
  
  if npinz > 1:
    dz = (zmax-zmin)/(npinz-1)
  else:
    dz = pinw
  #endif
  
  if (fexist('urad_phase.fld')):
    nfld = ncread("nfld","x:y:z:iegam:egam:s0:s1:s2:s3:p:exr:exi:eyr:eyi:ezr:ezi:bxr:bxi:byr:byi:bzr:bzi:nx:ny:nz","urad_phase.fld")
    s0sumfld=nsum(nfld,"s0",isilent=1)
    s0sumfldh=nsum(nfld,"s0",sely0,isilent=1)
    ezrmaxfld=nmax(nfld,"ezr")
    ezrmaxfldh=nmax(nfld,"ezr",sely0)
    ezimaxfld=nmax(nfld,"ezi")
    ezimaxfldh=nmax(nfld,"ezi",sely0)
    s0maxfld=nmax(nfld,"s0")
    s0maxfldh=nmax(nfld,"s0",sely0)
  else:
    Quit("Keine Daten!")
  #endif
  
  if (fexist('urad_phase.fdf')):
    nfdf = ncread("nfdf","iegam:iy:iz:egam:x:y:z:exr:exi:eyr:eyi:ezr:ezi:bxr:bxi:byr:byi:bzr:bzi:s0:s1:s2:s3","urad_phase.fdf")
  if (fexist('urad_phase.dum')):
    ndum = ncread("ndum","iegam:iy:iz:egam:x:y:z:exr:exi:eyr:eyi:ezr:ezi:bxr:bxi:byr:byi:bzr:bzi:s0:s1:s2:s3","urad_phase.dum")
  if (fexist('urad_phase.fdpf')):
    nfdpf = ncread("nfdpf","iegam:iy:iz:egam:x:y:z:exr:exi:eyr:eyi:ezr:ezi:bxr:bxi:byr:byi:bzr:bzi:s0:s1:s2:s3","urad_phase.fdpf")
  if (fexist('urad_phase.geo')):
    ngeo = ncread("ngeo","xo:yo:zo:iegam:egam:s0:s1:s2:s3:nx:ny:nz:xg:yg:zg","urad_phase.geo")
          
  if (fexist('urad_phase.bun')):
    nbun = ncread("nbun","jbun:isub:ibu:bunchx:rxi1:ryi1:rzi1:ypi1:zpi1:rxin:ryin:rzin:ypin:zpin:eel:deel:x:y:z:iegam:egam:spec:s0:s1:s2:s3:p:fb28:dt:exr:exi:eyr:eyr:ezr:ezi:bxr:bxi:byr:byi:bzr:bzi","urad_phase.bun")
    nelec = nbun.ibu.max()
  #endif
      
  if (fexist('urad_phase.fdp')):
        
    nfdp = ncread("nfdp","x:y:z:iegam:egam:s0:s1:s2:s3:exr:exi:eyr:eyi:ezr:ezi:bxr:bxi:byr:byi:bzr:bzi:nx:ny:nz","urad_phase.fdp")
    
    s0sumfdp=nsum(nfdp,"s0",isilent=1)
    s0sumfdph=nsum(nfdp,"s0",sely0,isilent=1)
    ezrmaxfdp=nmax(nfdp,"ezr")
    ezrmaxfdph=nmax(nfdp,"ezr",sely0)
    ezimaxfdp=nmax(nfdp,"ezi")
    ezimaxfdph=nmax(nfdp,"ezi",sely0)
    s0maxfdp=nmax(nfdp,"s0")
    s0maxfdph=nmax(nfdp,"s0",sely0)
    
    lfdp = len(nfdp)
    negam = nfdp.iegam.max()
    nzy = int(lfdp/negam)
    
    yminph = -pinhprop/2.0
    ymaxph = -yminph
    
    if npinyprop > 1:
      dyph = (ymaxph-yminph)/(npinyprop-1)
    else:
      dyph = pinhprop
    #endif
    
    zminph = -pinwprop/2.0
    zmaxph = -zminph
    if npinzprop > 1:
      dzph = (zmaxph-zminph)/(npinzprop-1)
    else:
      dzph = pinwprop
    #endif
    
    ninfo(nfdp)
          
  #endif
      
  if (fexist('urad_phase.wig')):
    nwig = ncread("nwig","iz:iy:itz:ity:x:y:z:ty:tz:iegam:egam:wzzr:wzzi:wzyr:wzyi:wyzr:wyzi:wyyr:wyyi","urad_phase.wig")
    nywig = nwig.iy.max()
    yminwig = nwig.y.min()
    ymaxwig = nwig.y.max()
    if nywig > 1:
      dywig = (ymaxwig-yminwig)/(nywig-1)
    else:
      dywig = 1.0
    #endif
    yminwig -= dywig/2.0
    ymaxwig += dywig/2.0
    nzwig = nwig.iz.max()
    zminwig = nwig.z.min()
    zmaxwig = nwig.z.max()
    if nzwig > 1:
      dzwig = (zmaxwig-zminwig)/(nzwig-1)
    else:
      dzwig = 1.0
    #endif
    if nywig > 1:
      dywig = (ymaxwig-yminwig)/(nywig-1)
    else:
      dywig = 1.0
    #endif
    zminwig -= dzwig/2.0
    zmaxwig += dzwig/2.0
    nty = nwig.ity.max()
    tymin = nwig.ty.min()
    tymax = nwig.ty.max()
    if nty > 1:
      dty = (tymax-tymin)/(nty-1)
    else:
      dty = 1.0
    #endif
    tymin -= dty/2.0
    tymax += dty/2.0
    ntz = nwig.itz.max()
    tzmin = nwig.tz.min()
    tzmax = nwig.tz.max()
    if ntz > 1:
      dtz = (tzmax-tzmin)/(ntz-1)
    else:
      dtz = 1.0
    #endif
    tzmin -= dtz/2.0
    tzmax += dtz/2.0
    ninfo(nwig)
  #endif
  
  if (fexist('urad_phase.ebm')):
    nebm = ncread("nebm","iegam:iz:iy:egam:z:y:exr:exi:eyr:eyi:ezr:ezi:bxr:bxi:byr:byi:bzr:bzi:nen:s0:s1:s2:s3:nx:ny:nz","urad_phase.ebm")
  #endif
  
  nlist()
  
  if args[1] == 'none':
    pass
  
  elif args[1] == 'hcfluxden': #key hcfluxden

    lolo()
    npll(nfld,"egam:s0")
    txyz("S0, nfld","E [eV]")
    
  elif args[1] == 'default': #key default
    
    zone(3,2)
    
    npll(nfld,"z:s0",sely0)
    txyz("S0, nfld," + sely0,"z[mm]")
    
    nex()
    optnstat()
    
    npll(nfld,"z:ezr",sely0)
    npllgs(nfld,"z:ezi",sely0)
    txyz("Ez, nfld," + sely0,"z[mm]")
    
    optstat()
    nex()
    fn = str('/') + str(s0sumfldh)
    npll(nfld,"z:nz*1000*s0"+fn,sely0)    
    txyz("nz, nfld, norm." + sely0,"z[mm]")
    
    nex()
    
    npll(nfdp,"z:s0",sely0)
    txyz("S0, nfdp," + sely0,"z[mm]")
    
    nex()
    optnstat()
    
    fn = str('/') + str(ezrmaxfdph)
    npll(nfdp,"z:ezr"+fn,sely0)
    npllgs(nfdp,"z:ezi"+fn,sely0)    
    txyz("Ez, nfdp, norm." + sely0,"z[mm]")
    
    nex()
    optstat()
    fn = str('/') + str(s0sumfdph)
    npll(nfdp,"z:nz*1000*s0"+fn,sely0)
    txyz("nz, nfdp, norm." + sely0,"z[mm]")
    
  #endif args[1]
  
#endif UradPhase

if not ntuples:
    for n in range(len(Nhead)):
        snam = Nhead[n][1]
        exec(snam + ' = nget("' + snam + '")')
    #endfor
#endif

if histos:
  for hh in H1head:
    snam = hh[0]
    exec(snam + ' = hget("' + snam + '")')
  #endfor
  for hh in H2head:
    snam = hh[0]
    exec(snam + ' = hget("' + snam + '")')
  #endfor
#endif

#wans()
