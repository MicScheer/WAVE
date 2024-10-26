# +PATCH,//WAVE/PYTHON
# +DECK,msh_ipylogon  ,T=PYTHON.

print("\n--- urad_phase ipylogon.py ---\n")
#plt.style.use('seaborn-dark')

import sys,os

args=sys.argv; nargs = len(args)
pwd = os.getcwd()

import m_hbook as m
from m_hbook import *
from numpy import *

#reakpoint()
args=sys.argv; nargs = len(args)
seed(0)

NL = '\n'
BL = ' '

ntuples = 0
histos = 1

global Irunmin,Irunmax

def get_phase_error(run):
  Fout = open("a/wave.out." + run,"r")
  fout = Fout.readlines()
  ifound = 0
  for line in fout:
    if ifound:
      phwav = float(line.split()[0])
      phdeg = float(line.split()[1])
      Fout.close()
      return phwav,phdeg
    #endif
    if re.search('deg.\)',line):
      ifound = 1
    #endif
  #endfor
#enddef

def ns0peak(nt):

  iemax = -1
  emax = 0.0
  s0max = -1.0

  res = npeaks(nt,"e:s0/1.0e6",isilent=1,iretval=1)

  try:
    s0list = res[2]
    s0max = 0.0
    for i in range(len(s0list)):
      s0 = s0list[i]
      if s0 > s0max:
        s0max = s0
        is0 = i
      #endif
    #endfor
    iemax = res[0][i]
    emax = res[1][i]
  except:
    print("*** Warning in ns0peak: Bad return from npeaks ***\n")
#    ninfo(nt)
#    Quit("\n*** Aborted in ns0peak(nt)")
  #endtry
  return iemax,emax,s0max
#enddef

def phase_error( key = "stokes", \
fs0="wave_stokes_selected.dat",\
fs0f="wave_stokesf_selected.dat", \
fs0e="wave_stokese_selected.dat", \
fs0ef="wave_stokesef_selected.dat", \
fs0ref="/home/scheer/spectra/wav/a/wave_halbach_pencil_on-axis.dat.1060", \
fs0reff="/home/scheer/spectra/wav/a/wave_halbach_emit_on-axis.dat.949", \
fs0refe="/home/scheer/spectra/wav/a/wave_halbach_espread_on-axis.dat.950", \
fs0refef="/home/scheer/spectra/wav/a/wave_halbach_emit_espread_on-axis.dat.951", \
fall="SRI22_phase-errors.dat"
):

#fs0ref="/home/scheer/spectra/wav/a/wave_halbach_pencil_on-axis.dat.476", \
#fs0reff="/home/scheer/spectra/wav/a/wave_halbach_emit_on-axis.dat.477", \
#fs0refe="/home/scheer/spectra/wav/a/wave_halbach_espread_on-axis.dat.478", \
#fs0refef="/home/scheer/spectra/wav/a/wave_halbach_emit_espread_on-axis.dat.479", \

  global Nhead,Ntup,Nind,Nntup,Irunmin,Irunmax

  nl = "\n"
  print("Ref:",nl,fs0ref,nl,fs0refe,nl,fs0reff,nl,fs0refef,nl)
  print(key,nl,fs0,nl,fs0e,nl,fs0f,nl,fs0ef)

  Irunmin = 1e20
  Irunmax = -1e20

  iruns0 = 0
  iruns0e = 0
  iruns0f = 0
  iruns0ef = 0
  iruns0ref = 0

  vs0ene = None
  vs0 = None
  vs0f = None
  vs0e = None
  vs0ef = None

  optnstat()

  if nexist("ns0"): ndelete("ns0")
  if nexist("ns0e"): ndelete("ns0e")
  if nexist("ns0f"): ndelete("ns0f")
  if nexist("ns0ef"): ndelete("ns0ef")

  if nexist("ns0ref"): ndelete("ns0ref")
  if nexist("ns0refe"): ndelete("ns0refe")
  if nexist("ns0reff"): ndelete("ns0reff")
  if nexist("ns0refef"): ndelete("ns0refef")

  ns0ref = ncread("ns0ref","e:s0:s1:s2:s3",fs0ref,silent=1,skiphead=2)
  iemaxr,emaxr,s0ref = ns0peak(ns0ref)
  vs0refene = ns0ref['e']
  vs0ref = ns0ref['s0']/1.e6
  nener = len(vs0ref)

  ns0reff = ncread("ns0reff","e:s0:s1:s2:s3",fs0reff,silent=1,skiphead=2)
  iemaxrf,emaxrf,s0reff = ns0peak(ns0reff)
  vs0refenef = ns0reff['e']
  vs0reff = ns0reff['s0']/1.e6

  ns0refe = ncread("ns0refe","e:s0:s1:s2:s3",fs0refe,silent=1,skiphead=2)
  iemaxre,emaxre,s0refe = ns0peak(ns0refe)
  vs0refenee = ns0refe['e']
  vs0refe = ns0refe['s0']/1.e6

  ns0refef = ncread("ns0refef","e:s0:s1:s2:s3",fs0refef,silent=1,skiphead=2)
  iemaxref,emaxref,s0refef = ns0peak(ns0refef)
  vs0refeneef = ns0refef['e']
  vs0refef = ns0refef['s0']/1.e6

  phwav = -9999.
  phdeg = -9999.

  try:
    F = open(fs0,"r")
    runs0 = F.readline().strip().split()[0].strip()
    F.close()
    iruns0 = int(runs0)
    phwav,phdeg = get_phase_error(runs0)
    if iruns0 < Irunmin: Irunmin = iruns0
    if iruns0 > Irunmax: Irunmax = iruns0
    ns0 = ncread("ns0","e:s0:s1:s2:s3",fs0,silent=1,skiphead=2)
    if len(ns0) != nener: print("\n*** Warning: Strange number of energies for ns0 ***")
  except:
    print("*** File",fs0,"not found ***")
    Quit("Schlecht")
    return iruns0,iruns0e,iruns0f,iruns0ef,iruns0ef
  #endtry

  try:
    F = open(fs0e,"r")
    runs0e = F.readline().strip().split()[0].strip()
    F.close()
    iruns0e = int(runs0e)
    if iruns0e < Irunmin: Irunmin = iruns0e
    if iruns0e > Irunmax: Irunmax = iruns0e
    ns0e = ncread("ns0e","e:s0:s1:s2:s3",fs0e,silent=1,skiphead=2)
    if len(ns0e) != nener: print("\n*** Warning: Strange number of energies for ns0e ***")
  except: pass

  try:
    F = open(fs0f,"r")
    runs0f = F.readline().strip().split()[0].strip()
    F.close()
    iruns0f = int(runs0f)
    if iruns0f < Irunmin: Irunmin = iruns0f
    if iruns0f > Irunmax: Irunmax = iruns0f
    ns0f = ncread("ns0f","e:s0:s1:s2:s3",fs0f,silent=1,skiphead=2)
    if len(ns0f) != nener: print("\n*** Warning: Strange number of energies for ns0f ***")
  except: pass

  try:
    F = open(fs0ef,"r")
    runs0ef = F.readline().strip().split()[0].strip()
    iruns0ef = int(runs0ef)
    F.close()
    if iruns0ef < Irunmin: Irunmin = iruns0ef
    if iruns0ef > Irunmax: Irunmax = iruns0ef
    ns0ef = ncread("ns0ef","e:s0:s1:s2:s3",fs0ef,silent=1,skiphead=2)
    if len(ns0ef) != nener: print("\n*** Warning: Strange number of energies for ns0ef ***")
  except: pass

  nplc(ns0,"e:s0/1.0e6")

  iemax,emax,s0max = ns0peak(ns0)
  if iemax < 0: s0max = -1.

  vs0ene = ns0['e']
  vs0 = ns0['s0']/1.e6
  vs0e = vs0 * 0.0
  vs0f = vs0 * 0.0
  vs0ef = vs0 * 0.0

  eharmf = 0.0
  eharme = 0.0
  eharmef = 0.0
  s0harmf  =  1.0
  s0harme  =  1.0
  s0harmef =  1.0

  legend("S0max, rref   " + g3(s0max) + "   " + g5(s0max/s0ref))

  cruns = str(Irunmin) + "-" + str(Irunmax)
  tit = "On-axis flux-density for " + key + " (Runs " + cruns + ")"
  txyz(tit,"E$_{ph}$ [eV]","N$_{ph}$/s/mm$^{2}$/100mA/0.1%BW")

  Fred = open("a/real_beam_" + key.strip() + "_" + cruns + ".dat","w")

  print("\nS0ref:   " + g5(s0ref))
  print("S0max:   " + g5(s0max) + "   " + g5(s0max/s0ref))
  Fred.write("S0ref: " + g5(s0ref) + "\n")
  Fred.write("S0max: " + g5(s0max) + "   " + g5(s0max/s0max)  + "   " + g5(s0max/s0ref) + "\n")

  if nexist("ns0f"):
    optnstat()
    nplcgs(ns0f,"e:s0/1.0e6")
    iemaxf,emaxf,s0maxf = ns0peak(ns0f)
    if iemax > 0:
      s0harmf = ns0f.s0[iemax]/1.e6
      eharmf = ns0f.e[iemax]
    vs0enef = ns0f['e']
    vs0f = ns0f['s0']/1.e6
    rdf = s0maxf/s0max
    rdhf = s0harmf/s0max
    line1 = "S0max_f,  Emax_f, rdf  : " + g3(s0maxf)  + BL + g5(emaxf)  \
    + BL + g3(rdf)
    line2 = "S0harm_f, Emax_p,  rdhf : " + g3(s0harmf) + BL + g5(eharmf) \
    + BL + g3(rdhf)
    legend(line1 + "\n" + line2 + "   " + g5(s0max/s0ref))
    print("\n" + line1 + "\n" + line2 + "   " + g5(s0max/s0ref))
    Fred.write("\n" + line1 + "\n" + line2 + "\n")
  #endif

  if nexist("ns0e"):
    optnstat()
    nplcbs(ns0e,"e:s0/1.0e6")
    iemaxe,emaxe,s0maxe = ns0peak(ns0e)
    if iemax > 0:
      s0harme = ns0e.s0[iemax]/1.e6
      eharme = ns0e.e[iemax]
    vs0enee = ns0e['e']
    vs0e = ns0e['s0']/1.e6
    rde = s0maxe/s0max
    rdhe = s0harme/s0max
    line1 = "S0max_e,  Emax_e, rde  : " + g3(s0maxe)  + BL + g5(emaxe)  \
    + BL + g3(rde)
    line2 = "S0harm_e, Emax_p,  rdhe : " + g3(s0harme) + BL + g5(eharme) \
    + BL + g3(rdhe)
    legend(line1 + "\n" + line2 + "   " + g5(s0max/s0ref))
    print("\n" + line1 + "\n" + line2 + "   " + g5(s0max/s0ref))
    Fred.write("\n" + line1 + "\n" + line2 + "\n")
  #endif

  if nexist("ns0ef"):
    optnstat()
    nplccs(ns0ef,"e:s0/1.0e6")
    iemaxef,emaxef,s0maxef = ns0peak(ns0ef)
    if iemax > 0:
      s0harmef = ns0ef.s0[iemax]/1.e6
      eharmef = ns0ef.e[iemax]
    vs0eneef = ns0ef['e']
    vs0ef = ns0ef['s0']/1.e6
    rdef = s0maxef/s0max
    rdhef = s0harmef/s0max
    line1 = "S0max_ef,  Emax_ef, rdef  : " + g3(s0maxef)  + BL + g5(emaxef)  \
    + BL + g3(rdef)
    line2 = "S0harm_ef, Emax_p,  rdhef : " + g3(s0harmef) + BL + g5(eharmef) \
    + BL + g3(rdhef)
    legend(line1 + "\n" + line2 + "   " + g5(s0max/s0ref))
    print("\n" + line1 + "\n" + line2 + "   " + g5(s0max/s0ref))
    Fred.write("\n" + line1 + "\n" + line2 + "\n")
  #endif

  Fred.close()

  pp("a/real_beam_" + key.strip() + "_" + cruns + ".pdf")
  legend()
  pp("a/real_beam_" + key.strip() + "_" + cruns + "_legend.pdf")

  nsig = nget("nsig")

  Fall = open(fall,"a")

  s0maxp =s0max
  s0refp =s0ref

  ksig = int(key.split("_")[1]) - 1
  phsig = nsig.phsig[ksig]

  #  Type beam  run  S0p  S0_pen S0H_pen S0H S0p/S0_pen S0p/S0H_pen S0p/S0H
  res = key + " pencil  " + runs0 + BL + g5(phsig) + BL + g5(phwav)  + BL + g5(phdeg) + BL \
  + g5(s0maxp) + BL + g5(s0maxp) + BL \
  + g5(s0refp) + BL + g5(s0refp) + BL \
  + g5(s0maxp/s0maxp) + BL + g5(s0maxp/s0refp) + BL + g5(s0maxp/s0refp)

  print(res)
  Fall.write(res+nl)

  #  Type beam  run  S0f  S0_pen S0H_en S0H S0f/S0_pen S0f/S0H_pen S0f/S0H
  res = key + " emit    " + runs0f + BL + g5(phsig) + BL + g5(phwav)  + BL + g5(phdeg) + BL \
  + g5(s0maxf) + BL + g5(s0maxp) + BL \
  + g5(s0refp) + BL + g5(s0reff) + BL \
  + g5(s0maxf/s0maxp) + BL + g5(s0maxf/s0refp) + BL + g5(s0maxf/s0reff)

  print(res)
  Fall.write(res+nl)

  #  Type beam  run  S0e  S0_pen S0H_pen S0H S0e/S0_pen S0e/S0H_pen S0e/S0H
  res = key + " espread " + runs0e + BL + g5(phsig) + BL + g5(phwav)  + BL + g5(phdeg) + BL \
  + g5(s0maxe) + BL + g5(s0maxp) + BL \
  + g5(s0refp) + BL + g5(s0refe) + BL \
  + g5(s0maxe/s0maxp) + BL + g5(s0maxe/s0refp) + BL + g5(s0maxe/s0refe)

  print(res)
  Fall.write(res+nl)

  #  Type beam  run  S0ef  S0_pen S0H_pen S0H S0ef/S0_pen S0ef/S0H_pen S0ef/S0H
  res = key + " em+esp  " + runs0ef + BL + g5(phsig) + BL + g5(phwav)  + BL + g5(phdeg) + BL \
  + g5(s0maxef) + BL + g5(s0maxp) + BL \
  + g5(s0refp) + BL + g5(s0refef) + BL \
  + g5(s0maxef/s0maxp) + BL + g5(s0maxef/s0refp) + BL + g5(s0maxef/s0refef)

  print(res)

  Fall.write(res+nl)

  Fall.close()

  optnstat()

  return vs0ene,vs0,vs0f,vs0e,vs0ef

#enddef

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

if nargs > 1:

    if args[1] == "wave":
        import waveplot as w
        from waveplot import *
    #endif

    if args[1] == "power_cpmu17":
        import waveplot as w
        from waveplot import *
        if nexist("n2000") == 0:
            args[1] = "tribs_CPMU17"
            print("*** Kein Spektrum berechnet, schalte auf ,,tribs_CPMU17'' um ***")
        #endif
    #endif

    if args[1] == "none":
      h2 = hbook2("Hid","htit",nx=1000,ny=1000)
      Quit("Ende")
    elif args[1] == "cpmu17_mess":

        var="ener:H:K:L:Epoch:ringc:kth11:flux:kth13:kth14:kth25:kth26:kth27:kth28:S2kt0:S2kt1:pzhei:pzpit:pzrol:gapSV:gapRB:tG800:T_Cr1:T_Cr2:scant:Sec1:Sec2"

        n29=ncread("n29",var,"emil_2021_3429.dat",skiphead=199)
        n30=ncread("n30",var,"emil_2021_3430.dat",skiphead=199)
        n32=ncread("n32",var,"emil_2021_3432.dat",skiphead=199)
        n33=ncread("n33",var,"emil_2021_3433.dat",skiphead=199)
        n38=ncread("n38",var,"emil_2021_3438.dat",skiphead=199)
        n63=ncread("n63",var,"emil_2021_3463.dat",skiphead=199)

        nplc(n29,"ener:flux")
        nplcs(n30,"ener:flux")
        nplcs(n32,"ener:flux")
        nplcs(n33,"ener:flux")
        nplcs(n38,"ener:flux")
#        nplcbs(n63,"ener:flux")

#        x=vcre(7,0.,10.)
#        y=1.3+sin(x)+sin(2*x)+0.3*cos(3*x)
#        #y=cos(x)
#        z=cos(x)+cos(2*x)+cos(3*x)+cos(4*x)
#        vwritexyz(x,y,z,"test.dat")
        pass

    elif args[1] == "spectra":

            nwt=ncread("nwt","e:s0:s1:s2:s3", \
            "wave_halbach_tab_pencil.dat" \
            ,skiphead=2)

            nsph=ncread("nsph","e:s0", \
            "~/spectra/spectra_halbach_pencil.dat" \
            ,skiphead=1)

            nspt=ncread("nspt","e:s0", \
            "~/spectra/spectra_halbach_tab_pencil.dat" \
            ,skiphead=1)

            optnstat()
            nplc(nwt,"e:s0/1.e6",legend="WAVE")
            nplcbs(nspt,"e:s0",legend="SPECTRA, tab.")
            nplcgs(nsph,"e:s0",legend="SPECTRA, mod.")
            legend()
            txyz("U32, 7th Harm.","E$_{ph}$ [eV]","N$_{ph}$/s/mm$^{2}$/100mA/0.1%BW")

            print("Reduction for tab:",nspt.s0.max()/(nwt.s0.max()/1.e6))
            print("Reduction for mod:",nsph.s0.max()/(nwt.s0.max()/1.e6))

            pp("wave_spectra_fluxden.pdf")

    elif args[1] == "test":
      n=ncread("n","x:y","paul.dat")
      npll(n,"x:y")
      nplcbs(n,"x+0.1:y")

    elif args[1] == "zone":

      plt.rcParams['axes.grid'] = True

      fig,axs = plt.subplots(3,1,sharex=True)
      plt.subplots_adjust(top=0.92, bottom=0.08, left=0.10, right=0.95,
                          hspace=0.1,wspace=0.35)
      axs[0].plot([1,2],[3,4])
      axs[1].plot([0,1],[2,3])
      axs[2].plot([0,1],[2,3])

      plt.show(block=False)

    elif args[1] == "booster":
      
      import waveplot as w
      from waveplot import *
      
      for n in range(len(Nhead)):
        snam = Nhead[n][1]
        exec(snam + ' = nget("' + snam + '")')
      #endfor
      
      kplot = 1
      
      if kplot == 2:
        
        n2 = ncread("n2","ebeam:flux","booster_flux_2eV.dat")
        npl(n2,"ebeam:flux")
        txyz("Flux for 2eV, pinhole 10mm x 10mm in 10m","Ebeam [GeV]","Flux")
        pp("booster_flux_2eV.pdf")
        
      if kplot == 1:
        
        zone(2,2)
        nby()
        nextzone()
        ndistpin()
        nextzone()
        #lolo()
        npl(n3601,"ener:fldn/1.e6")
        txyz("Central Flux-density","Eph")
        nextzone()
        npl(n3600,"ener:flux")
        txyz("Flux","Eph")
        pp("booster_"+str(Webea)+"GeV_overview.pdf")
        win()
        ndistpin('Fd',2)
        pp("booster_"+str(Webea)+"GeV_dist_2eV.pdf")
        nstat(n3600,"flux","iene==2")
        #f2 = open("booster_flux_2eV.dat","a")
        #f2.write(str(Webea)+"   "+pg4(w.Nmean)+"\n")
        #f2.close()
        wans()
      #endif kplot
      
    elif args[1] == "phase":
      
      import waveplot as w
      from waveplot import *
      
      for n in range(len(Nhead)):
        snam = Nhead[n][1]
        exec(snam + ' = nget("' + snam + '")')
      #endfor
      
      sy0 = "abs(y)<1.0e-10"
      siy = "iy==10"
      siziy = "iz==5 and iy==10"
      sc = str(clight1)
      
      nstat(n6001,"re_x*rb_x*"+sc,siziy)
      px = w.Nmean
      
      nstat(n6001,"sqrt(re_y**2+re_z**2)",siziy)
      reayz = w.Nmean
      nstat(n6001,"sqrt(rb_y**2+rb_z**2)*"+sc,siziy)
      rebyz = w.Nmean
      
      nstat(n6001,"sqrt(re_x**2+re_y**2+re_z**2)",siziy)
      rea = w.Nmean
      nstat(n6001,"sqrt(rb_x**2+rb_y**2+rb_z**2)*"+sc,siziy)
      reb = w.Nmean
      
      nstat(n6001,"re_x*rb_x*"+sc+"+re_z*rb_z*"+sc+"+re_y*rb_y*"+sc,siziy)
      p = w.Nmean
      #        nstat(n6001,"re_x*rb_x*"+sc,siziy)
      #        nstat(n6001,"re_z*rb_z*"+sc+"+re_y*rb_y*"+sc,siziy)

      print(p7g(reayz-rea))
      print(p7g(rebyz-reb))
      print(p7g(rea))
      print(p7g(reb))
      print(p5g(p/reayz/rebyz),p5g(acosd(p/reayz/rebyz)))
      print(p5g(p/rea/reb),p5g(acosd(p/rea/reb)))
      #print(p5g(px))
      #print(p5g(p-px))

      kplot = 5
      
      if kplot==5:

        zone(2,2)
        
        nplt(n6001,"z:re_z",siy,plopt='L')
        nnplt(n6001,"z:im_z",siy,plopt='L')
        
        nnplt(n6001,"z:3.0e8*rb_y",siy,plopt='L')
        nnplt(n6001,"z:3.0e8*ib_y",siy,plopt='L')
        
        nstat(n6000,"sqrt(re_x**2+re_y**2+re_z**2)",siziy)
        rea = w.Nmean
        nstat(n6000,"sqrt(rb_x**2+rb_y**2+rb_z**2)*"+sc,siziy)
        reb = w.Nmean
        
        nstat(n6000,"re_x*rb_x*"+sc+"+re_z*rb_z*"+sc+"+re_y*rb_y*"+sc,siziy)
        p = w.Nmean
        
        print(p7g(rea))
        print(p7g(reb))
        print(p5g(p/rea/reb),p5g(acosd(p/rea/reb)))

      elif kplot==4:

        zone(2,2)
        
        #nplt( n6001,"z:re_x",siziy,plopt='L')
        #nplt(n6001,"z:re_y",siziy,plopt='L')
        nplt(n6001,"z:re_z",siziy,plopt='L')
        #nnplt(n6001,"z:im_x",siziy,plopt='L')
        #nnplt(n6001,"z:im_y",siziy,plopt='L')
        nnplt(n6001,"z:im_z",siziy,plopt='L')
        
        #nnplt(n6001,"z:3.0e8*rb_x",siziy,plopt='L')
        nnplt(n6001,"z:3.0e8*rb_y",siziy,plopt='L')
        #nnplt(n6001,"z:3.0e8*rb_z",siziy,plopt='L')
        #nnplt(n6001,"z:3.0e8*ib_x",siziy,plopt='L')
        nnplt(n6001,"z:3.0e8*ib_y",siziy,plopt='L')
        #nnplt(n6001,"z:3.0e8*ib_z",siziy,plopt='L')
        
        nscan(n6001,"re_x:re_y:re_z:re_x*rb_x*"+sc+"+re_z*rb_z*"+sc+"+re_y*rb_y*"+sc+"",siziy)
        nscan(n6001,"rb_x*"+sc+":rb_y*"+sc+":rb_z*"+sc+":re_x*rb_x*"+sc+"+re_z*rb_z*"+sc+"+re_y*rb_y*"+sc+"",siziy)
        print("\n")
        nscan(n6001,"im_x:im_y:im_z:re_x*rb_x*"+sc+"+re_z*rb_z*"+sc+"+re_y*rb_y*"+sc+"",siziy)
        nscan(n6001,"ib_x*"+sc+":ib_y*"+sc+":ib_z*"+sc+":re_x*rb_x*"+sc+"+re_z*rb_z*"+sc+"+re_y*rb_y*"+sc+"",siziy)
        
      elif kplot==3:
        
        zone(1,2)

        setxstat(0.2)
        npl(n6001,"z:(rb_y**2+ib_y**2)*9.e16",sy0)
        setxstat(0.8)
        nplmgs(n6001,"z:re_z**2+im_z**2",sy0)

        setxstat(0.2)
        nnpl(n6000,"z:(rb_y**2+ib_y**2)*9.e16",sy0)
        #nnpl(n6000,"z:(rb_z**2+ib_z**2)",sy0)
        setxstat(0.8)
        nplmgs(n6000,"z:re_z**2+im_z**2",sy0)
        
        ninfo(n6000)
        nstat(n6000,"(rb_y**2+ib_y**2)*9.e16",sy0)
        nstat(n6000,"re_z**2+im_z**2",sy0)

      elif kplot==2:
        
        zone(2,4)
        
        #nplt( n6001,"z:re_x",sy0,plopt='L')
        #nplt(n6001,"z:re_y",sy0,plopt='L')
        #nnplt(n6001,"z:im_y",sy0,plopt='L')
        
        #nnpll(n6001,"z:re_x*rb_x*"+sc+"+re_z*rb_z*"+sc+"+re_y*rb_y*"+sc+"",sy0,plopt='L')
        #nnplt(n6001,"z:rb_x*"+sc+"",sy0,plopt='L')
        #nnplt(n6001,"z:rb_y*"+sc+"",sy0,plopt='L')
        #nnplt(n6001,"z:rb_z*"+sc+"",sy0,plopt='L')
        #nnplt(n6001,"z:rb_y*"+sc+"/re_z+1",sy0,plopt='L')
        #Quit()
        
        #nstat(n6001,"re_z/(-rb_y)/3.0e8",sy0)
        #nstat(n6001,"re_y/(rb_z)/3.0e8",sy0)
        nstat(n6001,"re_x*rb_x*"+sc+"+re_z*rb_z*"+sc+"+re_y*rb_y*"+sc+"",sy0)
        #Quit()

        #nnplt(n6000,"z:re_x",sy0,plopt='L')
        nplt( n6000,"z:re_y",sy0,plopt='L')
        nnplt(n6000,"z:im_y",sy0,plopt='L')
        nnplt(n6000,"z:re_z",sy0,plopt='L')
        nnplt(n6000,"z:im_z",sy0,plopt='L')
        
        #nnpll(n6000,"z:re_x*rb_x*"+sc+"+re_z*rb_z*"+sc+"+re_y*rb_y*"+sc+"",sy0,plopt='L')
        #nnplt(n6000,"z:rb_x*"+sc+"",sy0,plopt='L')
        nnplt(n6000,"z:rb_y*"+sc+"",sy0,plopt='L')
        nnplt(n6000,"z:ib_y*"+sc+"",sy0,plopt='L')
        nnplt(n6000,"z:rb_z*"+sc+"",sy0,plopt='L')
        nnplt(n6000,"z:ib_z*"+sc+"",sy0,plopt='L')
        
        #nstat(n6000,"re_z/(-rb_y)/3.0e8",sy0)
        #nstat(n6000,"re_y/(rb_z)/3.0e8",sy0)
        nstat(n6000,"re_x*rb_x*"+sc+"+re_z*rb_z*"+sc+"+re_y*rb_y*"+sc+"")
        
      elif kplot==1:
        ninfo(n6001)
        zone(2,2)
        ndistpinh()
        nextzone()
        ndistpinh('AzR')
        nextzone()
        ndistphaseh()
        nextzone()
        ndistphaseh('AzR')
        #Quit()
      #endif
          
    elif args[1] == "default" or args[1]=="d":
      
      import waveplot as w
      from waveplot import *
      
      for n in range(len(Nhead)):
        snam = Nhead[n][1]
        exec(snam + ' = nget("' + snam + '")')
      #endfor
      
      kdef = 1
      
      if (kdef == 1 and fexist('urad_phase.wig')):
        nwig = ncread("nwig","iz:iy:itz:ity:x:y:z:ty:tz:iegam:egam:wzzr:wzzi:wzyr:wzyi:wyzr:wyzi:wyyr:wyyi","urad_phase.wig")      
        ninfo(nwig)
        #zone(3,1)
        #npl(nwig,"iz:itz:wzzr:wzzr","y==0 and ty==0")
        #npl(nwig,"z:y:wzzr:wzzr","tz==0 and ty==0")
#        npl(nwig,"z:tz:wzzr:wzzr","y==0 and ty==0")
        #nnpl(nwig,"tz:ty:wzzr:wzzr","y==0 and z==0")
      #endif kdef
      
      iwave = 0
      
      if iwave == 12:
        zone(2,2)
        ndistpin('azr')
        nextzone()
        ndistpin('azi')
        nextzone()
        ndistpin('ayr')
        nextzone()
        ndistpin('ayi')
        #nextzone()
      #endif
      
      if iwave == 11:
        #zone(1,2)
        set_y_stat(0.7)
        npll("n3700","z:re_z","abs(y)<1.0e-10")
        owf('../stage/mhb')
#        sethistcolor('b')
        set_x_stat(0.7)
        npllgs("n3700","z:re_z","abs(y)<1.0e-10")
      #endif
      
      if iwave == 10:
        
        nfld = ncread("nfld","x:y:z:iegam:egam:s0:s1:s2:s3:p:exr:exi:eyr:eyi:ezr:ezi:bxr:bxi:byr:byi:bzr:bzi:nx:ny:nz","~/urp/urad_phase.fld")
        
        s0mx=nmax(n3700,"spec","y==0")
        us0mx=nmax(nfld,"s0","y==0")
        emx=nmax(n3700,"re_z","y==0")
        uemx=nmax(nfld,"ezr","y==0")
        
        zone(4,1)
        
        npll("n3700","z*1000.:spec/"+str(s0mx),"y==0")
        nspec=ncread("nspec","z:spec",getwavedump())
        set_y_stat(0.3)
        npllls("nfld","z:s0/"+str(us0mx),"y==0")
        ns0=ncread("ns0","z:s0",getwavedump())
        
        nextzone()        
        optnstat()
        
        npll("n3700","z*1000.:re_z/"+str(emx),"y==0")
        nre_z=ncread("nre_z","z:re_z",getwavedump())
        npllls("nfld","z:ezr/"+str(uemx),"y==0")
        nezr=ncread("nezr","z:ezr",getwavedump())
        
        set_y_stat(0.8)
        ezr=nezr.ezr
        re_z=nre_z.re_z
        s0=ns0.s0
        spec=nspec.spec
        ds0=ns0.s0-nspec.spec
        rs0=ns0.s0/nspec.spec
        z=ns0.z
        
        nextzone()        
        optstat()        
        vplxy(z,rs0,'line')
        
        nextzone()        
        vplxy(z,ds0,'line')
        
      #endif
      
      if iwave == 9:        
        zone(2,2)
        ndistpin()
        nextzone()
        ndistpinh('ezr')
        nextzone()
        ndistphase()
        nextzone()
        ndistphaseh('ezr')
      #endif

      if iwave == 8 and os.path.exists("WAVE.mhb"):        
        mhb_to_pylist("/home/scheer/wav/wig/WAVE.mhb")
        optnstat()        
        mhb_cd_up()
        iy = 1
        npll("n3700","z*1000.:re_z","iy==" + str(iy) + " and iene==1")
#        nplmrs("n3700","z*1000.:re_z","iy==" + str(iy) + " and iene==1")
        mhb_cd_down()
        npllbs("n3700","z*1000.:re_z","iy==" + str(iy) + " and iene==1")
#        nplmbs("n3700","z*1000.:re_z","iy==" + str(iy) + " and iene==1")
      #endif
      
      if iwave == 7 and os.path.exists("WAVE.mhb"):        
        mhb_to_pylist("/home/scheer/wav/wig/WAVE.mhb")
        optnstat()        
        zone(1,2)
        mhb_cd_up()
        iy = 1
        npll("n3700","z*1000.:re_z","iy==" + str(iy) + " and iene==1")
        nplmrs("n3700","z*1000.:re_z","iy==" + str(iy) + " and iene==1")
        mhb_cd_down()
        nplmbs("n3700","z*1000.:re_z","iy==" + str(iy) + " and iene==1")
        nextzone()
        mhb_cd_up()
        iy = int(n3700.iy.max()/2+1)
        npll("n3700","z*1000.:re_z","iy==" + str(iy) + " and iene==1")
        nplmrs("n3700","z*1000.:re_z","iy==" + str(iy) + " and iene==1")
        mhb_cd_down()
        nplmbs("n3700","z*1000.:re_z","iy==" + str(iy) + " and iene==1")
      #endif
      
      if iwave == 6 and os.path.exists("WAVE.mhb"):
        
        mhb_to_pylist("/home/scheer/wav/wig/WAVE.mhb")
        mhb_cd_up()
        
        optnstat()
        zone(1,3)
        
        npll("n3700","z*1000.:re_z","iy==1 and iene==1")
        mhb_cd_down()
        npllbs("n3700","z*1000.:re_z","iy==1 and iene==1")
        
        nextzone()
        mhb_cd_up()
        npll("n3700","z*1000.:im_z","iy==1 and iene==1")
        mhb_cd_down()
        npllbs("n3700","z*1000.:im_z","iy==1 and iene==1")
        
        nextzone()
        optstat()
        mhb_cd_up()
        npll("n3700","z*1000.:re_z**2+im_z**2","iy==1 and iene==1")
        mhb_cd_down()
        set_x_stat(0.8)
        npllbs("n3700","z*1000.:re_z**2+im_z**2","iy==1 and iene==1")
        
      #endif
      
      if iwave == 5 and os.path.exists("WAVE.mhb"):
        optnstat()
        zone(1,4)
        npll(n3700,"z*1000.:re_z","iy==1 and iene==1")
        # nextzone()
        npllgs(n3700,"z*1000.:im_z","iy==1 and iene==1")
        nextzone()
        optstat()
        npll(n3700,"z*1000.:re_z**2+im_z**2","iy==1 and iene==1")
        optnstat()
        nextzone()
        npll(n6000,"z*1000.:re_z","iy==1 and ie==1")
        # nextzone()
        npllgs(n6000,"z*1000.:im_z","iy==1 and ie==1")
        nextzone()
        optstat()
        npll(n6000,"z*1000.:re_z**2+im_z**2","iy==1 and ie==1")
        optnstat()
      #endif
      
      if iwave == 4 and os.path.exists("WAVE.mhb"):
        
        optnstat()
        
        ndump("n3700","z*1000.:re_z:im_z:spec","iy==1 and iene==1","w2.dat")
        ndump("n6001","z*1000.:re_z:im_z:spec","iy==1 and ie==1","w2ph.dat")
        
        mhb_to_pylist("/home/scheer/wav/wig/WAVE.mhb")
        
        ndump("n3700","z*1000.:re_z:im_z:spec","iy==1 and iene==1","wig.dat")
        ndump("n6001","z*1000.:re_z:im_z:spec","iy==1 and ie==1","wigph.dat")
        
        nw2=ncread("nw2","z:re:im:spec","w2.dat")
        nwig=ncread("nwig","z:re:im:spec","wig.dat")
        nw2ph=ncread("nw2ph","z:re:im:spec","w2.dat")
        nwigph=ncread("nwigph","z:re:im:spec","wig.dat")
        
        zone(1,2)
        
        npll(nw2,"z:re**2+im**2")
        npllls(nw2ph,"z:re**2+im**2")
        #      npllbs(nwig,"z:re**2+im**2")
        #      npllcs(nwigph,"z:re**2+im**2")
        
        nextzone()
        npll(nw2,"z:re")
        npllls(nw2ph,"z:re")
        #      npllbs(nwig,"z:re")
        #      npllcs(nwig,"z:re")
      #endif
      
      if iwave == 3 and os.path.exists("WAVE.mhb"):
          
        import waveplot as w
        from waveplot import *

        set_console_title("wavesDefault")
        for n in range(len(Nhead)):
          snam = Nhead[n][1]
          exec(snam + ' = nget("' + snam + '")')
        #endfor
        
        for hh in H1head:
          snam = hh[0]
          exec(snam + ' = hget("' + snam + '")')
        #endfor
        
        for hh in H2head:
          snam = hh[0]
          exec(snam + ' = hget("' + snam + '")')
        #endfor

        #reakpoint()
        nbeam('f')
        #hflux('f',"S")
        
      #endif
      
      if iwave == 2 and os.path.exists("WAVE.mhb"):

        import waveplot as w
        from waveplot import *

        set_console_title("wavesDefault")
        for n in range(len(Nhead)):
          snam = Nhead[n][1]
          exec(snam + ' = nget("' + snam + '")')
        #endfor
        
        for hh in H1head:
          snam = hh[0]
          exec(snam + ' = hget("' + snam + '")')
        #endfor
        
        for hh in H2head:
          snam = hh[0]
          exec(snam + ' = hget("' + snam + '")')
        #endfor

        optstat()
        #hcfluxden()
        #npl(n30,"ener:spec/1.e6","","","Sprof")
        npl(n30,"ener:spec/1.e6","nel==1")

        get_console(w.Console)
        
      #endif

      if iwave == 1 and os.path.exists("WAVE.mhb"):

        import waveplot as w
        from waveplot import *

        set_console_title("wavesDefault")
        for n in range(len(Nhead)):
          snam = Nhead[n][1]
          exec(snam + ' = nget("' + snam + '")')
        #endfor
        
        for hh in H1head:
          snam = hh[0]
          exec(snam + ' = hget("' + snam + '")')
        #endfor
        
        for hh in H2head:
          snam = hh[0]
          exec(snam + ' = hget("' + snam + '")')
        #endfor
          
        optstat()
        
        if Wibun:
          zone(1,3)
          nx = hcopy1d("h148000","h148")
          nproj1(n30,"ener","spec/1.e6","nel==1",idh="h148")
          h148 = hget("h148")
          hplot1d("h148")
        else: 
          zone(1,2)
          hplot1("h148000")
        #endif
        
        nextzone()
        hplot1("h48000")

        fdxopt = 0; fdyopt = 0; fdxrms = 0
        fdxoptr = 109.98967856671656
        fdyoptr = 101503309178039.47
        fdxrmsr = 1.4193676101680432

        fxopt = 0; fdopt = 0; fdrms = 0
        
        fxoptr = 107.39587225780122
        fyoptr = 339804478603828.06
        fxrmsr = 1.914615402497087
        
        if hexist("h148000"):
          if Wibun:
            fdsumn,fdsumy,fdxmean,fdxrms,fdxopt,fdyopt = hstat1d("h148")
          else:
            fdsumn,fdsumy,fdxmean,fdxrms,fdxopt,fdyopt = hstat1d("h148000")
          #endif          
          print("Fd:",fdxopt,fdyopt,fdxrms)
        #endif
        
        if hexist("h48000"):
          fsumn,fsumy,fxmean,fxrms,fxopt,fyopt = hstat1d("h48000")
          print("F: ",fxopt,fyopt,fxrms)
        #endif
        
        if hexist("h148000"):
          try:
            if Wneib <= 1: print("Dev. Fd:",g5(1-fdxopt/fdxoptr),g5(1-fdyopt/fdyoptr),g5(1-fdxrms/fdxrmsr))
            else: print("Dev. Fd:",g5(1-fdxopt/fdxoptr),g5(fdyopt/fdyoptr),g5(1-fdxrms/fdxrmsr))
          except: pass
        #endif
        
        if hexist("h48000"):
          try: 
            if Wneib < 2: print("Dev. F:",g5(1-fxopt/fxoptr),g5(1-fyopt/fyoptr),g5(1-fxrms/fxrmsr))
            else: print("Dev. F:",g5(1-fxopt/fxoptr),g5(fyopt/fyoptr),g5(1-fxrms/fxrmsr))
          except: pass
        #endif
        
        if Wibun:
          nextzone()
          nbeam('f')
        #endif
        
        #Quit()
        get_console(w.Console)
        
      #endif

      idefault = -1
      
      if idefault == 5:
        nzm = ncread("nzm","e:fd","a/Wave_dipole_BESSY-III_z=-840.0mm.dat")
        nz0 = ncread("nz0","e:fd","a/Wave_dipole_BESSY-III_z=0mm.dat")
        nzp = ncread("nzp","e:fd","a/Wave_dipole_BESSY-III_z=1306.9mm.dat")

        winl()
        sel = "fd > 1000."
        npll(nz0,"e:fd",sel)
        legend("z=0.0m")
        npllgs(nzm,"e:fd",sel)
        legend("z=-0.84m")
        npllbs(nzp,"e:fd",sel)
        legend("z=-1.31m")
        legend()
        txyz("Flux-density of BESSY III Dipole",xTit,yTit)
        nplmrs(nz0,"e:fd",sel)
        nplmbs(nzp,"e:fd",sel)
        nplmgs(nzm,"e:fd",sel)
        pp("Wave_dipole_BESSY-III.pdf")
      #endif
      
      if idefault == 4:

        optstat()
        ner = ncread("ner","be:phe:fd","wave_serie_phase-error_20220505.out")
        fdmax=ner.fd.max()
        svar = "phe:fd/" + str(fdmax)
        npl(ner,svar)

        nerall = ncread("nerall","be:phe:eopt:fdopt:ene:fd","wave_serie_phase-error_20220505_all.out")
        ninfo(nerall)

        optnstat()

        npll(nerall,"ene+1.7956976237883282:fd/8.25556e+10","abs(phe-4.995800)<1.e-6")
        npllgs(nerall,"ene+2500-2498.2541630992987:fd/8.25556e+10","be==0.004 and abs(phe-4.999086)<1.e-6")
        npllcs(nerall,"ene+2500-2500.9044799726225:fd/8.25556e+10","be==0.004 and abs(phe-4.996088)<1.e-6")
        pp("phase-error_5deg_1000m_norm_shifted.pdf")

      elif idefault == 3:

        optstat()
        ner = ncread("ner","be:phe:fd","wave_serie_phase-error_20220504.out")
        fdmax=ner.fd.max()
        svar = "phe:fd/" + str(fdmax)
        npl(ner,svar)

        nerall = ncread("nerall","be:phe:eopt:fdopt:ene:fd","wave_serie_phase-error_all_20220504.out")
        ninfo(nerall)

        npll(nerall,"ene-0.39574412052434127:fd/8.34694e+14","be == 0.004 and abs(phe-5)<0.1 and phe<4.925")
        npllgs(nerall,"ene+3.1456757159221524:fd/8.34694e+14","be == 0.004 and abs(phe-5)<0.1 and abs(phe-4.95)<0.001")
        npllcs(nerall,"ene+2.2143208185048024:fd/8.34694e+14","be == 0.004 and abs(phe-5)<0.1 and phe>5.075")
        npllbs(nerall,"ene+3.087656223977774:fd/8.34694e+14","be == 0.004 and abs(phe-5)<0.1 and abs(phe-5)<0.005")
        pp("phase-error_5deg_be_0.004_norm_shifted.pdf")

      elif idefault == 2:

        n6507=ncread("n6507","e:s0:s1:s2:s3","wave_stokes_selected.dat.6507")
        n6508=ncread("n6508","e:s0:s1:s2:s3","wave_stokes_selected.dat.6508")
        n6509=ncread("n6509","e:s0:s1:s2:s3","wave_stokes_selected.dat.6509")
        nlist()

        optnstat()
        npl(n6507,"e:s0/1.e6")
        npllgs(n6508,"e:s0/1.e6")
        npllcs(n6509,"e:s0/1.e6")

      #endif idefault

    elif args[1] == "ipac":
      
      import waveplot as w
      from waveplot import *
      
      zone(2,2)
      optnstat()
      optrun(False)
      w.Krun=False
      w.Kdate=False
      
      ntrack()
      nextzone()
      hflux()
      nextzone()
      ndistpin()
      nextzone()
      ndistpow()
      
      pp("ipac2023_wave.pdf")
      
    elif args[1] == "urs":
      
      optnstat()
      
#      no=ncread("no","i:x:z:zp:s:t:bx:bz:r:rnx:rny:rnz:nnbx:nnby:nnbz:\
#phi:om:azr:azi","~/wav/w2/urad_spline.out")
#      np=ncread('np',"t:ifr:azr:azi","paul.dat")
      no=ncread("no","ix:ifr:r1:r2:r3:phi:om:t:x:z:azr:azi","urad_spline.out")
      nf=ncread("nf","ifr:e:ayr:ayi:azr:azi:s","~/wav/w2/urad_spline.amp")
      
#      nnnb=ncread("nnnb","t2:x2:z2:rnx:rny:rnz:bx:by:bz:rnnbx:rnnby:rnnbz:\
#t0ph:cdumi:rargx:rargy:rargz:dum3r:dum3i:phase:om:azr:azi:phi:phir:r","~/wav/w1/urad_nnb.dat")
      
      nlist()
      ninfo(no)
      
      kplot = 0
      
      if kplot == -1:
        print("")
        
      elif kplot == 3:
        
        optstat()
        npll(nf,"e:azr**2+azi**2")
        
      elif kplot == 2:
        
        npll(nnnb,"x2:z2")
        nplmgs(no,"x:z")
        
      elif kplot == 1:
        zone(2,1)
        npll(no,"t:nnbz*sin(om*phi)")
        npl(no,"t:nnbz*cos(om*phi)","","","CS")
        npllgs(no,"t:nnbz*sin(om*phi)")
        npl(no,"t:nnbz*sin(om*phi)","","","CS",color='g')
        nex()
        optstat()
        npll(nf,"e:s")
        ninfo(nf)
      #endif
      
    elif args[1] == 'wiga':
      
      namp=ncread("namp","iz:iy:jz:jy:nm:np:zw:yw:zp:yp:er:ei","urad_phase.amp")
      ninfo(namp)
                    
      if (fexist('urad_phase.fdp')):        
        nfdp = ncread("nfdp","x:y:z:iegam:egam:s0:s1:s2:s3:exr:exi:eyr:eyi:ezr:ezi:bxr:bxi:byr:byi:bzr:bzi:nx:ny:nz","urad_phase.fdp")        
      #endif
      nlist()
      
    elif args[1] == "nbun":
        
        if (fexist('urad_phase.bun')):
            nbun = ncread("nbun","jbun:isub:ibu:bunchx:rxi1:ryi1:rzi1:ypi1:zpi1:rxin:ryin:rzin:ypin:zpin:eel:deel:x:y:z:iegam:egam:spec:s0:s1:s2:s3:p:fb28:dt:exr:exi:eyr:eyr:ezr:ezi:bxr:bxi:byr:byi:bzr:bzi","urad_phase.bun")
            #ninfo(nbun)
            zone(3,1)
            npl(nbun,"rzin:zpin")
            nex()
            npl(nbun,"rzin")
            nex()
            npl(nbun,"z:s0",plopt='prof')
            [sumn,sumy,xmean,sigz,xopt,yopt] = hstat("HnPlot")
            ystat(0.1)
            npl(nbun,"y:s0",plopt='Sprof',color='g')
            hy = hget("HnPlot")
            [sumn,sumy,xmean,sigy,xopt,yopt] = hstat("HnPlot")
            print(g3(sigz),g3(sigy),g3(sigz/sigy))
        else:
            print("Nbun nicht vorhanden!")
        #endif
        
    elif args[1] == "test_tanaka" or args[1] == "tt":
        
        ntau = ncread("ntau","x:y:z:t:betx:bety:betz:tau","test_tanaka.tau")
        nfld = ncread("nfld","epho:z:y:exr:eyr:ezr:exi:eyi:ezi","test_tanaka.fld")
        
        npl(nfld,"epho:(ezr**2+ezi**2)")
        
    elif args[1] == "urad" or args[1] == "pinh":
      
      if (fexist('urad_phase.fld')):
        nfld = ncread("nfld","x:y:z:iegam:egam:s0:s1:s2:s3:p:exr:exi:eyr:eyi:ezr:ezi:bxr:bxi:byr:byi:bzr:bzi:nx:ny:nz","urad_phase.fld")
        ninfo(nfld)
      else:
        Quit("Keine Daten!")
      #endif

      #reakpoint()
      ymina = sqrt(min(nfld.y**2))
      sel = 'y == ' + str(ymina)
      npll(nfld,"z:s0",sel)
      txyz("Flux density, "+sel,"z[mm]")
      
    elif args[1] == "urad" or args[1] == "ezr":
      
      if (fexist('urad_phase.fld')):
        nfld = ncread("nfld","x:y:z:iegam:egam:s0:s1:s2:s3:p:exr:exi:eyr:eyi:ezr:ezi:bxr:bxi:byr:byi:bzr:bzi:nx:ny:nz","urad_phase.fld")
        ninfo(nfld)
      else:
        Quit("Keine Daten!")
      #endif

      #reakpoint()
      ymina = sqrt(min(nfld.y**2))
      sel = 'y == ' + str(ymina)
      npll(nfld,"z:ezr",sel)
      txyz("EzR, "+sel,"z[mm]")
      
    elif args[1] == "urad" or args[1] == "u" or args[1] == "ud":
      
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

      #reakpoint()
      
      if (fexist('urad_phase.fld')):
        nfld = ncread("nfld","x:y:z:iegam:egam:s0:s1:s2:s3:p:exr:exi:eyr:eyi:ezr:ezi:bxr:bxi:byr:byi:bzr:bzi:nx:ny:nz","urad_phase.fld")
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
      
      #reakpoint()
      
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
      
      #nref0=ncread("nref0","z:s0","ref0.dat")
      #nref=ncread("nref","z:s0","ref_ff.dat")
          
      nlist()
      #wait(2)
      #reakpoint()
      
      #s0max = nfld.s0.max()
      #s3max = nfld.s3.max()
      #print(pg5(s0max))
      #if not isnan(s0max): kplot = 1
      
      #reakpoint()
      
      if args[1] == "ud":
        
       s0max = nfld.s0.max()
       emax = nfld.query("s0=="+str(s0max)).egam.max()
       selgam = "abs(egam-" + str(emax) + ")<1.0e-10"
       
       beth=9.8
       emith=7.685e-9
       betv=3.29
       emitv=0.01537e-9

       L = 0.051*78
       lam=evnm(emax)/1.e9
       
       sigrkim=sqrt(lam*L)/4/pi        
       sige=sqrt(beth*emith)
       sigep=sqrt(emith/beth)
       sigkim = sqrt(sigrkim**2+sige**2)
       sigrpkim=sqrt(lam/L)
       sigrpwalker=sqrt(lam/L/2)
       sigpkim = sqrt(sigrpkim**2+sigep**2)
       sigpwalker = sqrt(sigrpwalker**2+sigep**2)
       
       #reakpoint()
       kplot = 2
       
       if kplot == 5:
         zone(2,2)
         npl(nebm,"z:y","","s0",plopt='boxes')
         nex()
         npl(nebm,"z:y","","ezr",plopt='boxes')
         nex()
         npl(nebm,"nz*1000")
         nex()
         npl(nebm,"ny*1000")
         
       if kplot == 4:
         
         optdate()
         zone(2,1)
         nprof(nfdp,"z:s0")
         hz = hget('HnPlot')
         hplot(hz,'M')
         m.GtitFontSize='medium'
         gtit('X=' + str(pinxprop) + ', PinX=' + str(pinx) + ', PinW=' + str(pinw) + \
         ', PinH=' + str(pinh) + ', nh=' + str(npinz) + ', nv=' + str(npiny) + \
         '\n Modepin=' + str(modepin) + ', Ifold=' + str(ifold) + ', Nelec=' + str(nelec)
         + ', Ihbunch_' + str(ihbunch))
         txyz('z:S0, projet.')
         nex()
         nprof(nfdp,"y:s0")
         hy = hget('HnPlot')
         hplot(hy,'M')
         txyz('y:S0, projet.')
         pp("urad_phase_modepin_" + str(modepin) + '_ifold_' + str(ifold) + \
         '_ihbunch_' + str(ihbunch) + '.pdf')
         
       if kplot == 3:
         zone(2,1)
         npll(nfld,"z:ezr","y==0")
#         ystat(0.5)
#         npllgs(nfld,"z:ezi","y==0")
         nex()
         ystat(0.8)
         npll(nfdp,"z:ezr","y==0")
#         ystat(0.5)
#         npllgs(nfdp,"z:ezi","y==0")
         
       if kplot == 1:
         
         if ifold == 0: 
           zone(2,1)
           npll(nfld,"z:ezr/" + str(nfld.ezr.max()),"y==0")
           ystat(0.5)
           npllgs(nfld,"z:ezi/" + str(nfld.ezr.max()),"y==0")
           nex()
           ystat(0.8)
           npll(nfld,"y:ezr/" + str(nfld.ezr.max()),"z==0")
           ystat(0.5)
           npllgs(nfld,"y:ezi/" + str(nfld.ezr.max()),"z==0")
           #nf=ncread("nf","z:ezr:ezi","/hs/urp/ez_fold.dat")
           #optnstat()
           #npllbs(nf,"z:ezr/" + str(nf.ezr.max()))
           #npllcs(nf,"z:ezi/" + str(nf.ezr.max()))
           ndump(nfld,"z:ezr:ezi","y==0","ez.dat")
         else:
           zone(2,1)
           npll(nfdf,"z:ezr","y==0")
           ystat(0.5)
           npllgs(nfdf,"z:ezi","y==0")
           nex()
           ystat(0.8)
           npll(nfdf,"y:ezr","z==0")
           ystat(0.5)
           npllgs(nfdf,"y:ezi","z==0")
           ndump(nfdf,"z:ezr:ezi","y==0","ez_fold.dat")
         #endif
         pass
       
       if kplot == 2 and ifold == 1:
        
        m.GtitFontSize='small'

        if ifieldprop == 0:
          
          zone(2,2)
          
          hbook2('Hpin','distribution in pinhole',npinz,zmin,zmax,npiny,ymin,ymax)
          nproj2(nfld,"z:y","s0",selgam,idh='Hpin',ioverwrite=0)
          hplave('Hpin','boxes')
          txyz("S0 in Pinhole","z","y","S0")
          
          nex()
          hbook2('HpinF','folded distribution in pinhole',npinz,zmin,zmax,npiny,ymin,ymax)
          nproj2(nfdf,"z:y","s0",selgam,idh='HpinF',ioverwrite=1)
          hplave('HpinF','boxes')
          txyz("S0 in Pinhole","z","y","S0_folded")
          
          nex()
          hbook1('HpinzS0F','z:S0_folded',npinz,zmin,zmax)
          nproj1(nfdf,"z","s0",selgam,idh='HpinzS0F',ioverwrite=0)
          hpl('HpinzS0F','L')
          txyz("S0, proj., folded","z","S0")
          rmsS0Fz = hrms(HpinzS0PF,isilent=0)
          
          nex()
          hbook1('HpinyS0F','y:S0',npiny,ymin,ymax)
          nproj1(nfdf,"y","s0",selgam,idh='HpinyS0F',ioverwrite=0)
          hpl('HpinyS0F','L')
          txyz("S0, proj., folded","y","S0")
          
        else:
          
          zone(2,2)
          
          hbook2('HpinF','folded distribution in pinhole',npinz,zmin,zmax,npiny,ymin,ymax)
          nproj2(nfdf,"z:y","s0",selgam,idh='HpinF',ioverwrite=1)
          hplave('HpinF','boxes')
          txyz("S0_folded in Pinhole","z","y","",tity=1.0)
          
          nex()
          hbook2('HpinS0propF','Prop. S0, folded',npinzprop,zminph,zmaxph,npinyprop,yminph,ymaxph)
          nproj2(nfdpf,"z:y","s0",selgam,idh='HpinS0propF',ioverwrite=0)
          hplave('HpinS0propF','boxes')
          txyz("S0_folded in Source","z","y","",tity=1.0)
          
          nex()
          hbook1('HpinzS0PF','z:S0_folded',npinzprop,zminph,zmaxph)
          nproj1(nfdpf,"z","s0",selgam,idh='HpinzS0PF',ioverwrite=0)
          hpl('HpinzS0PF','L')
          txyz("S0, prop, proj., folded","z","S0")
          rmsS0Fz = hrms(HpinzS0PF,isilent=0)
          print("RMS S0,z:",rmsS0Fz)
          
          nex()
          hbook1('HpinyS0PF','y:S0_folded',npinyprop,yminph,ymaxph)
          nproj1(nfdpf,"y","s0",selgam,idh='HpinyS0PF',ioverwrite=0)
          hpl('HpinyS0PF','L')
          txyz("S0, prop, proj., folded","y","S0")
          rmsS0Fy = hrms(HpinyS0PF,isilent=0)
          print("RMS S0,y:",rmsS0Fy)
          
        #endif ifieldprop == 0:
        
        gtit('X=' + str(pinxprop) + ', PinX=' + str(pinx) + ', PinW=' + str(pinw) + \
        ', PinH=' + str(pinh) + ', nh=' + str(npinz) + ', nv=' + str(npiny) + \
        '\n Modepin=' + str(modepin) + ', Ifold=' + str(ifold) + ', Nelec=' + str(nelec) \
        + ', Ihbunch=' + str(ihbunch))
        
        fbase = ("urad_phase_modepin_" + str(modepin) + '_ifold_' + str(ifold) + \
        '_ihbunch_' + str(ihbunch))
        
        pp(fbase + '.pdf')
        
        os.system('cp fld ' + fbase + '.fld')
        if ifold: os.system('cp fdf ' + fbase + '.fdf')
        if ihbunch: os.system('cp bun ' + fbase + '.bun')
        if ifieldprop: os.system('cp fdp ' + fbase + '.fdp')
        
        
       if kplot == 2 and ifold != 1:
         
        #reakpoint()
         
        zone(2,3)
        
        #if ifieldprop == 2: zone(3,1)
        #else: zone(2,1)
        #reakpoint()
        
        if nelec > 1 and modepin == 1:

          if npinz > 1 and npiny > 1 and ihbunch > 0:
            #reakpoint()
            #npl(nbun,"z:y:s0:s0")
            hbook2('Hpin','distribution in pinhole (nbun)',max(21,npinz),zmin,zmax,\
            max(21,npiny),ymin,ymax)
            nproj2(nbun,"z:y","s0",selgam,idh='Hpin',ioverwrite=0)
            hplave('Hpin','boxes')
            txyz("S0 in Pinhole")
          else:
            npl(nebm,"z:y","","s0",plopt='boxes')
          #endif npinz, npiny
          
          #nargs = 3
          #args.append('phwig')
          
        else: #modepin
          
          if npinz > 1 and npiny > 1:
            hbook2('Hpin','distribution in pinhole',npinz,zmin,zmax,npiny,ymin,ymax)
            nproj2(nfld,"z:y","s0",selgam,idh='Hpin',ioverwrite=1)
            hplave('Hpin','boxes')
            #npl(nfld,"z:y:s0:s0")
            txyz("nfld,z:y:fd (" + selgam +")")
            rmsS0zy = hrms('Hpin',isilent=1)
            print("RMS S0,z,y:",pg5(rmsS0zy[0]),pg5(rmsS0zy[1]))
          #endif npinz, npiny
        
        #endif nbun
        
        if nexist('nfdf') and modepin != 1 and npinz > 1 and npiny > 1:
          nex()
          hbook2('HpinF','folded distribution in pinhole',npinz,zmin,zmax,npiny,ymin,ymax)
          nproj2(nfdf,"z:y","s0",selgam,idh='HpinF',ioverwrite=1)
          hplave('HpinF','boxes')
          #npl(nfld,"z:y:s0:s0")
          txyz("nfld,z:y:fdf (" + selgam +")")
        #endif npinz, npiny
        
        if nexist("nfdp"):
          nex()
          if nfdp.y.max() > nfdp.y.min(): 
            hbook2('HpinS0prop','Prop. S0',npinzprop,zminph,zmaxph,npinyprop,yminph,ymaxph)
            nproj2(nfdp,"z:y","s0",selgam,idh='HpinS0prop',ioverwrite=0)
            #hplave('HpinS0prop','inter')
            hplave('HpinS0prop','boxes')
            #npl(nfdp,"z:y","","s0")
            #npl(nfdp,"z:y:s0:s0")
            txyz("S0, Source")
            rmsS0propzy = hrms('HpinS0prop',isilent=1)
            print("RMS S0prop,z,y:",pg5(rmsS0propzy[0]),pg5(rmsS0propzy[1]))
          else:
              npl(nfdp,"z:s0",selgam)
              txyz("nfdp, " + selgam,"z","S0")
            #endif
          #endif
          
          #reakpoint()
          if ifieldprop == 2:
            nex()
            if nfdp.y.max() > nfdp.y.min(): 
              hbook2('HpinEprop','Prop. EE*',npinzprop,zminph,zmaxph,npinyprop,yminph,ymaxph)
              nproj2(nfdp,"z:y","(ezr*ezr+ezi*ezi+eyr*eyr+eyi*eyi)",selgam,idh='HpinEprop',ioverwrite=0)
              #hplave('HpinS0prop','inter')
              hplave('HpinS0prop','boxes')
              #npl(nfdp,"z:y:(ezr*ezr+ezi*ezi+eyr*eyr+eyi*eyi):(ezr*ezr+ezi*ezi+eyr*eyr+eyi*eyi)")
            else:
              npl(nfdp,"z:(ezr*ezr+ezi*ezi+eyr*eyr+eyi*eyi)")
            #endif
          elif ifieldprop == 1:
            if nexist("nfdpf"):
              nex()
              if nfdpf.y.max() > nfdpf.y.min(): 
                hbook2('HpinS0propF','PropF. S0',npinzprop,zminph,zmaxph,npinyprop,yminph,ymaxph)
                nproj2(nfdpf,"z:y","s0",selgam,idh='HpinS0propF',ioverwrite=0)
                #hplave('HpinS0prop','inter')
                hplave('HpinS0propF','boxes')
                #npl(nfdp,"z:y","","s0")
                #npl(nfdp,"z:y:s0:s0")
                txyz("nfdpf, " + selgam,"z","y","S0")
              else:
                npl(nfdpf,"z:s0")
                txyz(selgam,"z","S0")
              #endif
            #endif
          # endif prop
          #wc()
          if not ifold: pp("u51.pdf")
          else: pp("u51_fold.pdf")
          if ifieldprop == 1 and ifold == 1:
            winl()
            zone(2,2)
            hbook1('HpinS0PropZ','Prop. z:S0',npinzprop,zminph,zmaxph)
            nproj1(nfdp,"z","s0",selgam,idh='HpinS0PropZ',ioverwrite=0)
            hpl('HpinS0PropZ')
            nex()
            hbook1('HpinS0PropY','Prop. y:S0',npinyprop,yminph,ymaxph)
            nproj1(nfdp,"y","s0",selgam,idh='HpinS0PropY',ioverwrite=0)
            hpl('HpinS0PropY','L')
            #nprof(nfdp,"y:s0")
            nex()
            hbook1('HpinS0PropZF','PropF. z:S0',npinzprop,zminph,zmaxph)
            nproj1(nfdpf,"z","s0",selgam,idh='HpinS0PropZF',ioverwrite=0)
            hpl('HpinS0PropZF')
            hbook1('HpinS0PropYF','PropF. y:S0',npinyprop,yminph,ymaxph)
            nproj1(nfdpf,"y","s0",selgam,idh='HpinS0PropYF',ioverwrite=0)
            hpl('HpinS0PropYF','sameline')
          else:
            #reakpoint()
            
            hbook1('HpinS0Z','z:S0',npinz,zmin,zmax)
            nproj1(nfld,"z","s0",selgam,idh='HpinS0Z',ioverwrite=0)
            nex()
            hpl('HpinS0Z','L')
            #reakpoint()            
            txyz("HpinS0Z, proj.","z[mm]")
            
            hbook1('HpinS0Y','y:S0',npiny,ymin,ymax)
            nproj1(nfld,"y","s0",selgam,idh='HpinS0Y',ioverwrite=0)
            nex()
            hpl('HpinS0Y','L')
            txyz("HpinS0Y, proj.","y[mm]")
            
            hbook1('HpinS0PropZ','Prop. z:S0',npinzprop,zminph,zmaxph)
            nproj1(nfdp,"z","s0",selgam,idh='HpinS0PropZ',ioverwrite=0)
            
            nex()
            hpl('HpinS0PropZ','L')
            txyz("HpinS0ProZ, proj.","z[mm]")
            
            hbook1('HpinS0PropY','Prop. y:S0',npinyprop,yminph,ymaxph)            
            nproj1(nfdp,"y","s0",selgam,idh='HpinS0PropY',ioverwrite=0)
            print("RMS S0prop,z,y:",pg5(rmsS0propzy[0]),pg5(rmsS0propzy[1]))
            
            if ifold:
              xstat(0.6)
              setlinecolor('g')
              hpl('HpinS0PropY','SL')
            else:
              nex()
              hpl('HpinS0PropY','L')
              txyz("HpinS0ProY, proj.","y[mm]")
            #endif
        #endif
        
        if hexist('HpinS0Z'):
          rmsS0z = hrms('HpinS0Z',isilent=1)
          print("RMS, S0,z:",pg5(rmsS0z))
        #endif
        if hexist('HpinS0Y'):
          rmsS0y = hrms('HpinS0Y',isilent=1)
          print("RMS, S0,y:",pg5(rmsS0y))
        #endif
        if hexist('HpinS0PropZ'):
          rmsS0propz = hrms('HpinS0PropZ',isilent=1)
          print("RMS, S0prop,z:",pg5(rmsS0propz))
        #endif
        if hexist('HpinS0PropY'):
          rmsS0propy = hrms('HpinS0PropY',isilent=1)
          print("RMS, S0prop,y:",pg5(rmsS0propy))
        #endif
                  
        if not ifold: pp("u51_zy.pdf")
        else: pp("u51_zy_fold.pdf")
        
        gtit('X=' + str(pinxprop) + ', PinX=' + str(pinx) + ', PinW=' + str(pinw) + \
        ', PinH=' + str(pinh) + ', nh=' + str(npinz) + ', nv=' + str(npiny) + \
        '\n Modepin=' + str(modepin) + ', Ifold=' + str(ifold) + ', Nelec=' + str(nelec) \
        + ', Ihbunch=' + str(ihbunch))
        
        fbase = ("urad_phase_modepin_" + str(modepin) + '_ifold_' + str(ifold) + \
        '_ihbunch_' + str(ihbunch))
        
        pp(fbase + '.pdf')
        
        os.system('cp fld ' + fbase + '.fld')
        if ifold: os.system('cp fdf ' + fbase + '.fdf')
        if ihbunch: os.system('cp bun ' + fbase + '.bun')
        if ifieldprop: os.system('cp fdp ' + fbase + '.fdp')
        
        
        if ifieldprop:
          sigr0 = nrms(nfdp,"z:s0","y==0",isilent=1)
          sigr = nrms(nfdp,"z:s0",isilent=1)
          sig = sqrt(sigr**2+(sige*1000.0)**2)
          print("\nSigE:",sige*1000)
          print("SigR_Kim, Sig_kim:",pg3(sigrkim*1000),pg3(sigkim*1000))
          print("SigR0, SigR, Sig:",pg3(sigr0),pg3(sigr),pg3(sig))
          print("SigR'_Kim, Sig'_Kim:",pg3(sigrpkim*1000),pg3(sigpkim*1000))
          print("SigR'_Walker, Sig'_Walker:",pg3(sigrpwalker*1000),pg3(sigpwalker*1000))
        #endif
         
       #endif kplot
         
      #endif "ud"

      if nargs > 2:
          
        if args[2] == 'phwig' :

          kwig = 0
          if fexist('urad_phase.wig'):
            nwig = ncread("nwig","iz:iy:itz:ity:x:y:z:ty:tz:iegam:egam:wzzr:wzzi:wzyr:wzyi:wyzr:wyzi:wyyr:wyyi","urad_phase.wig")      
            ninfo(nwig)
            kwig = 1
          #endif
            
          kwkn = 0
          if fexist('urad_phase.wkn'):
            nwkn = ncread("nwkn","iz:iy:izp:iyp:iegam:egam:x:z:y:zp:yp:wzzr:wzzi:wzyr:wzyi:wyzr:wyzi:wyyr:wyyi","urad_phase.wkn")
            ninfo(nwkn)
            kwkn = 1
          #endif
          
          if kwig*kwkn:
            zone(3,1)
            npl(nwig,"z:tz","y==0 and ty==0","wzzr",'inter')
            txyz("WzzR, y==0, ty==0","z[mm]","tz[mrad]")
            nex()
            npl(nwkn,"z:zp","","wzzr")
            txyz("WzzR","z[mm]","z'[mm]")
            nex()
            npl(nwkn,"z","y==0 and yp==0 and zp==0","wzzr")
          elif kwig:
            npl(nwig,"z:tz","y==0 and ty==0","wzzr",'inter')
            txyz("WzzR, y==0, ty==0","z[mm]","tz[mrad]")
          elif kwkn:
            zone(2,1)
            npl(nwkn,"z:zp","","wzzr")
            txyz("WzzR","z[mm]","z'[mm]")
            nex()
            npl(nwkn,"z","y==0 and yp==0 and zp==0","wzzr")
          #endif
          
#          npl(nwkn,"z:y","","wzzr")
#          txyz("WzzR","z[mm]","y[mm]")

        if args[2] == 'phwigalt' :    
          
          nwig = ncread("nwig","iz:iy:itz:ity:x:y:z:ty:tz:iegam:egam:wzzr:wzzi:wzyr:wzyi:wyzr:wyzi:wyyr:wyyi","urad_phase.wig")      
          ninfo(nwig)
          
          zone(2,2)
          
          npl(nwig,"z:tz","y==0 and ty==0","wzzr",'inter')
          txyz("WzzR, y==0, ty==0","z[mm]","tz[mrad]")
          
          nex()          
          npl(nwig,"y:ty","z==0 and tz==0","wzzr",'inter')
          txyz("WzzR, z=0, z=0","y[mm]","y'[mrad]")
 
          nex()
          npl(nwig,"z:wzzr","tz==0 and y==0 and ty==0")
#          txyz("z:wzzr, z'=0, y=0, y'=0")
#          
#          nex()          
#          npl(nwig,"y:wzzr","tz==0 and z==0 and ty==0",plopt='prof')
#          txyz("y:wzzr, z'=0, z=0, y'=0")
          
        elif args[2] == 'fd':
          npll(nfld,"egam:s0")
        #endif
        
        if args[2] == 'default4' :
          optnstat()
          npll(nfdp,"z:ezr","y==0")
          npllgs(nfdp,"y:ezr","z==0")
          txyz("z, y","EzR")
          
#Baustelle      
        if args[2] == 'default3' :
           
            zone(3,2)
            optnstat()
            
            npll(nfld,"z:ezr",'y==0')
            npllgs(nfld,"y:ezr",'z==0')
            txyz("nfld","z,y","EzR")
            nex()
            npll(nfld,"z:ezi",'y==0')
            npllgs(nfld,"y:ezi",'z==0')
            txyz("nfld","z,y","EzI")
            nex()
            optstat()
            npll(nfld,"z:ezr**2+ezi**2+eyr**2+eyi**2",'y==0')
            optnstat()
            npllgs(nfld,"y:ezr**2+ezi**2+eyr**2+eyi**2",'z==0')
            txyz("nfld","z,y","S0")

            nex()
            
            npll(nfdp,"z:ezr",'y==0')
            npllgs(nfdp,"y:ezr",'z==0')
            txyz("nfdp","z,y","EzR")
            nex()
            npll(nfdp,"z:ezi",'y==0')
            npllgs(nfdp,"y:ezi",'z==0')
            txyz("nfdp","z,y","EzI")
            nex()
            optstat()
            npll(nfdp,"z:ezr**2+ezi**2+eyr**2+eyi**2",'y==0')
            optnstat()
            npllgs(nfdp,"y:ezr**2+ezi**2+eyr**2+eyi**2",'z==0')
            txyz("nfdp","z,y","S0")

            
        if args[2] == 'default2' :
            zone(3,2)
            npl(nfld,"egam:s0")
            nex()
            nplp(nbun,"egam:s0")
          
        if args[2] == 'default1' :
            
          nref0=ncread("nref0","z:s0","ref0.dat")
          nref=ncread("nref","z:s0","ref_ff.dat")
          #npll(nref0,"z:s0",color='g')
          npll(nref,"z:s0",color='g')
          
          ystat(0.5)
          #nplmgs(nbun,"z:s0",'jbun==1')
          #ystat(0.1)
          npllbs(nbun,"z:s0",'jbun==2')
          
        if args[2] == 'default' :
            
          nref0=ncread("nref0","z:s0","ref0.dat")
          nref=ncread("nref","z:s0","ref_ff.dat")
          
          #reakpoint()
          s0max = nfld.s0.max()
          emax = nfld.query("s0=="+str(s0max)).egam.max()
          selgam = "abs(egam-" + str(emax) + ")<1.0e-10"
                    
          ny = 1
          
          if nexist('nfdp'): ny += 1
          if nexist('nwig'): ny += 1
          
          #zone(2,ny)
          
          if nexist('nbun'):
              zone(3,1)
          else:
              zone(2,1)
          #endif
          
          if nexist('nebm'):
              npl(nebm,"z:y:ezr**2+ezi**2:ezr**2+ezi**2")
              txyz("Fbunch","z","y","S0")
              nex()
              npl(nebm,"z:ezr**2+ezi**2",plopt="prof")
              txyz("Fbunch","z","S0")
              nex()
              npl(nebm,"y:ezr**2+ezi**2",plopt="prof")
              txyz("Fbunch","y","S0")
          else:
              #npl(nfld,"z:y:s0:s0")
              npl(nfld,"z:y","","s0")
              txyz("Fdense","z","y","S0")
              nex()
              #npl(nfld,"z:s0",plopt="prof")
              npl(nfld,"z:s0","y==0")
              ystat(0.5)
              npllls(nref0,"z:s0")
              ystat(0.2)
              npllgs(nref,"z:s0")
              txyz("nfld, y=0","z","S0")
              #              ystat(0.2)
              #              npl(nfld,"z:ezr**2+ezi**2",plopt="Sprof",color='g')
              #              txyz("Fdense","z","S0")
              #              npl(nfld,"z:ezr",plopt="Sprof",color='g')
              #              npl(nfld,"z:ezr",plopt="prof",color='g')
              #              h=hget("HnPlot")
              #              optnstat()
              #              vplxy(h.x,h.ave,'SL')
              #              txyz("Fdense","z","EzR")
              if nexist('nbun'):
                  nex()
                  #              ystat(0.8)
                  #              npl(nbun,"z:s0","ibu==1")
                  #              nplmgs(nbun,"z:s0","ibu==2")
                  npl(nbun,"z:s0",plopt='prof')
                  ystat(0.5)
                  npllgs(nref0,"z:s0")
                  ystat(0.2)
                  npllgs(nref,"z:s0")
                  #              npl(nbun,"z:s0","y==0")
                  txyz("nbun, y==0","z","S0")
              #endif
#              optstat()
#              npl(nbun,"z:ezr",plopt='prof')
#              h=hget("HnPlot")
#              optnstat()
#              vplxy(h.x,h.ave,'SL')
#              txyz("nbun","z","EzR")
#              npl(nfld,"y:ezr**2+ezi**2",plopt="prof")
#              txyz("Fdense","y","S0")
          #endif
#         optnstat()
          if nexist('nfdp'):
              nex()
              npl(nfdp,"z:y:ezr**2+ezi**2:ezr**2+ezi**2")
              txyz("Fprop","z","y","S0")
              nex()
              npl(nfdp,"z:ezr**2+ezi**2",plopt="prof")
              txyz("Fprop","z","S0")
          #endif
          if nexist('nwig'):
              
              hbook2('HWigZZrZZPA',"Wzz zzp",nzwig,zminwig,zmaxwig,ntz,tzmin,tzmax)
              nproj2(nwig,"z:tz","wzzr**2+wzzi**2",selgam +' and y==0 and ty==0',idh='HWigZZrZZPA',ioverwrite=0)
              
              hbook2('HWigZZrZZP',"WzzR zzp",nzwig,zminwig,zmaxwig,ntz,tzmin,tzmax)
              nproj2(nwig,"z:tz","wzzr",selgam +' and y==0 and ty==0',idh='HWigZZrZZP',ioverwrite=0)
              
              hbook2('HWigZZiZZP',"WzzI zzp",nzwig,zminwig,zmaxwig,ntz,tzmin,tzmax)
              nproj2(nwig,"z:tz","wzzi",selgam +' and y==0 and ty==0',idh='HWigZZiZZP',ioverwrite=0)
          
              nex()
              hplave('HWigZZrZZP','inter')
              txyz("y=0, ty=0","z","tz","WzzR")
              nex()
              hplave('HWigZZrZZPA','inter')
              txyz("y=0, ty=0","z","tz","Wzz")
          #endif
#          npll(nbun,"z:y:ezr**2+ezi**2:ezr**2+ezi**2")
        #endif
        
        if args[2] == 's0pin':
          
          s0max = nfld.s0.max()
          emax = nfld.query("s0=="+str(s0max)).egam.max()
          selgam = "abs(egam-" + str(emax) + ")<1.0e-10"
          #reakpoint()

          if nexist(nbun):
            
            if npinz > 1 and npiny > 1:
              #reakpoint()
              #npl(nbun,"z:y:s0:s0")
              hbook2('Hpin','distribution in pinhole (nbun)',npinz,zmin,zmax,npiny,ymin,ymax)
              nproj2(nbun,"z:y","s0",selgam,idh='Hpin',ioverwrite=0)
              hplave('Hpin','boxes')
              #txyz("nbun,z:y:fd (" + selgam +")")
            #endif npinz, npiny
            
          else: #nbun
            
            if npinz > 1 and npiny > 1:
              hbook2('Hpin','distribution in pinhole',npinz,zmin,zmax,npiny,ymin,ymax)
              nproj2(nfld,"z:y","s0",selgam,idh='Hpin',ioverwrite=0)
              hplave('Hpin')
              txyz("nfld,z:y:fd (" + selgam +")")
            elif npinz > 1:
              #            hbook1('Hpin','distribution in pinhole',npinz,zmin,zmax)
              #            nproj1(nfld,"z","s0",selgam,idh='Hpin',ioverwrite=0)
              #            hplave('Hpin')
              if nexist("nbun"):
                npl(nbun,"z:s0")
              else:
                npll(nfld,"z:s0")
              #endif
              txyz(selgam,"z","S0")
            elif npiny > 1:
              hbook1('Hpin','distribution in pinhole',npiny,ymin,ymax)
              nproj1(nfld,"y","s0",selgam,idh='Hpin',ioverwrite=0)
            elif len(nbun) > 0:
              hbook2('Hpin','distribution in pinhole',npinz,zmin,zmax,npiny,ymin,ymax)
              nproj2( nbun,"z:y","s0",selgam,idh='Hpin',ioverwrite=0)
            #endif
          #endif nbun
        #endif s0pin
        
        if args[2] == 's0pinh':
          #reakpoint()            
          s0max = nfld.s0.max()
          emax = nfld.query("s0=="+str(s0max)).egam.max()
          selgam = "abs(egam-" + str(emax) + ")<1.0e-10"
          #reakpoint()
          
          if nexist("nbun"):
            npl(nbun,"z:s0","y==0")
          else:
            npll(nfld,"z:s0","y==0")
          #endif
          txyz(selgam,"z","S0")
          
        elif args[2] == 'azi' or args[2] == 'ezi':
                    
          s0max = nfld.s0.max()
          emax = nfld.query("s0=="+str(s0max)).egam.max()
          selgam = "abs(egam-" + str(emax) + ")<1.0e-10"
          
          if npinz > 1 and npiny > 1:
            hbook2('HpinPh','distribution in pinhole',npinz,zmin,zmax,npiny,ymin,ymax)
            nproj2(nfld,"z:y","ezi",selgam,idh='HpinPh',ioverwrite=0)
            hplave('HpinPh')
            txyz(selgam,"z","y","EzI")          
          elif npinz > 1:
#            hbook1('HpinPh','distribution in pinhole',npinz,zmin,zmax)
#            nproj1(nfld,"z","azr",selgam,idh='HpinPh',ioverwrite=0)
            npll(nfld,"z:ezi")
            txyz(selgam,"z","EzI")          
          elif npiny > 1:
            hbook1('HpinPh','distribution in pinhole',npiny,ymin,ymax)
            nproj1(nfld,"y","azi",selgam,idh='HpinPh',ioverwrite=0)
          elif len(nbun) > 0:
            hbook2('HpinPh','distribution in pinhole',npinz,zmin,zmax,npiny,ymin,ymax)
            nproj2( nbun,"z:y","azi",selgam,idh='HpinPh',ioverwrite=0)
          #endif
          
        elif args[2] == 'azr' or args[2] == 'ezr':
            
          s0max = nfld.s0.max()
          emax = nfld.query("s0=="+str(s0max)).egam.max()
          selgam = "abs(egam-" + str(emax) + ")<1.0e-10"
          
          if npinz > 1 and npiny > 1:
            hbook2('HpinPh','distribution in pinhole',npinz,zmin,zmax,npiny,ymin,ymax)
            nproj2(nfld,"z:y","ezr",selgam,idh='HpinPh',ioverwrite=0)
            #hplave('HpinPh')
            hpl('HpinPh')
            txyz(selgam,"z","y","EzR")          
          elif npinz > 1:
#            hbook1('HpinPh','distribution in pinhole',npinz,zmin,zmax)
#            nproj1(nfld,"z","azr",selgam,idh='HpinPh',ioverwrite=0)
            npll(nfld,"z:ezr")
            txyz(selgam,"z","EzR")          
          elif npiny > 1:
            hbook1('HpinPh','distribution in pinhole',npiny,ymin,ymax)
            nproj1(nfld,"y","azr",selgam,idh='HpinPh',ioverwrite=0)
          elif len(nbun) > 0:
            hbook2('HpinPh','distribution in pinhole',npinz,zmin,zmax,npiny,ymin,ymax)
            nproj2( nbun,"z:y","azr",selgam,idh='HpinPh',ioverwrite=0)
          #endif
                    
        elif args[2] == 's0propall':
            
          s0max = nfld.s0.max()
          emax = nfld.query("s0=="+str(s0max)).egam.max()
          selgam = "abs(egam-" + str(emax) + ")<1.0e-10"
          
          #reakpoint()
          
          if npinz > 1 and npiny > 1:
            zone(3,2)
            npl(nfld,"z:y","","s0")
#            hbook2('Hpin','distribution in pinhole',npinz,zmin,zmax,npiny,ymin,ymax)
#            nproj2(nfld,"z:y","s0",selgam,idh='Hpin',ioverwrite=0)
#            hplave('Hpin')
            txyz("nfld","z","y","S0")
            nex()
            npl(nfld,"z:y","","ezr")
            txyz("nfld","z","y","S0")
            nex()
            npl(nfld,"z","y==0","ezr")
#            hbook2('Hpin','distribution in pinhole',npinz,zmin,zmax,npiny,ymin,ymax)
#            nproj2(nfld,"z:y","s0",selgam,idh='Hpin',ioverwrite=0)
#            hplave('Hpin')
            txyz("nfld","z","y","EzR (y==0)")
            nex()
#            hbook2('HpinS0prop','Prop. S0',npinzprop,zminph,zmaxph,npinyprop,yminph,ymaxph)
#            nproj2(nfdp,"z:y","s0",selgam,idh='HpinS0prop',ioverwrite=0)
#            hplave('HpinS0prop')
            npl(nfdp,"z:y","","s0")
            txyz(selgam,"z","y","S0")
#            zone(2,2,2,"s")
            nex()
            npl(nfdp,"z:s0",selgam,plopt='prof')
#            zone(2,2,4,"s")
            nex()
            npl(nfdp,"y:s0",selgam,plopt='prof')            
          elif npinz > 1:
            npll(nfdp,"z:s0")
            txyz(selgam,"z","S0")          
          elif npiny > 1:
            hbook1('HpinS0prop','Prop. S0',npinyprop,yminph,ymaxph)
            nproj1(nfdp,"y","s0",selgam,idh='HpinS0prop',ioverwrite=0)
          elif len(nbun) > 0:
            hbook2('HpinS0prop','Prop. S0',npinzprop,zminph,zmaxph,npinyprop,yminph,ymaxph)
            nproj2( nbun,"z:y","s0",selgam,idh='HpinS0prop',ioverwrite=0)
          #endif
          
        elif args[2] == 's0prop':
          s0max = nfld.s0.max()
          emax = nfld.query("s0=="+str(s0max)).egam.max()
          selgam = "abs(egam-" + str(emax) + ")<1.0e-10"
          npl(nfdp,"z:y","","s0")
          txyz(selgam,"z","y","S0")
          
        elif args[2] == 's0proph':
          s0max = nfld.s0.max()
          emax = nfld.query("s0=="+str(s0max)).egam.max()
          selgam = "y==0 and abs(egam-" + str(emax) + ")<1.0e-10"
          npl(nfdp,"z:s0",selgam)
          txyz(selgam,"z","S0")
          
        elif args[2] == 'eprop':
                    
          s0max = nfdp.s0.max()
          emax = nfdp.query("s0=="+str(s0max)).egam.max()
          selgam = "abs(egam-" + str(emax) + ")<1.0e-10"
          
          zone(2,2)
          #reakpoint()
          
          if npinzprop > 1 and npinyprop > 1:
            hbook2('HpinEzRprop','Prop. EzR',npinzprop,zminph,zmaxph,npinyprop,yminph,ymaxph)
            nproj2(nfdp,"z:y","ezr",selgam,idh='HpinEzRprop',ioverwrite=0)
            hplave('HpinEzRprop')
            txyz(selgam,"z","y","EzR")
            nex()
            hbook2('HpinEzIprop','Prop. EzI',npinzprop,zminph,zmaxph,npinyprop,yminph,ymaxph)
            nproj2(nfdp,"z:y","ezi",selgam,idh='HpinEzIprop',ioverwrite=0)
            hplave('HpinEzIprop')
            txyz(selgam,"z","y","EzI")
            nex()
            npll(nfdp,"z:ezr","y==0")
            ystat(0.2)
            npllgs(nfdp,"z:ezi","y==0")
            nex()
            ystat(0.6)
            npll(nfdp,"y:ezr","z==0")
            ystat(0.2)
            npllgs(nfdp,"y:ezi","z==0")
          elif npinzprop > 1:
            npll(nfdp,"z:ezr")
            txyz(selgam,"z","EzR")          
            nex()
            npll(nfdp,"z:ezi")
            txyz(selgam,"z","EzI")          
          elif npinyprop > 1:
            hbook1('HpinEzRprop','Prop. EzR',npinyprop,yminph,ymaxph)
            nproj1(nfdp,"y","ezr",selgam,idh='HpinEzRprop',ioverwrite=0)
            nex()
            hbook1('HpinEzIprop','Prop. EzI',npinyprop,yminph,ymaxph)
            nproj1(nfdp,"y","ezi",selgam,idh='HpinEzIprop',ioverwrite=0)
          elif len(nbun) > 0:
            hbook2('HpinEzRprop','Prop. EzR',npinzprop,zminph,zmaxph,npinyprop,yminph,ymaxph)
            nproj2( nbun,"z:y","ezr",selgam,idh='HpinEzRprop',ioverwrite=0)
            nex()
            hbook2('HpinEzIprop','Prop. EzI',npinzprop,zminph,zmaxph,npinyprop,yminph,ymaxph)
            nproj2( nbun,"z:y","ezi",selgam,idh='HpinEzIprop',ioverwrite=0)
          #endif
          
        elif args[2] == 'wig' or args[2] == 'wigner':    
          
          nwig = ncread("nwig","iz:iy:itz:ity:x:y:z:ty:tz:iegam:egam:wzzr:wzzi:wzyr:wzyi:wyzr:wyzi:wyyr:wyyi","urad_phase.wig")      
          ninfo(nwig)
          
          zone(2,2)
          
          npl(nwig,"z:tz","y==0 and ty==0","wzzr",'inter')
          txyz("WzzR, y==0, ty==0","z[mm]","tz[mrad]")
          
          nex()          
          npl(nwig,"y:ty","z==0 and tz==0","wzzr",'inter')
          txyz("WzzR, z=0, z=0","y[mm]","y'[mrad]")
          
          nex()          
          npl(nwig,"z:tz","y==0 and ty==0","wzzr**2+wzzi**2",'inter')
          txyz("S0, y==0, ty==0","z[mm]","tz[mrad]")
          
          nex()          
          npl(nwig,"y:ty","z==0 and tz==0","wzzr**2+wzzi**2",'inter')
          txyz("S0, z=0, z=0","y[mm]","y'[mrad]")
          
          winr()
          
          zone(2,2)
          
          nprof(nwig,"z:wzzr**2+wzzi**2")
          txyz("","Z[mm]","S0")

          nex()
          nprof(nwig,"tz:wzzr**2+wzzi**2")
          txyz("","Tz[mm]","S0")
          
          nex()
          nprof(nwig,"y:wzzr**2+wzzi**2")
          txyz("","Y[mm]","S0")

          nex()
          nprof(nwig,"ty:wzzr**2+wzzi**2")
          txyz("","Ty[mm]","S0")

          winl()
          
          zone(2,2)
          
          npl(nwig,"z:wzzr","tz==0 and y==0 and ty==0")
          ndump(nwig,"z:wzzr","tz==0 and y==0 and ty==0","wig_z_wzzr.dat")
          txyz("y,ty,tz = 0","Z[mm]","Wzz")

          nex()
          npl(nwig,"tz:wzzr","z==0 and y==0 and ty==0")
          ndump(nwig,"tz:wzzr","z==0 and y==0 and ty==0","wig_tz_wzzr.dat")
          txyz("z,y,ty = 0","Tz[mm]","Wzz")
          
          nex()
          npl(nwig,"y:wzzr","ty==0 and z==0 and tz==0")
          ndump(nwig,"y:wzzr","tz==0 and z==0 and ty==0","wig_y_wzzr.dat")
          txyz("ty,z,tz = 0","Y[mm]","Wzz")

          nex()
          npl(nwig,"ty:wzzr","z==0 and y==0 and tz==0")
          ndump(nwig,"ty:wzzr","z==0 and y==0 and tz==0","wig_ty_wzzr.dat")
          txyz("z,y,tz = 0","Ty[mm]","Wzz")

        elif args[2] == 'wigzzp':
          #reakpoint()
          s0max = nfld.s0.max()
          emax = nfld.query("s0=="+str(s0max)).egam.max()
          selgam = "abs(egam-" + str(emax) + ")<1.0e-10"
                    
          hbook2('HWigZZrZZPA',"Wzz zzp",nzwig,zminwig,zmaxwig,ntz,tzmin,tzmax)
          nproj2(nwig,"z:tz","wzzr**2+wzzi**2",selgam +' and y==0 and ty==0',idh='HWigZZrZZPA',ioverwrite=0)
          
          hbook2('HWigZZrZZP',"WzzR zzp",nzwig,zminwig,zmaxwig,ntz,tzmin,tzmax)
          nproj2(nwig,"z:tz","wzzr",selgam +' and y==0 and ty==0',idh='HWigZZrZZP',ioverwrite=0)
          
          hbook2('HWigZZiZZP',"WzzI zzp",nzwig,zminwig,zmaxwig,ntz,tzmin,tzmax)
          nproj2(nwig,"z:tz","wzzi",selgam +' and y==0 and ty==0',idh='HWigZZiZZP',ioverwrite=0)
          
          zone(3,2)

          #optnstat()
          npll(nfdp,"z:ezr","y==0")
          optnstat()
          npllgs(nfdp,"y:ezr","z==0")
          txyz(selgam,"z, y","EzR")
          optstat()
          nex()
#          hplave('HWigZZrZZP','inter')
#          txyz(selgam,"z","tz","WzzR y,ty=(0,0)")
          npl(nwig,"z:wzzr","tz==0 and y==0 and ty==0")
          txyz("WzzR tz=0, y=0, ty=0","z","WzzR")
          nex()
          npl(nwig,"tz:wzzr","z==0 and y==0 and ty==0")
#          nscan(nwig,"tz:wzzr","z==0 and y==0 and ty==0")
          txyz("WzzR z=0, y=0, ty=0","tz","WzzR")
          nex()
          hplave('HWigZZrZZP','inter')
          txyz("y=0, ty=0","z","tz","WzzR")
          nex()
          hplave('HWigZZiZZP','inter')
          txyz("y=0, ty=0","z","tz","WzzI")
          nex()
          hplave('HWigZZrZZPA','inter')
          txyz("y=0, ty=0","z","tz","Wzz")
          
        elif args[2] == 'wig':
                    
          s0max = nfld.s0.max()
          emax = nfld.query("s0=="+str(s0max)).egam.max()
          selgam = "abs(egam-" + str(emax) + ")<1.0e-10"
          
          zone(2,2)
                    
          hbook2('HWigZZiZY','WzzI zy',nzwig,zminwig,zmaxwig,nywig,yminwig,ymaxwig)
          nproj2(nwig,"z:y","wzzi",selgam +' and tz==0 and ty==0',idh='HWigZZiZY',ioverwrite=0)
          hplave('HWigZZiZY')
          txyz(selgam,"z","y","WzzI Theta=(0,0)")
          
          nex()
          hbook2('HWigZYiZY','WzyI zy',nzwig,zminwig,zmaxwig,nywig,yminwig,ymaxwig)
          nproj2(nwig,"z:y","wzyi",selgam +' and tz==0 and ty==0',idh='HWigZYiZY',ioverwrite=0)
          hplave('HWigZYiZY')
          txyz(selgam,"z","y","WzyI Theta=(0,0)")
          
          nex()
          hbook2('HWigZZiTzy','WzzI Tzy',ntz,tzmin,tzmax,nty,tymin,tymax)
          nproj2(nwig,"tz:ty","wzzi",selgam +' and z==0 and y==0',idh='HWigZZiTzy',ioverwrite=0)
          hplave('HWigZZiTzy')
          txyz(selgam,"z","y","WzzI (z,y)=(0,0)")
          
          nex()
          hbook2('HWigZYiTzy','WzyI Tzy',ntz,tzmin,tzmax,nty,tymin,tymax)
          nproj2(nwig,"tz:ty","wzyi",selgam +' and z==0 and y==0',idh='HWigZYiTzy',ioverwrite=0)
          hplave('HWigZYiTzy')
          txyz(selgam,"z","y","WzyI (z,y)=(0,0)")
          
        elif args[2] == 'fd':
          npll(nfld,"egam:s0")
        #endif
        
      else:
        
        kplot = 0
        
        #nsto = ncread("nsto","x:y:z:iegam:egam:s0:s1:s2:s3","urad_phase.sto")
        if nexist("nflx") == 0:
            nflx = ncread("nflx","iegam:egam:s0:s1:s2:s3","urad_phase.flx")
        if nexist("nfld") == 0:
            nfld = ncread("nfld","x:y:z:iegam:egam:s0:s1:s2:s3:p:exr:exi:eyr:eyi:ezr:ezi:bxr:bxi:byr:byi:bzr:bzi:nx:ny:nz","urad_phase.fld")
        if nexist("nbun") == 0 and fexist('urad_phase.bun'):
            nbun = ncread("nbun","jbun:isub:ibu:bunchx:rxi1:ryi1:rzi1:ypi1:zpi1:rxin:ryin:rzin:ypin:zpin:eel:deel:x:y:z:iegam:egam:spec:s0:s1:s2:s3:p:fb28:dt:axr:axi:ayr:ayr:azr:azi","urad_phase.bun")
        
        #ninfo(nsto)
        #ninfo(nfld)
        #ninfo(nflx)
        #ninfo(nbun)
        #Quit(pg5(nfld.s0.max()))

        if kplot == 3:
          zone(2,1)
#          npl(nfld,"z:y:s0:s0")
          npl(nfld,"z:s0","y==0")
          nextzone()
#          npl(nfdp,"z:y:s0:s0")
          npl(nfdp,"z:s0","y==0")
        #endif

        if kplot == -1:
          zone(1,2)
          npll(nfld,"z:s0","abs(y)<1.0e-6")
          nextzone()
          npll(nfld,"z:ezr","abs(y)<1.0e-6")
        #endif

        if kplot == 2:
          npll(nfld,"egam:s0")
          setxstat(0.8)
          npllgs(nfld,"egam:s3")
          #nplmgs(nbun,"egam:s0","egam>0")
        #endif
        
        if kplot == 1:
          
          s0max = nfld.s0.max()
          emax = nfld.query("s0=="+str(s0max)).egam.max()
          selgam = "abs(egam-" + str(emax) + ")<1.0e-10"
          
          if npinz > 1 and npiny > 1:
            hbook2('Hpin','distribution in pinhole',npinz,zmin,zmax,npiny,ymin,ymax)
            nproj2(nfld,"z:y","s0",selgam,idh='Hpin',ioverwrite=0)
          elif npinz > 1:
            hbook1('Hpin','distribution in pinhole',npinz,zmin,zmax)
            nproj1(nfld,"z","s0",selgam,idh='Hpin',ioverwrite=0)
          elif npiny > 1:
            hbook1('Hpin','distribution in pinhole',npiny,ymin,ymax)
            nproj1(nfld,"y","s0",selgam,idh='Hpin',ioverwrite=0)
          elif len(nbun) > 0:
            hbook2('Hpin','distribution in pinhole',npinz,zmin,zmax,npiny,ymin,ymax)
            nproj2( nbun,"z:y","s0",selgam,idh='Hpin',ioverwrite=0)
          #endif
          
          if nfld.iegam.max() == 1:
            if npinz > 1 and npiny > 1:
              hbook2('HpinEzR','Ez_real in pinhole',npinz,zmin,zmax,npiny,ymin,ymax)
              nproj2(nfld,"z:y","ezr",selgam,idh='HpinEzR',ioverwrite=0)
            elif npinz > 1:
              hbook1('HpinEzR','Ez_real in pinhole',npinz,zmin,zmax)
              nproj1(nfld,"z","ezr",selgam,idh='HpinEzR',ioverwrite=0)
            elif npiny > 1:
              hbook1('HpinEzR','Ez_real in pinhole',npiny,ymin,ymax)
              nproj1(nfld,"y","ezr",selgam,idh='HpinEzR',ioverwrite=0)
            elif len(nbun) > 0:
              hbook2('HpinEzR','Ez_real in pinhole',npinz,zmin,zmax,npiny,ymin,ymax)
              nproj2( nbun,"z:y","ezr",selgam,idh='HpinEzR',ioverwrite=0)
            #endif
            if npinz > 1 and npiny > 1:
              hbook2('HpinEyR','Ey_real in pinhole',npinz,zmin,zmax,npiny,ymin,ymax)
              nproj2(nfld,"z:y","eyr",selgam,idh='HpinEyR',ioverwrite=0)
            elif npinz > 1:
              hbook1('HpinEyR','Ey_real in pinhole',npinz,zmin,zmax)
              nproj1(nfld,"z","eyr",selgam,idh='HpinEyR',ioverwrite=0)
            elif npiny > 1:
              hbook1('HpinEyR','Ey_real in pinhole',npiny,ymin,ymax)
              nproj1(nfld,"y","eyr",selgam,idh='HpinEyR',ioverwrite=0)
            elif len(nbun) > 0:
              hbook2('HpinEyR','Ey_real in pinhole',npinz,zmin,zmax,npiny,ymin,ymax)
              nproj2( nbun,"z:y","eyr",selgam,idh='HpinEyR',ioverwrite=0)
            #endif
          #endif
          
          print("\n","Eg:",emax, " eV")
          
          if kplot > 0 and npinz*npiny > 1:
            if nfld.iegam.max() > 1:
              zone(2,2)
              nprof(nfld,"egam:s0")
              txyz('nsto,"egam:s0"')
              nextzone()
              npl(nfld,"egam:s0","abs(y)<1.0e-6 and abs(z)<1.0e-6")
              txyz('nsto,"egam:s0","abs(y)<1.0e-6 and abs(z)<1.0e-6"')
              nextzone()
              npl(nflx,"egam:s0")
              txyz("nflx,'egam:s0'")
              nextzone()
              hplave('Hpin')
              txyz("nfld,z:y:fd (" + selgam +")")
            else:
              zone(2,2)
              hplave('HpinEzR')
              txyz("Ez_Real")
              nextzone()
              hplave('HpinEyR')
              txyz("Ey_Real")
              nextzone()
              hplave('Hpin')
              txyz("Flux_dens.")
              nex()
            #endif
            
          elif kplot > 0:
            if len(nbun):
              optnstat()
              npl(nbun,"egam:s0",plopt='prof')
              optstat()
              nplmgs(nbun,"egam:s0","ibu==1")
            else:
              npl(nfld,"egam:s0")
            #endif
          #endif
          
        #endif kplot > 0
      
      #endif nargs

      get_console()
      wans('Hit Q or q to quit:')
      
    elif args[1] == "overview" or args[1] == "over" or args[1] == "ov":
      
      if os.path.exists("WAVE.mhb"):
        import waveplot as w
        from waveplot import *
        WaveOverview()
      #endif os.path.exists("WAVE.mhb"):

    elif args[1] == "ndistpin" or args[1] == "pin":
      
      if os.path.exists("WAVE.mhb"):
        
        import waveplot as w
        from waveplot import *
        
        if nargs > 2:
          ndistpin(args[2])
        else:
          ndistpin()
        #endif
        
      #endif os.path.exists("WAVE.mhb"):

    elif args[1] == "ndistphase" or args[1] == "phpin" or args[1] == "pinph":
      
      if os.path.exists("WAVE.mhb"):
        
        import waveplot as w
        from waveplot import *
        
        if nargs > 2:
          ndistphase(args[2])
        else:
          ndistphase()
        #endif
        
      #endif os.path.exists("WAVE.mhb"):

    elif args[1] == "ndistphaseh" or args[1] == "phpinh" or args[1] == "pinhph":
      
      if os.path.exists("WAVE.mhb"):
        
        import waveplot as w
        from waveplot import *
        
        if nargs > 2:
          ndistphaseh(args[2])
        else:
          ndistphaseh()
        #endif
        
      #endif os.path.exists("WAVE.mhb"):

    elif args[1] == "ndistpinv" or args[1] == "pinv" or args[1] == "vcut":
      if os.path.exists("WAVE.mhb"):
        import waveplot as w
        from waveplot import *
        ndistpinv()
      #endif os.path.exists("WAVE.mhb"):

    elif args[1] == "ndistpinh" or args[1] == "pinh" or args[1] == "hcut":
      if os.path.exists("WAVE.mhb"):
        import waveplot as w
        from waveplot import *
        if nargs > 2:
          ndistpinh(args[2])
        else:
          ndistpinh()
        #endif
      #endif os.path.exists("WAVE.mhb"):

    elif args[1] == "serie":

        #if os.path.exists("WAVE.mhb"):

            #import waveplot as w
            #from waveplot import *

            fil = "wave_serie_pencil_selected.dat"
            F = open(fil,"r")
            runs0 = F.readline().strip().split()[0].strip()
            ns0 = ncread("ns0","e:s0:s1:s2:s3",fil,silent=1,skiphead=2)
            F.close()

            fil = "wave_serie_espread_selected.dat"
            F = open(fil,"r")
            runs0e = F.readline().strip().split()[0].strip()
            ns0e = ncread("ns0e","e:s0:s1:s2:s3",fil,silent=1,skiphead=2)
            F.close()

            fil = "wave_serie_emit_selected.dat"
            F = open(fil,"r")
            runs0f = F.readline().strip().split()[0].strip()
            ns0f = ncread("ns0f","e:s0:s1:s2:s3",fil,silent=1,skiphead=2)
            F.close()

            fil = "wave_serie_emit_espread_selected.dat"
            F = open(fil,"r")
            runs0ef = F.readline().strip().split()[0].strip()
            ns0ef = ncread("ns0ef","e:s0:s1:s2:s3",fil,silent=1,skiphead=2)
            F.close()

            optnstat()

            nplc(ns0,"e:s0/1.0e6")
            #s0max = ns0.s0.max()/1.0e6
            iemax,emax,s0max = ns0peak(ns0)

            legend("S0,      " + g3(s0max))

            tit = "Ideal Undulator, Real Beam, (" + runs0 + "-" + runs0ef + ")"
            txyz(tit,"E$_{ph}$ [eV]","N$_{ph}$/s/mm$^{2}$/100mA/0.1%BW")

            Fred = open("real_beam_serie_selected.dat","w")

            print("S0max:   " + g5(s0max))
            Fred.write("S0max: " + g5(s0max) + "\n")

            if nexist("ns0e"):
              nplcbs(ns0e,"e:s0/1.0e6")
              iemaxe,emaxe,s0maxe = ns0peak(ns0e)
              if iemax > 0:
                s0harme = ns0e.s0[iemax]/1.e6
              rde = s0maxe/s0max
              rdhe = s0harme/s0max
              line1 = "S0max_e, rde  : " + g3(s0maxe) + BL + g3(rde)
              line2 = "S0harm_e, rdhe: " + g3(s0harme) + BL + g3(rdhe)
              legend(line1 + "\n" + line2)
              print("\n" + line1 + "\n" + line2)
              Fred.write("\n" + line1 + "\n" + line2 + "\n")
            if nexist("ns0f"):
              nplcgs(ns0f,"e:s0/1.0e6")
              iemaxf,emaxf,s0maxf = ns0peak(ns0f)
              if iemax > 0:
                s0harmf = ns0f.s0[iemax]/1.e6
              rdf = s0maxf/s0max
              rdhf = s0harmf/s0max
              line1 = "S0max_f, rdf  : " + g3(s0maxf) + BL + g3(rdf)
              line2 = "S0harm_f, rdhf: " + g3(s0harmf) + BL + g3(rdhf)
              legend(line1 + "\n" + line2)
              print("\n" + line1 + "\n" + line2)
              Fred.write("\n" + line1 + "\n" + line2 + "\n")
            if nexist("ns0ef"):
              nplccs(ns0ef,"e:s0/1.0e6")
              iemaxef,emaxef,s0maxef = ns0peak(ns0ef)
              if iemax > 0:
                s0harmef = ns0ef.s0[iemax]/1.e6
              rdef = s0maxef/s0max
              rdhef = s0harmef/s0max
              line1 = "S0max_ef, rdef  : " + g3(s0maxef) + BL + g3(rdef)
              line2 = "S0harm_ef, rdhef: " + g3(s0harmef) + BL + g3(rdhef)
              legend(line1 + "\n" + line2)
              print("\n" + line1 + "\n" + line2)
              Fred.write("\n" + line1 + "\n" + line2 + "\n")
            #endif

            legend()
            Fred.close()

            pp("real_beam_serie_selected.pdf")

        #endif os.path.exists("WAVE.mhb"):

    elif args[1] == "ErnteAlt" or args[1] == "EA" :

      fall = "SRI22_phase-errors_Johannes.dat"
      nph = ncread("nph","key:beam:run:phsig:phwav:phdeg:s0:s0p:s0hp:s0h:rdp:rdhp:rdb",fall)

      ninfo(nph)

      optnstat()

      for btype in ['emit','espread','em+esp']:

        circ()
        phpl = "phsig:s0/s0hp"
        npl(nph,phpl,"beam=='pencil' and key=='bend_01'",legend='bend, single e-')
        nplmbs(nph,phpl,"beam=='pencil' and key=='cos_01'",legend='cos, single e-')
        nplmgs(nph,phpl,"beam=='pencil' and key=='sin_01'",legend='sin, single e-')
        nplmcs(nph,phpl,"beam=='pencil' and key=='taper_01'",legend='taper, single e-')

        bull()
        phpl = "phsig:rdb"
        sel = "beam=='" + btype + "' and key=='bend_01'"

        nplmrs(nph,phpl,sel + " and key=='bend_01'",legend='bend, beam')
        nplmbs(nph,phpl,sel + " and key=='cos_01'",legend='cos, beam')
        nplmgs(nph,phpl,sel + " and key=='sin_01'",legend='sin, beam')
        nplmcs(nph,phpl,sel + " and key=='taper_01'",legend='taper, beam')
        legend()
        txyz("Effect of Field Errors (" + btype + ")","phase error","Normalized of flux-density")

        for key in ['bend','cos','sin','taper']:
          for beam in ['pencil',btype]:

            bull()
            if beam == 'pencil': circle()

            for n in range(1,14):

              phpl = "phsig:rdb"
              if beam == 'pencil': phpl = "phsig:s0p/s0hp"

              if n < 10: keyn = key + "_0" + str(n)
              else: keyn = key + "_" + str(n)

              sel = "beam=='" + beam + "'  and key=='" + keyn + "'"

              if key == 'bend':
                nplmrs(nph,phpl,sel)
              elif key == 'cos':
                nplmbs(nph,phpl,sel)
              elif key == 'sin':
                nplmgs(nph,phpl,sel)
              elif key == 'taper':
                nplmcs(nph,phpl,sel)
              #endif

            #endfor
          #endfor
        #endfor

        fs = fall.split(".")
        fpdf = ""
        for f in fs[:-1]: fpdf += f

        pp(fpdf+"_"+btype+".pdf")
        pp(fpdf+"_"+btype+".png")

      #endfor

    elif args[1] == "Ernte" or args[1] == "E" :

      fall = "SRI22_phase-errors_Johannes.dat"
      nph = ncread("nph","key:beam:run:phsig:phwav:phdeg:s0:s0p:s0hp:s0h:rdp:rdhp:rdb",fall)

      ninfo(nph)

      optnstat()

      for btype in ['emit','espread','em+esp']:

        circ()
        phpl = "phwav*7.:s0/s0hp"
        npl(nph,phpl,"beam=='pencil' and key=='bend_01'",legend='bend, single e-')
        nplmbs(nph,phpl,"beam=='pencil' and key=='cos_01'",legend='cos, single e-')
        nplmgs(nph,phpl,"beam=='pencil' and key=='sin_01'",legend='sin, single e-')
        nplmcs(nph,phpl,"beam=='pencil' and key=='taper_01'",legend='taper, single e-')

        bull()
        phpl = "phwav*7.:rdb"
        sel = "beam=='" + btype + "' and key=='bend_01'"

        nplmrs(nph,phpl,sel + " and key=='bend_01'",legend='bend, beam')
        nplmbs(nph,phpl,sel + " and key=='cos_01'",legend='cos, beam')
        nplmgs(nph,phpl,sel + " and key=='sin_01'",legend='sin, beam')
        nplmcs(nph,phpl,sel + " and key=='taper_01'",legend='taper, beam')
        legend()
        txyz("IVUE32, 7th Harmonic, 2500 keV, Effect of Field Errors (" + btype + ")","n x phase error [rad]","Normalized of flux-density")

        for key in ['bend','cos','sin','taper']:
          for beam in ['pencil',btype]:

            bull()
            if beam == 'pencil': circle()

            for n in range(1,14):

              phpl = "phwav*7.:rdb"
              if beam == 'pencil': phpl = "phwav*7.:s0p/s0hp"

              if n < 10: keyn = key + "_0" + str(n)
              else: keyn = key + "_" + str(n)

              sel = "beam=='" + beam + "'  and key=='" + keyn + "'"

              if key == 'bend':
                nplmrs(nph,phpl,sel)
              elif key == 'cos':
                nplmbs(nph,phpl,sel)
              elif key == 'sin':
                nplmgs(nph,phpl,sel)
              elif key == 'taper':
                nplmcs(nph,phpl,sel)
              #endif

            #endfor
          #endfor
        #endfor

        fs = fall.split(".")
        fpdf = ""
        for f in fs[:-1]: fpdf += f

        pp(fpdf+"_"+btype+".pdf")
        pp(fpdf+"_"+btype+".png")

        ndump(nph,"phwav*7.:s0/s0hp:phsig:key","beam=='pencil'","SRI22_phase-errors_Johannes_pencil.dat")
        ndump(nph,"phwav*7.:rdb:phsig:key","beam=='emit'","SRI22_phase-errors_Johannes_emit.dat")
        ndump(nph,"phwav*7.:rdb:phsig:key","beam=='espread'","SRI22_phase-errors_Johannes_espread.dat")
        ndump(nph,"phwav*7.:rdb:phsig:key","beam=='em+esp'","SRI22_phase-errors_Johannes_emit_espread.dat")

      #endfor

    elif args[1] == "Johannes" or args[1] == "J":

      nsig = ncread("nsig","phsig","n_sigma_phi-corrected.dat")

      fall = "SRI22_phase-errors_Johannes.dat"
      Fall = open(fall,"w")
      Fall.write("* Type beam  run  PhSig PhErr_r PhErr_d S0  S0_pen S0H_pen S0H S0/S0_pen S0/S0H_pen S0/S0H\n")
      date = time.asctime(time.localtime(time.time()))
      Fall.write("* " + date + "\n")
      Fall.close()

      Fprot = open("serie_Johannes.pro","w")

      for key in ['bend','cos','sin','taper']:
        for n in range(1,14):

          if n < 10: keyn = key + "_0" + str(n)
          else: keyn = key + "_" + str(n)

          print(NL,keyn,NL)

          flis = key + ".lis"
          stat = os.system("rm " + flis + " 2>/dev/null")
          stat = os.system("ls -1 a/wave_*" + keyn + "* >> " + flis)
          os.system("cat " + flis)
          sleep(1)
          Flis = open(flis,"r")
          fdo = Flis.readlines()
          Flis.close()

          if len(fdo) == 4:
            nf = 0
            for f in fdo:
              nf += 1
              f = f.strip()
              for i in range(5):
                if re.search("pencil",f):
                  fs0p = f
                  runp = f.split(".")[-1]
                elif re.search("emit_espread",f):
                  fs0ef =f
                  runef = f.split(".")[-1]
                elif re.search("emit",f):
                  fs0f =f
                  runf = f.split(".")[-1]
                elif re.search("espread",f):
                  fs0e =f
                  rune = f.split(".")[-1]
                else: Quit("*** Bad key-word in " + f)
              #endfor
            #endfor

            vs0ene,vs0p,vs0f,vs0e,vs0ef = phase_error(keyn,fs0p,fs0f,fs0e,fs0ef,fall=fall)

            phdeg,phwav = get_phase_error(runp)
            write(Fprot,keyn,phdeg,phwav)

          #endif 4 files
      #endfor key

      Fprot.close()
      Quit()

    elif args[1] == "b0erroramprep":

      fserie = "/home/scheer/wav/stage/serie_amprep.out"
      namp = ncread("namp","run:iseed:key:pherr:e:s0",fserie)
      ninfo(namp)

      fserie = "/home/scheer/wav/work/serie_b0error_pencil.out"
      nb00 = ncread("nb00","run:iseed:b0err:pherr:e:s0",fserie)
      ninfo(nb00)

      fserie = "/home/scheer/wav/work/serie_b0error.out"
      nb0 = ncread("nb0","run:iseed:b0err:pherr:e:s0",fserie)
      ninfo(nb0)

      smpencil = nstat(namp,"s0","key=='pencil' and pherr==0",isilent=1)
      smemit = nstat(namp,"s0","key=='emit' and pherr==0",isilent=1)
      smespread = nstat(namp,"s0","key=='espread' and pherr==0",isilent=1)
      smemiesp = nstat(namp,"s0","key=='emiesp' and pherr==0",isilent=1)

      smbunch = nstat(nb0,"s0","b0err==0",isilent=1)
      smb00= nstat(nb00,"s0","b0err==0",isilent=1)

      optnstat()

      nharm=7
      sharmrad = str(nharm * pi/180.)

      verror=ncopv(namp,"pherr","iseed==1 and key=='pencil'")
      vwalker = exp(-(verror/180.*pi*nharm)**2)
      verrscl = verror * nharm * pi/180.
      vplxy(verrscl,vwalker,"spline",color='black',label='Walker')

      key = 'all'

      if key == 'pencil':
        circ()
        sm=str(smpencil[2])
        nplmrs(namp,sharmrad + "*pherr:s0/"+sm,"key=='pencil'",legend="s0")
        bull()
        npl(namp,sharmrad + "*pherr:s0/"+sm,"","","sameprof",legend="s0_mean + sig_mean",color='g')
      #endif

      if key == 'all':

        sm=str(smb00[2])
        bull()
        npl(nb00,sharmrad + "*pherr:s0/"+sm,"","","same",legend="pencil",color='black')

        sm=str(smbunch[2])
        npl(nb0,sharmrad + "*pherr:s0/"+sm,"","","same",legend="bunch",color='magenta')

        sm=str(smpencil[2])
        npl(namp,sharmrad + "*pherr:s0/"+sm,"key=='pencil'","","sameprof",legend="pencil",color='r')
        sm=str(smemit[2])
        circ()
        npl(namp,sharmrad + "*pherr:s0/"+sm,"key=='emit'","","sameprof",legend="emit",color='g')
        sm=str(smespread[2])
        npl(namp,sharmrad + "*pherr:s0/"+sm,"key=='espread'","","sameprof",legend="espread",color='b')
        sm=str(smemiesp[2])
        bull()
        npl(namp,sharmrad + "*pherr:s0/"+sm,"key=='emiesp'","","sameprof",legend="emit + espread",color='c')

      #endif

      legend()
      txyz("Effects of Phase Errors on the 7th Harmonic","n x PhErr [rad]","Rel. on-axis Flux-density")

    elif args[1] == "b0error":

      fserie = "serie_b0error.out"

      nb0 = ncread("nb0","run:iseed:b0err:pherr:e:s0",fserie)
      ninfo(nb0)

      optnstat()
      #nplt(nb0,"b0err:pherr","","","prof")
      nplt(nb0,"pherr:s0","","","prof")

    elif args[1] == "amprep":

      fserie = "serie_amprep.out"

      namp = ncread("namp","run:iseed:key:pherr:e:s0",fserie)
      ninfo(namp)

      smpencil = nstat(namp,"s0","key=='pencil' and pherr==0",isilent=1)
      smemit = nstat(namp,"s0","key=='emit' and pherr==0",isilent=1)
      smespread = nstat(namp,"s0","key=='espread' and pherr==0",isilent=1)
      smemiesp = nstat(namp,"s0","key=='emiesp' and pherr==0",isilent=1)

      optnstat()

      nharm=7
      sharmrad = str(nharm * pi/180.)

      verror=ncopv(namp,"pherr","iseed==1 and key=='pencil'")
      vwalker = exp(-(verror/180.*pi*nharm)**2)
      verrscl = verror * nharm * pi/180.
      vplxy(verrscl,vwalker,"spline",color='black',label='Walker')

      key = 'all'

      if key == 'pencil':
        circ()
        sm=str(smpencil[2])
        nplmrs(namp,sharmrad + "*pherr:s0/"+sm,"key=='pencil'",legend="s0")
        bull()
        npl(namp,sharmrad + "*pherr:s0/"+sm,"","","sameprof",legend="s0_mean + sig_mean",color='g')
      #endif

      if key == 'all':
        sm=str(smpencil[2])
        npl(namp,sharmrad + "*pherr:s0/"+sm,"key=='pencil'","","sameprof",legend="pencil",color='r')
        sm=str(smemit[2])
        circ()
        npl(namp,sharmrad + "*pherr:s0/"+sm,"key=='emit'","","sameprof",legend="emit",color='g')
        sm=str(smespread[2])
        npl(namp,sharmrad + "*pherr:s0/"+sm,"key=='espread'","","sameprof",legend="espread",color='b')
        sm=str(smemiesp[2])
        bull()
        npl(namp,sharmrad + "*pherr:s0/"+sm,"key=='emiesp'","","sameprof",legend="emit + espread",color='c')
      #endif

      legend()
      txyz("Effects of Phase Errors on the 7th Harmonic","n x PhErr [rad]","Rel. on-axis Flux-density")

    elif args[1] == "amprep_alt":

      fall = "amprep/amprep.dat"
      famp = open(fall,"w")

      nfiles = 0
      nkeys = 0

      for fkey in ["pencil","emit","espread","emiesp"]:

        nkeys += 1
        files = glob.glob("amprep/*" + fkey + "*selected*")
        nfiles += len(files)

        for af in files:
#          print(af)
          ff = af.split("/")[1]
          try:
            fs = ff.split("_")
            key = fs[1]
            pherr = fs[2]
            seed = fs[4].split(".")[0]
          except:
            Quit(af)
          #endtry
#          if key == 'pencil': pheall.append(float(pherr))
          fdat = open(af,"r")
          lines = fdat.readlines()
          fdat.close()
          il=0
          for l in lines:
            il += 1
            if il <= 2: continue
            famp.write(key + " " + pherr + " " + seed + " " + l)
          #endfor
        #endfor key

      #endfor fkey

      famp.close()

      namp=ncread("namp","key:phe:seed:e:s0:s1:s2:s3",fall)
      ninfo(namp)

      nseeds = int(namp.seed.max())
      emin = namp.e.min()
      pheall = ncopv(namp,"phe","e == " + str(emin) + " and seed == 1 " + " and key == 'pencil'")
      nerr = len(pheall)

      ierr0 = -1
      for ierr in range(nerr):
        if pheall[ierr] == 0:
          ierr0 = ierr
          break
        #endif
      #endfor

      if ierr0 == -1: Quit("*** Fehler: Zero error case is missing ***")

      s0maxpencil = np.zeros([nerr,nseeds])
      s0maxemit = np.zeros([nerr,nseeds])
      s0maxespread = np.zeros([nerr,nseeds])
      s0maxemiesp = np.zeros([nerr,nseeds])

      iplot = 1
      if iplot > 1:
        optnstat()
        zone(2,1)
      #endif

      iph = 0

      for ierr in range(nerr):

        spherr = str(pheall[ierr])

        for kseed in range(nseeds):

          iseed = kseed + 1
          seed = str(iseed)

          sel = "seed == " + seed + " and phe == " + spherr + " and key == 'pencil'"
          s0maxpencil[ierr,kseed] = namp.query(sel).s0.max()
          if iplot > 1:
            nplls(namp,"e:s0/1.e6",sel)
            if iph == 0: legend('pencil')
          #endif

          sel = "seed == " + seed + " and phe == " + spherr + " and key == 'emit'"
          s0maxemit[ierr,kseed] = namp.query(sel).s0.max()
          if iplot > 1:
            npllgs(namp,"e:s0/1.e6",sel)
            if iph == 0: legend('emit')
          #endif

          sel = "seed == " + seed + " and phe == " + spherr + " and key == 'espread'"
          s0maxespread[ierr,kseed] = namp.query(sel).s0.max()
          if iplot > 1:
            npllbs(namp,"e:s0/1.e6",sel)
            if iph == 0: legend('espread')
          #endif

          sel = "seed == " + seed + " and phe == " + spherr + " and key == 'emiesp'"
          s0maxemiesp[ierr,kseed] = namp.query(sel).s0.max()
#          if spherr == '0.0' or spherr == '1.0':
#            print(ierr,seed,kseed,spherr,s0maxemiesp[ierr,kseed])
          if iplot > 1:
            npllcs(namp,"e:s0/1.e6",sel)
            if iph == 0:
              legend('emit + espread')
              legend()
            #endif
          #endif

          iph += 1
          #if iph > 5: break

        #endfor iseed

      #endfor pherr

      if iplot > 1: nextzone()

      fv = open("amprep/recover.dat","w")

      if iplot:

        circ()

        verror = vcre(nerr)
        vpencil = vcre(nerr)
        vemit = vcre(nerr)
        vespread = vcre(nerr)
        vemiesp = vcre(nerr)

        iph = 0

#        print("\n")
        for iseed in range(nseeds):

          for ierr in range(nerr):

            verror[ierr] = pheall[ierr]

            vpencil[ierr] = s0maxpencil[ierr,iseed] / s0maxpencil[ierr0,iseed]
            vemit[ierr] = s0maxemit[ierr,iseed] / s0maxemit[ierr0,iseed]
            vespread[ierr] = s0maxespread[ierr,iseed] / s0maxespread[ierr0,iseed]
            vemiesp[ierr] = s0maxemiesp[ierr,iseed] / s0maxemiesp[ierr0,iseed]

            fwrite(fv,iseed,ierr,verror[ierr], \
            vpencil[ierr],vemit[ierr],vespread[ierr],vemiesp[ierr])

            if vemiesp[ierr] > 1: print(ierr,iseed,verror[ierr], \
                                        s0maxemiesp[ierr,iseed], \
                                        s0maxemiesp[ierr0,iseed])

          #endfor pherr

          if iph == 0:
            vplxy(verror,vpencil,label='pencil')
          #endif

          vplxy(verror+0.2,vemit,'same',color='g')
          if iph == 0: legend('emit')
          vplxy(verror-0.2,vespread,'same',color='b')
          if iph == 0: legend('espread')
          vplxy(verror,vemiesp,'same',color='c')
          if iph == 0:
            legend('emit + espread')
            legend()
            txyz("Effects of Phase Errors on Brilliance","phase error [degree]","rel. reduction")
          #endif

          vplxy(verror,vpencil,'same',color='r')

          iph += 1
          #if iph > 3: break
          #break

        #endfor iseed

        fv.close()
        ns0 = ncread("ns0","is:ie:err:pen:emi:esp:emiesp","amprep/recover.dat")
        ninfo(ns0)

        winr()
        optnstat()
        bull()

        verrsort=vsortx(verror)
        vwalker = exp(-(verrsort/180.*pi*7)**2)
        vplxy(verrsort,vwalker,"spline",color='black',label='Walker')

        npl(ns0,"err:pen","","","sameprof",legend='pencil')
        npl(ns0,"err:emi","","","sameprof",color='g',legend='emit')
        npl(ns0,"err:esp","","","sameprof",color='b',legend='espread')
        npl(ns0,"err:emiesp","","","sameprof",color='c',legend='emit + espread')

        legend()

        txyz("Effects of Phase Errors on Brilliance","phase error [degree]","rel. reduction")

      #endif iplot


#      sel = "phe == 0 and key == 'espread'"
#      npllb(namp,"e:s0/1.e6",sel,legend='espread')

    elif args[1] == "bend_01":

      ipencil = 567
      iemit = ipencil + 1
      iespread = ipencil + 2
      ief = ipencil + 3

      vs0ene,vs0,vs0f,vs0e,vs0ef = \
      phase_error("bend_01",
      "a/wave_bend_01_pencil_on-axis.dat." + str(ipencil),
      "a/wave_bend_01_emit_on-axis.dat." + str(iemit),
      "a/wave_bend_01_espread_on-axis.dat." + str(iespread),
      "a/wave_bend_01_emit_espread_on-axis.dat." + str(ief),
                  fall="real_beam_bend_01.dat")

    elif args[1] == "phase_bend":

      ipencil = 340
      ipencil = 482
      iemit = ipencil + 1
      iespread = ipencil + 2
      ief = ipencil + 3

      vs0ene,vs0,vs0f,vs0e,vs0ef = \
      phase_error("bend",
      "a/wave_bend_pencil_on-axis.dat." + str(ipencil),
      "a/wave_bend_emit_on-axis.dat." + str(iemit),
      "a/wave_bend_espread_on-axis.dat." + str(iespread),
      "a/wave_bend_emit_espread_on-axis.dat." + str(ief),
                  fall="real_beam_bend.dat")

    elif args[1] == "phase_bold":

        #if os.path.exists("WAVE.mhb"):

            try: nb = ncread("nb","e:s0:s1:s2:s3", \
            "a/wave_bend_pencil_on-axis.dat.340",skiphead=2)
            except: pass
            try: nbf = ncread("nbf","e:s0:s1:s2:s3", \
            "a/wave_bend_emit_on-axis.dat.341",skiphead=2)
            except: pass
            try: nbe = ncread("nbe","e:s0:s1:s2:s3", \
            "a/wave_bend_espread_on-axis.dat.342",skiphead=2)
            except: pass
            try: nbef = ncread("nbef","e:s0:s1:s2:s3", \
            "a/wave_bend_emit_espread_on-axis.dat.343",skiphead=2)
            except: pass

            optnstat()

            try: nplc(nb,"e:s0/1.0e6")
            except: pass
            try: nplcbs(nbe,"e:s0/1.0e6")
            except: pass
            try: nplcgs(nbf,"e:s0/1.0e6")
            except: pass
            try: nplccs(nbef,"e:s0/1.0e6")
            except: pass

            try: s0max = nb.s0.max()/1.0e6
            except: pass
            try: s0maxe = nbe.s0.max()/1.0e6
            except: pass
            try: s0maxf = nbf.s0.max()/1.0e6
            except: pass
            try: s0maxef = nbef.s0.max()/1.0e6
            except: pass

            try: text(0.6,0.9,"S0 Maximum: "+g3(s0max),halign='left')
            except: pass
            try:
                rde = s0maxe/s0max
                text(0.6,0.8,"Reduction e-spread: "+g3(rde),halign='left')
            except: pass
            try:
                rdf = s0maxf/s0max
                text(0.6,0.7,"Reduction emit.: "+g3(rdf),halign='left')
            except: pass
            try:
                rdef = s0maxef/s0max
                text(0.6,0.6,"Reduction both: "+g3(rdef),halign='left')
            except: pass

            txyz("Phase Error (bend.dat), Real Beam, 7th Harm.","E$_{ph}$ [eV]","N$_{ph}$/s/mm$^{2}$/100mA/0.1%BW")


            Fred = open("real_beam_bend_340-343.dat","w")

            try:
                print("S0max: " + g5(s0max))
                Fred.write("S0max: " + g5(s0max) + "\n")
            except: pass
            try:
                print("S0max_f, rdf: " + g5(s0maxf) + BL + g5(rdf))
                Fred.write("S0max_f, rdf: " + g5(s0maxf) + BL + g5(rdf) + "\n")
            except: pass
            try:
                print("S0max_e, rde: " + g5(s0maxe) + BL + g5(rde))
                Fred.write("S0max_e, rde: " + g5(s0maxe) + BL + g5(rde) + "\n")
            except: pass
            try:
                print("S0max_ef, rdef: " + g5(s0maxef) + BL + g5(rdef))
                Fred.write("S0max_ef, rdef: " + g5(s0maxef) + BL + g5(rdef) + "\n")
            except: pass

            Fred.close()

            pp("real_beam_bend_340-343.pdf")

        #endif

    elif args[1] == "phase_s":

        #if os.path.exists("WAVE.mhb"):

            try: nsi = ncread("nsi","e:s0:s1:s2:s3", \
            "a/wave_sin_pencil_selected.dat.432",skiphead=2)
            except: pass
            try: nsie = ncread("nsie","e:s0:s1:s2:s3", \
            "a/wave_sin_espread_selected.dat.343",skiphead=2)
            except: pass
            try: nsif = ncread("nsif","e:s0:s1:s2:s3", \
            "a/wave_sin_emit_selected.dat.342",skiphead=2)
            except: pass
            try: nsief = ncread("nsief","e:s0:s1:s2:s3", \
            "a/wave_sin_emit_espread_selected.dat.435",skiphead=2)
            except: pass

            optnstat()

            try: nplc(nsi,"e:s0/1.0e6")
            except: pass
            try: nplcbs(nsie,"e:s0/1.0e6")
            except: pass
            try: nplcgs(nsif,"e:s0/1.0e6")
            except: pass
            try: nplccs(nsief,"e:s0/1.0e6")
            except: pass

            try: s0max = nsi.s0.max()/1.0e6
            except: pass
            try: s0maxe = nsie.s0.max()/1.0e6
            except: pass
            try: s0maxf = nsif.s0.max()/1.0e6
            except: pass
            try: s0maxef = nsief.s0.max()/1.0e6
            except: pass

            try: text(0.6,0.9,"S0 Maximum: "+g3(s0max),halign='left')
            except: pass
            try:
                rde = s0maxe/s0max
                text(0.6,0.8,"Reduction e-spread: "+g3(rde),halign='left')
            except: pass
            try:
                rdf = s0maxf/s0max
                text(0.6,0.7,"Reduction emit.: "+g3(rdf),halign='left')
            except: pass
            try:
                rdef = s0maxef/s0max
                text(0.6,0.6,"Reduction both: "+g3(rdef),halign='left')
            except: pass

            txyz("Phase Error (sin.dat), Real Beam, 7th Harm.","E$_{ph}$ [eV]","N$_{ph}$/s/mm$^{2}$/100mA/0.1%BW")

            Fred = open("real_beam_sin_432-435.dat","w")

            try:
                print("S0max: " + g5(s0max))
                Fred.write("S0max: " + g5(s0max) + "\n")
            except: pass
            try:
                print("S0max_f, rdf: " + g5(s0maxf) + BL + g5(rdf))
                Fred.write("S0max_f, rdf: " + g5(s0maxf) + BL + g5(rdf) + "\n")
            except: pass
            try:
                print("S0max_e, rde: " + g5(s0maxe) + BL + g5(rde))
                Fred.write("S0max_e, rde: " + g5(s0maxe) + BL + g5(rde) + "\n")
            except: pass
            try:
                print("S0max_ef, rdef: " + g5(s0maxef) + BL + g5(rdef))
                Fred.write("S0max_ef, rdef: " + g5(s0maxef) + BL + g5(rdef) + "\n")
            except: pass

            Fred.close()

            pp("real_beam_342-345sin.pdf")

        #endif

    elif args[1] == "phase_c":

        #if os.path.exists("WAVE.mhb"):

            try: nsi = ncread("nsi","e:s0:s1:s2:s3", \
            "a/wave_cos_pencil_selected.dat.336",skiphead=2)
            except: pass
            try: nsie = ncread("nsie","e:s0:s1:s2:s3", \
            "a/wave_cos_espread_selected.dat.338",skiphead=2)
            except: pass
            try: nsif = ncread("nsif","e:s0:s1:s2:s3", \
            "a/wave_cos_emit_selected.dat.337",skiphead=2)
            except: pass
            try: nsief = ncread("nsief","e:s0:s1:s2:s3", \
            "a/wave_cos_emit_espread_selected_339.dat",skiphead=2)
            except: pass

            optnstat()

            try: nplc(nsi,"e:s0/1.0e6")
            except: pass
            try: nplcbs(nsie,"e:s0/1.0e6")
            except: pass
            try: nplcgs(nsif,"e:s0/1.0e6")
            except: pass
            try: nplccs(nsief,"e:s0/1.0e6")
            except: pass

            try: s0max = nsi.s0.max()/1.0e6
            except: pass
            try: s0maxe = nsie.s0.max()/1.0e6
            except: pass
            try: s0maxf = nsif.s0.max()/1.0e6
            except: pass
            try: s0maxef = nsief.s0.max()/1.0e6
            except: pass

            try: text(0.6,0.9,"S0 Maximum: "+g3(s0max),halign='left')
            except: pass
            try:
                rde = s0maxe/s0max
                text(0.6,0.8,"Reduction e-spread: "+g3(rde),halign='left')
            except: pass
            try:
                rdf = s0maxf/s0max
                text(0.6,0.7,"Reduction emit.: "+g3(rdf),halign='left')
            except: pass
            try:
                rdef = s0maxef/s0max
                text(0.6,0.6,"Reduction both: "+g3(rdef),halign='left')
            except: pass

            txyz("Phase Error (cos.dat), Real Beam, 7th Harm.","E$_{ph}$ [eV]","N$_{ph}$/s/mm$^{2}$/100mA/0.1%BW")

            Fred = open("real_beam_cos_336-339.dat","w")

            try:
                print("S0max: " + g5(s0max))
                Fred.write("S0max: " + g5(s0max) + "\n")
            except: pass
            try:
                print("S0max_f, rdf: " + g5(s0maxf) + BL + g5(rdf))
                Fred.write("S0max_f, rdf: " + g5(s0maxf) + BL + g5(rdf) + "\n")
            except: pass
            try:
                print("S0max_e, rde: " + g5(s0maxe) + BL + g5(rde))
                Fred.write("S0max_e, rde: " + g5(s0maxe) + BL + g5(rde) + "\n")
            except: pass
            try:
                print("S0max_ef, rdef: " + g5(s0maxef) + BL + g5(rdef))
                Fred.write("S0max_ef, rdef: " + g5(s0maxef) + BL + g5(rdef) + "\n")
            except: pass

            Fred.close()

            pp("real_beam_cos-336-339.pdf")

        #endif

    elif args[1] == "phase_h":

        #if os.path.exists("WAVE.mhb"):

            nh = ncread("nh","e:s0:s1:s2:s3", \
            "a/wave_halbach_pencil_selected.dat.347",silent=1,skiphead=2)
            nhf = ncread("nhf","e:s0:s1:s2:s3", \
            "a/wave_halbach_emit_selected.dat.348",silent=1,skiphead=2)
            nhe = ncread("nhe","e:s0:s1:s2:s3", \
            "a/wave_halbach_espread_selected.dat.349",silent=1,skiphead=2)
            nhef = ncread("nhef","e:s0:s1:s2:s3", \
            "a/wave_halbach_emit_espread_selected.dat.350",silent=1,skiphead=2)

            optnstat()
            nplc(nh,"e:s0/1.0e6")
            s0max = nh.s0.max()/1.0e6
            legend("S0,      " + g3(s0max) + ", 1.000")
            txyz("Ideal Undulator, Real Beam, 7th Harm.","E$_{ph}$ [eV]","N$_{ph}$/s/mm$^{2}$/100mA/0.1%BW")

            Fred = open("real_beam_halbach_347-350.dat","w")
            #Fred = open("real_beam_halbach_347-350.dat","w")

            print("S0max:   " + g5(s0max))
            Fred.write("S0max: " + g5(s0max) + "\n")

            if nexist("nhe"):
              nplcbs(nhe,"e:s0/1.0e6")
              s0maxe = nhe.s0.max()/1.0e6
              rde = s0maxe/s0max
              legend("S0_e,   " + g3(s0maxe) + ", " + g3(s0maxe/s0max))
              print("S0max_e, rde: " + g5(s0maxe) + BL + g5(rde))
              Fred.write("S0max_e, rde: " + g5(s0maxe) + BL + g5(rde) + "\n")
            if nexist("nhf"):
              nplcgs(nhf,"e:s0/1.0e6")
              s0maxf = nhf.s0.max()/1.0e6
              rdf = s0maxf/s0max
              legend("S0_f ,   " + g3(s0maxf) + ", " + g3(s0maxf/s0max))
              print("S0max_f, rdf: " + g5(s0maxf) + BL + g5(rdf))
              Fred.write("S0max_f, rdf: " + g5(s0maxf) + BL + g5(rdf) + "\n")
            if nexist("nhef"):
              nplccs(nhef,"e:s0/1.0e6")
              s0maxef = nhef.s0.max()/1.0e6
              rdef = s0maxef/s0max
              legend("S0_ef,   " + g3(s0maxef) + ", " + g3(s0maxef/s0max))
              Fred.write("S0max_ef, rdef: " + g5(s0maxef) + BL + g5(rdef) + "\n")
              print("Test:",rde*rdf/rdef)
            #endif
            legend()
            Fred.close()

            pp("real_beam_halbach_347-350.pdf")

    elif args[1] == "phase_ht":

        #if os.path.exists("WAVE.mhb"):

            nh = ncread("nh","e:s0:s1:s2:s3", \
            "wave_halba_tab_pencil_selected.dat.437",silent=1,skiphead=2)
            nhf = ncread("nhf","e:s0:s1:s2:s3", \
            "wave_halba_tab_emit_selected.dat.438",silent=1,skiphead=2)
            nhe = ncread("nhe","e:s0:s1:s2:s3", \
            "wave_halba_tab_espread_selected.dat.439",silent=1,skiphead=2)
            nhef = ncread("nhef","e:s0:s1:s2:s3", \
            "wave_halba_tab_emit_espread_selected.dat.340",silent=1,skiphead=2)

            optnstat()
            nplc(nh,"e:s0/1.0e6")
            s0max = nh.s0.max()/1.0e6
            text(0.6,0.9,"S0 Maximum: "+g3(s0max),halign='left')
            txyz("Ideal Undulator, Real Beam, 7th Harm.","E$_{ph}$ [eV]","N$_{ph}$/s/mm$^{2}$/100mA/0.1%BW")

            Fred = open("real_beam_halba_tab_437-440.dat","w")
            #Fred = open("real_beam_halbach_347-350.dat","w")

            print("S0max: " + g5(s0max))
            Fred.write("S0max: " + g5(s0max) + "\n")

            if nexist("nhe"):
              nplcbs(nhe,"e:s0/1.0e6")
              s0maxe = nhe.s0.max()/1.0e6
              rde = s0maxe/s0max
              text(0.6,0.8,"Reduction e-spread: "+g3(rde),halign='left')
              print("S0max_e, rde: " + g5(s0maxe) + BL + g5(rde))
              Fred.write("S0max_e, rde: " + g5(s0maxe) + BL + g5(rde) + "\n")
            if nexist("nhf"):
              nplcgs(nhf,"e:s0/1.0e6")
              s0maxf = nhf.s0.max()/1.0e6
              rdf = s0maxf/s0max
              text(0.6,0.7,"Reduction emit.: "+g3(rdf),halign='left')
              print("S0max_f, rdf: " + g5(s0maxf) + BL + g5(rdf))
              Fred.write("S0max_f, rdf: " + g5(s0maxf) + BL + g5(rdf) + "\n")
            if nexist("nhef"):
              nplccs(nhef,"e:s0/1.0e6")
              s0maxef = nhef.s0.max()/1.0e6
              rdef = s0maxef/s0max
              text(0.6,0.6,"Reduction both: "+g3(rdef),halign='left')
              print("S0max_ef, rdef: " + g5(s0maxef) + BL + g5(rdef))
              Fred.write("S0max_ef, rdef: " + g5(s0maxef) + BL + g5(rdef) + "\n")
            #endif

            Fred.close()

            pp("real_beam_halba_tab_347-340.pdf")

    elif args[1] == "phase_stokes":

      optnstat()
      vs0ene,vs0,vs0f,vs0e,vs0ef = \
      phase_error(fall = args[1] + ".dat")

    elif args[1] == "phase_berror":

      set_y_stat(0.2)

      fall = args[1] + ".dat"
      os.system("rm " + fall)

      vs0ene,vs0,vs0f,vs0e,vs0ef = \
      phase_error("berror",
      "wave_berror_pencil_selected.dat",
      "wave_berror_emit_selected.dat",
      "wave_berror_espread_selected.dat",
      "wave_berror_emit_espread_selected.dat",
      fall=fall)

      print("berror runs",iruns0, iruns0e, iruns0f, iruns0ef, iruns0ref)

    elif args[1] == "phase_taper":

      optnstat()
      vs0ene,vs0,vs0f,vs0e,vs0ef = \
      phase_error("taper",
      "wave_taper_pencil_selected.dat",
      "wave_taper_emit_selected.dat",
      "wave_taper_espread_selected.dat",
      "wave_taper_emit_espread_selected.dat",
      fall = args[1] + ".dat")

      print("taper runs",iruns0, iruns0e, iruns0f, iruns0ef, iruns0ref)

    elif args[1] == "phase_told":

            nt = ncread("nt","e:s0:s1:s2:s3", \
            "a/wave_taper_pencil_selected.dat.418",silent=1,skiphead=2)
            ntf = ncread("ntf","e:s0:s1:s2:s3", \
            "a/wave_taper_emit_selected.dat.419",silent=1,skiphead=2)
            nte = ncread("nte","e:s0:s1:s2:s3", \
            "a/wave_taper_espread_selected.dat.420",silent=1,skiphead=2)
            ntef = ncread("ntef","e:s0:s1:s2:s3", \
            "a/wave_taper_emit_espread_selected.dat.421",silent=1,skiphead=2)

            optnstat()
            nplc(nt,"e:s0/1.0e6")
            s0max = nt.s0.max()/1.0e6
            text(0.6,0.9,"S0 Maximum: "+g3(s0max),halign='left')
            txyz("Taper, Real Beam, 7th Harm.","E$_{ph}$ [eV]","N$_{ph}$/s/mm$^{2}$/100mA/0.1%BW")

            Fred = open("real_beam_taper_418-421.dat","w")
            #Fred = open("real_beam_taper_347-350.dat","w")

            print("S0max: " + g5(s0max))
            Fred.write("S0max: " + g5(s0max) + "\n")

            if nexist("nte"):
              nplcbs(nte,"e:s0/1.0e6")
              s0maxe = nte.s0.max()/1.0e6
              rde = s0maxe/s0max
              text(0.6,0.8,"Reduction e-spread: "+g3(rde),halign='left')
              print("S0max_e, rde: " + g5(s0maxe) + BL + g5(rde))
              Fred.write("S0max_e, rde: " + g5(s0maxe) + BL + g5(rde) + "\n")
            if nexist("ntf"):
              nplcgs(ntf,"e:s0/1.0e6")
              s0maxf = ntf.s0.max()/1.0e6
              rdf = s0maxf/s0max
              text(0.6,0.7,"Reduction emit.: "+g3(rdf),halign='left')
              print("S0max_f, rdf: " + g5(s0maxf) + BL + g5(rdf))
              Fred.write("S0max_f, rdf: " + g5(s0maxf) + BL + g5(rdf) + "\n")
            if nexist("ntef"):
              nplccs(ntef,"e:s0/1.0e6")
              s0maxef = ntef.s0.max()/1.0e6
              rdef = s0maxef/s0max
              text(0.6,0.6,"Reduction both: "+g3(rdef),halign='left')
              print("S0max_ef, rdef: " + g5(s0maxef) + BL + g5(rdef))
              Fred.write("S0max_ef, rdef: " + g5(s0maxef) + BL + g5(rdef) + "\n")
            #endif

            Fred.close()

            pp("real_beam_taper_418-421.pdf")

    elif args[1] == "bunch_h":

        if os.path.exists("WAVE.mhb"):

            nb = ncread("nb","e:s0:s1:s2:s3","wave_halbach_bunch_selected.dat",skiphead=2)
            nplc(nb,"e:s0")
            pp("real_beam_halbach_bunch.pdf")

    elif args[1] == "phase_err":

        nh = ncread("nh","e:s0:s1:s2:s3","wave_halbach__selected.dat",skiphead=2)
        if type(nh) == int: Quit()
        nhe = ncread("nhe","e:s0:s1:s2:s3","wave_halbach_e_selected.dat",skiphead=2)
        if type(nhe) == int: Quit()
        nhf = ncread("nhf","e:s0:s1:s2:s3","wave_halbach_f_selected.dat",skiphead=2)
        if type(nhf) == int: Quit()
        nhef = ncread("nhef","e:s0:s1:s2:s3","wave_halbach_ef_selected.dat",skiphead=2)
        if type(nhef) == int: Quit()

        nb = ncread("nb","e:s0:s1:s2:s3","wave_bend__selected.dat",skiphead=2)
        if type(nb) == int: Quit()
        nbe = ncread("nbe","e:s0:s1:s2:s3","wave_bend_e_selected.dat",skiphead=2)
        if type(nbe) == int: Quit()
        nbf = ncread("nbf","e:s0:s1:s2:s3","wave_bend_f_selected.dat",skiphead=2)
        if type(nbf) == int: Quit()
        nbef = ncread("nbef","e:s0:s1:s2:s3","wave_bend_ef_selected.dat",skiphead=2)
        if type(nbef) == int: Quit()

        optnstat()
        Fr = open("reduction.dat","w")

        nplc(nh,"e:s0","s0>1")
        nplcbs(nb,"e:s0","s0>1")
        text(0.6,0.8,"Max. Halbach:"+g4(nh.s0.max()),halign='left')
        text(0.6,0.75,"Max. bend.dat:"+g4(nb.s0.max()),halign='left')
        r = nb.s0.max()/nh.s0.max()
        text(0.6,0.7,"Reduction:"+g3(r),halign='left')
        print("single e: " + g5(r))
        Fr.write("single e: " + g5(r)+"\n")
        pp("reduction_single_e.pdf")

        nplc(nhe,"e:s0","s0>1")
        nplcbs(nbe,"e:s0","s0>1")
        text(0.6,0.8,"Max. Halbach:"+g4(nhe.s0.max()),halign='left')
        text(0.6,0.75,"Max. bend.dat:"+g4(nbe.s0.max()),halign='left')
        r = nbe.s0.max()/nhe.s0.max()
        text(0.6,0.7,"Reduction:"+g3(r),halign='left')
        print("espread: " + g5(r))
        Fr.write("espread: " + g5(r)+"\n")
        pp("reduction_espread.pdf")

        nplc(nhef,"e:s0","s0>1")
        nplcbs(nbef,"e:s0","s0>1")
        text(0.6,0.8,"Max. Halbach:"+g4(nhef.s0.max()),halign='left')
        text(0.6,0.75,"Max. bend.dat:"+g4(nbef.s0.max()),halign='left')
        r = nbef.s0.max()/nhef.s0.max()
        text(0.6,0.7,"Reduction:"+g3(r),halign='left')
        print("emit., espread: " + g5(r))
        Fr.write("emit. + espread: " + g5(r)+"\n")
        pp("reduction_emit_espread.pdf")
        Fr.close()

    elif args[1] == "wbmap":
        if nargs > 2:  fbm = args[2]
        else:          fbm = "wave_bmap.dat"
        nwb = ncread("nwb","x:y:z:bx:by:bz",fbm,skiphead=6)
        ninfo(nwb)

    elif args[1] == "umap":
        numap = ncread("numap","x:y:z:bx:by:bz:ifail:kfail","undumag.map",skiphead=4)
        ninfo(numap)

    elif args[1] == "s0dat":

      if nargs > 2:  fall = args[2]
      else:          fall = "wave_stokes_selected.dat"

      ns0 = ncread("ns0","e:s0:s1:s2:s3",fall)
      iemax,emax,s0max = ns0peak(ns0)
      npll(ns0,"e:s0/1.e6")

    elif args[1] == "s0datold":

      Irunmin = 1e20
      Irunmax = -1e20

      optstat()

      fil = "wave_stokes_selected.dat"
      F = open(fil,"r")
      runs0 = F.readline().strip().split()[0].strip()
      iruns0 = int(runs0)
      if iruns0 < Irunmin: Irunmin = iruns0
      if iruns0 > Irunmax: Irunmax = iruns0
      ns0 = ncread("ns0","e:s0:s1:s2:s3",fil,silent=1,skiphead=2)
      F.close()

      try:
        fil = "wave_stokese_selected.dat"
        F = open(fil,"r")
        runs0e = F.readline().strip().split()[0].strip()
        iruns0e = int(runs0e)
        if iruns0e < Irunmin: Irunmin = iruns0e
        if iruns0e > Irunmax: Irunmax = iruns0e
        ns0e = ncread("ns0e","e:s0:s1:s2:s3",fil,silent=1,skiphead=2)
        F.close()
      except: pass

      try:
        fil = "wave_stokesf_selected.dat"
        F = open(fil,"r")
        runs0f = F.readline().strip().split()[0].strip()
        iruns0f = int(runs0f)
        if iruns0f < Irunmin: Irunmin = iruns0f
        if iruns0f > Irunmax: Irunmax = iruns0f
        ns0f = ncread("ns0f","e:s0:s1:s2:s3",fil,silent=1,skiphead=2)
        F.close()
      except: pass

      try:
        fil = "wave_stokesef_selected.dat"
        F = open(fil,"r")
        runs0ef = F.readline().strip().split()[0].strip()
        iruns0ef = int(runs0ef)
        if iruns0ef < Irunmin: Irunmin = iruns0ef
        if iruns0ef > Irunmax: Irunmax = iruns0ef
        ns0ef = ncread("ns0ef","e:s0:s1:s2:s3",fil,silent=1,skiphead=2)
        F.close()
      except: pass

      nplc(ns0,"e:s0/1.0e6")
      #s0max = ns0.s0.max()/1.0e6
      iemax,emax,s0max = ns0peak(ns0)

      legend("S0,      " + g3(s0max))

      tit = "Ideal Undulator, Real Beam, (" + str(Irunmin) + "-" + str(Irunmax) + ")"
      txyz(tit,"E$_{ph}$ [eV]","N$_{ph}$/s/mm$^{2}$/100mA/0.1%BW")

      Fred = open("real_beam_folded_selected.dat","w")
      #Fred = open("real_beam_halbach_347-350.dat","w")

      print("S0max:   " + g5(s0max))
      Fred.write("S0max: " + g5(s0max) + "\n")

      if nexist("ns0e"):
        optnstat()
        nplcbs(ns0e,"e:s0/1.0e6")
        iemaxe,emaxe,s0maxe = ns0peak(ns0e)
        if iemax > 0:
          s0harme = ns0e.s0[iemax]/1.e6
        rde = s0maxe/s0max
        rdhe = s0harme/s0max
        line1 = "S0max_e, rde  : " + g3(s0maxe) + BL + g3(rde)
        line2 = "S0harm_e, rdhe: " + g3(s0harme) + BL + g3(rdhe)
        legend(line1 + "\n" + line2)
        print("\n" + line1 + "\n" + line2)
        Fred.write("\n" + line1 + "\n" + line2 + "\n")
      #endif

      if nexist("ns0f"):
        optnstat()
        nplcgs(ns0f,"e:s0/1.0e6")
        iemaxf,emaxf,s0maxf = ns0peak(ns0f)
        if iemax > 0:
          s0harmf = ns0f.s0[iemax]/1.e6
        rdf = s0maxf/s0max
        rdhf = s0harmf/s0max
        line1 = "S0max_f, rdf  : " + g3(s0maxf) + BL + g3(rdf)
        line2 = "S0harm_f, rdhf: " + g3(s0harmf) + BL + g3(rdhf)
        legend(line1 + "\n" + line2)
        print("\n" + line1 + "\n" + line2)
        Fred.write("\n" + line1 + "\n" + line2 + "\n")
      #endif

      if nexist("ns0ef"):
        optnstat()
        nplccs(ns0ef,"e:s0/1.0e6")
        iemaxef,emaxef,s0maxef = ns0peak(ns0ef)
        if iemax > 0:
          s0harmef = ns0ef.s0[iemax]/1.e6
        rdef = s0maxef/s0max
        rdhef = s0harmef/s0max
        line1 = "S0max_ef, rdef  : " + g3(s0maxef) + BL + g3(rdef)
        line2 = "S0harm_ef, rdhef: " + g3(s0harmef) + BL + g3(rdhef)
        legend(line1 + "\n" + line2)
        print("\n" + line1 + "\n" + line2)
        Fred.write("\n" + line1 + "\n" + line2 + "\n")
      #endif

      legend()
      Fred.close()

      pp("real_beam_folded.selected.pdf")

      optstat()

    elif args[1] == "phase_halbach":

      Fall = open("SRI22_phase-errors.dat","w")
  #  Type beam  run  S0*  S0_pen S0H_pen S0H S0*/S0_pen S0*/S0H_pen S0*/S0H
      Fall.close()

      vs0ene,vs0,vs0f,vs0e,vs0ef = \
      phase_error("halbach",
      "a/wave_halbach_pencil_on-axis.dat.347",
      "a/wave_halbach_emit_on-axis.dat.348",
      "a/wave_halbach_espread_on-axis.dat.349",
      "a/wave_halbach_emit_espread_on-axis.dat.350")
      Quit()

    elif args[1] == "phase_all":

      set_y_stat(0.2)

      fall = "a/SRI22_phase-errors.dat"
      Fall = open(fall,"w")
      Fall.write("* Type beam  run  S0  S0_pen S0H_pen S0H S0/S0_pen S0/S0H_pen S0/S0H\n")
      Fall.close()

      ipencil = 347
      ipencil = 476
      iemit = ipencil + 1
      iespread = ipencil + 2
      ief = ipencil + 3

      vs0ene,vs0,vs0f,vs0e,vs0ef = \
      phase_error("halbach",
      "a/wave_halbach_pencil_on-axis.dat." + str(ipencil),
      "a/wave_halbach_emit_on-axis.dat." + str(iemit),
      "a/wave_halbach_espread_on-axis.dat." + str(iespread),
      "a/wave_halbach_emit_espread_on-axis.dat." + str(ief),fall=fall)

      runs = str(Irunmin) + "-" + str(Irunmax)
      Fs0 = open("a/wave_halbach_S0_on-axis.dat." + runs,"w")
      for i in range(len(vs0)):
        Fs0.write(g5(vs0ene[i]) + BL + g5(vs0[i]) + BL \
        + g5(vs0f[i]) + BL + g5(vs0e[i]) + BL + g5(vs0ef[i]) + BL + "\n")
      #endfor
      Fs0.close()

      ipencil = 340
      ipencil = 482
      iemit = ipencil + 1
      iespread = ipencil + 2
      ief = ipencil + 3

      vs0ene,vs0,vs0f,vs0e,vs0ef = \
      phase_error("bend",
      "a/wave_bend_pencil_on-axis.dat." + str(ipencil),
      "a/wave_bend_emit_on-axis.dat." + str(iemit),
      "a/wave_bend_espread_on-axis.dat." + str(iespread),
      "a/wave_bend_emit_espread_on-axis.dat." + str(ief),fall=fall)

      runs = str(Irunmin) + "-" + str(Irunmax)
      Fs0 = open("a/wave_bend_S0_on-axis.dat." + runs,"w")
      for i in range(len(vs0)):
        Fs0.write(g5(vs0ene[i]) + BL + g5(vs0[i]) + BL \
        + g5(vs0f[i]) + BL + g5(vs0e[i]) + BL + g5(vs0ef[i]) + BL + "\n")
      #endfor
      Fs0.close()

      ipencil = 418
      ipencil = 486
      iemit = ipencil + 1
      iespread = ipencil + 2
      ief = ipencil + 3

      vs0ene,vs0,vs0f,vs0e,vs0ef = \
      phase_error("taper",
      "a/wave_taper_pencil_on-axis.dat." + str(ipencil),
      "a/wave_taper_emit_on-axis.dat." + str(iemit),
      "a/wave_taper_espread_on-axis.dat." + str(iespread),
      "a/wave_taper_emit_espread_on-axis.dat." + str(ief),fall=fall)

      runs = str(Irunmin) + "-" + str(Irunmax)
      Fs0 = open("a/wave_taper_S0_on-axis.dat." + runs,"w")
      for i in range(len(vs0)):
        Fs0.write(g5(vs0ene[i]) + BL + g5(vs0[i]) + BL \
        + g5(vs0f[i]) + BL + g5(vs0e[i]) + BL + g5(vs0ef[i]) + BL + "\n")
      #endfor
      Fs0.close()

      ipencil = 336
      ipencil = 491
      iemit = ipencil + 1
      iespread = ipencil + 2
      ief = ipencil + 3

      vs0ene,vs0,vs0f,vs0e,vs0ef = \
      phase_error("cos",
      "a/wave_cos_pencil_on-axis.dat." + str(ipencil),
      "a/wave_cos_emit_on-axis.dat." + str(iemit),
      "a/wave_cos_espread_on-axis.dat." + str(iespread),
      "a/wave_cos_emit_espread_on-axis.dat." + str(ief),fall=fall)

      runs = str(Irunmin) + "-" + str(Irunmax)
      Fs0 = open("a/wave_cos_S0_on-axis.dat." + runs,"w")
      for i in range(len(vs0)):
        Fs0.write(g5(vs0ene[i]) + BL + g5(vs0[i]) + BL \
        + g5(vs0f[i]) + BL + g5(vs0e[i]) + BL + g5(vs0ef[i]) + BL + "\n")
      #endfor
      Fs0.close()

      ipencil = 432
      ipencil = 495
      iemit = ipencil + 1
      iespread = ipencil + 2
      ief = ipencil + 3

      vs0ene,vs0,vs0f,vs0e,vs0ef = \
      phase_error("sin",
      "a/wave_sin_pencil_on-axis.dat." + str(ipencil),
      "a/wave_sin_emit_on-axis.dat." + str(iemit),
      "a/wave_sin_espread_on-axis.dat." + str(iespread),
      "a/wave_sin_emit_espread_on-axis.dat." + str(ief),fall=fall)

      runs = str(Irunmin) + "-" + str(Irunmax)
      Fs0 = open("a/wave_sin_S0_on-axis.dat." + runs,"w")
      for i in range(len(vs0)):
        Fs0.write(g5(vs0ene[i]) + BL + g5(vs0[i]) + BL \
        + g5(vs0f[i]) + BL + g5(vs0e[i]) + BL + g5(vs0ef[i]) + BL + "\n")
      #endfor
      Fs0.close()

      ipencil = 468
      ipencil = 499
      iemit = ipencil + 1
      iespread = ipencil + 2
      ief = ipencil + 3

      vs0ene,vs0,vs0f,vs0e,vs0ef = \
      phase_error("berror",
      "a/wave_berror_pencil_on-axis.dat." + str(ipencil),
      "a/wave_berror_emit_on-axis.dat." + str(iemit),
      "a/wave_berror_espread_on-axis.dat." + str(iespread),
      "a/wave_berror_emit_espread_on-axis.dat." + str(ief),fall=fall)

      runs = str(Irunmin) + "-" + str(Irunmax)
      Fs0 = open("a/wave_berror_S0_on-axis.dat." + runs,"w")
      for i in range(len(vs0)):
        Fs0.write(g5(vs0ene[i]) + BL + g5(vs0[i]) + BL \
        + g5(vs0f[i]) + BL + g5(vs0e[i]) + BL + g5(vs0ef[i]) + BL + "\n")
      #endfor
      Fs0.close()

      ipencil = 472
      ipencil = 503
      iemit = ipencil + 1
      iespread = ipencil + 2
      ief = ipencil + 3

      vs0ene,vs0,vs0f,vs0e,vs0ef = \
      phase_error("pherr5",
      "a/wave_pherr5_pencil_on-axis.dat." + str(ipencil),
      "a/wave_pherr5_emit_on-axis.dat." + str(iemit),
      "a/wave_pherr5_espread_on-axis.dat." + str(iespread),
      "a/wave_pherr5_emit_espread_on-axis.dat." + str(ief),fall=fall)

      runs = str(Irunmin) + "-" + str(Irunmax)
      Fs0 = open("a/wave_pherr5_S0_on-axis.dat." + runs,"w")
      for i in range(len(vs0)):
        Fs0.write(g5(vs0ene[i]) + BL + g5(vs0[i]) + BL \
        + g5(vs0f[i]) + BL + g5(vs0e[i]) + BL + g5(vs0ef[i]) + BL + "\n")
      #endfor
      Fs0.close()

      ipencil = 507
      iemit = ipencil + 1
      iespread = ipencil + 2
      ief = ipencil + 3

      vs0ene,vs0,vs0f,vs0e,vs0ef = \
      phase_error("berror_12deg",
      "a/wave_berror_12deg_pencil_on-axis.dat." + str(ipencil),
      "a/wave_berror_12deg_emit_on-axis.dat." + str(iemit),
      "a/wave_berror_12deg_espread_on-axis.dat." + str(iespread),
      "a/wave_berror_12deg_emit_espread_on-axis.dat." + str(ief),fall=fall)

      runs = str(Irunmin) + "-" + str(Irunmax)
      Fs0 = open("a/wave_berror_12deg_S0_on-axis.dat." + runs,"w")
      for i in range(len(vs0)):
        Fs0.write(g5(vs0ene[i]) + BL + g5(vs0[i]) + BL \
        + g5(vs0f[i]) + BL + g5(vs0e[i]) + BL + g5(vs0ef[i]) + BL + "\n")
      #endfor
      Fs0.close()

      Quit()

    elif args[1] == "66":

      n66=ncread("n66","iel:yi:zi:ypi:zpi:i:t:dph:r0x:r0y:r0z:x:y:z:ox:oy:oz:dox:doy:doz:ie:rex:iex:rda:ida","fort.66")
      n56=ncread("n56","ical:i:x:dph:ph:az:daz","fort.56")
      n67=ncread("n67","iel:yi:zi:ypi:zpi:i:t:dph:r0x:r0y:r0z:x:y:z:ox:oy:oz:dox:doy:doz:ie:rex:iex:rda:ida","~/wav/work/fort.67")
      n57=ncread("n57","ical:i:x:dph:ph:az:daz","~/wav/work/fort.57")

    elif args[1] == "stokes" or args[1] == "st":
      
      import waveplot as w
      from waveplot import *

      nsts = ncread("nsts","e:s0:s1:s2:s3","wave_stokes_selected.dat",skiphead=2)
      nst = ncread("nst","e:s0:s1:s2:s3","wave_stokes_flux.dat",skiphead=2)
      nstse = ncread("nstse","e:s0:s1:s2:s3","wave_stokese_selected.dat",skiphead=2)
      nste = ncread("nste","e:s0:s1:s2:s3","wave_stokese_flux.dat",skiphead=2)
      nstsf = ncread("nstsf","e:s0:s1:s2:s3","wave_stokesf_selected.dat",skiphead=2)
      nstf = ncread("nstf","e:s0:s1:s2:s3","wave_stokesf_flux.dat",skiphead=2)
      nstsef = ncread("nstsef","e:s0:s1:s2:s3","wave_stokesef_selected.dat",skiphead=2)
      nstef = ncread("nstef","e:s0:s1:s2:s3","wave_stokesef_flux.dat",skiphead=2)
      
      nlist()
      
      optnstat()
      sel = "e > " + str(Wflow) + " and e < " + str(Wfhig)
      #print(sel)
      
      #lilo()

      if nexist(nst):
        npllb(nst,"e:s0",sel)
        s0max = nst.s0.max()
        legend("S0,     " + g3(s0max) + ", 1.000")
      #endif
      
      #if nexist(nstf):
       # npllcs(nstf,"e:s0",sel)
        #s0fmax = nstf.s0.max()
        #legend("S0_f,  " + g3(s0fmax) + ", " + g3(s0fmax/s0max))
      #endif
      
      #if nexist(nste):
       # npllgs(nste,"e:s0",sel)
       # s0emax = nste.s0.max()
       # legend("S0_e,  " + g3(s0emax) + ", " + g3(s0emax/s0max))
      #endif
      
      if nexist(nstef):
        nplls(nstef,"e:s0",sel)
        s0efmax = nstef.s0.max()
        legend("S0_ef,  " + g3(s0efmax) + ", " + g3(s0efmax/s0max))
      #endif
      
      legend()
      
      wave_title()
      
      xTit="photon energy [eV]"      
      yTit = 'N$_{\gamma}$' + '/s/' + str(Wbw) + ' %BW/' + str(int(Wcurr*1000.+0.5)) + "mA"
      txyz(tpinhole(),xTit,yTit)
      
      pp("wave_stokesf_flux.pdf")
      
    elif args[1] == "hflux" or args[1] == "f":

        if os.path.exists("WAVE.mhb"):

            import waveplot as w
            from waveplot import *

            for n in range(len(Nhead)):
                snam = Nhead[n][1]
                exec(snam + ' = nget("' + snam + '")')
            #endfor

            hflux(clipe='no')
            optnstat()

            hs0=hget("h48000")

            s0max = hs0.y.max()
            legend("S0,     " + g3(s0max) + ", 1.000")

            #vwritexy(hs0.x,hs0.y,"s0.dat")

            if Wiefo:
              sethistcolor('g')
              hflux('s0e','same')
              hs0e=hget("h70000")
              s0emax = hs0e.y.max()
              legend("S0_e,  " + g3(s0emax) + ", " + g3(s0emax/s0max))
            #endif

            if Wifol:
              sethistcolor('cyan')
              if Wisto:
                hflux('s0f','same')
                hs0f=hget("h60000")
              else:
                hflux('ff','same')
                hs0f=hget("h49000")
              #endif
              s0fmax = hs0f.y.max()
              legend("S0_f,  " + g3(s0fmax)  + ", " + g3(s0fmax/s0max))
            #endif

            if Wifol*Wiefo:
              sethistcolor('r')
              hflux('s0ef','same')
              hs0ef=hget("h80000")
              s0efmax = hs0ef.y.max()
              legend("S0_ef, " + g3(s0efmax)  + ", " + g3(s0efmax/s0max))
            #endif

            legend()
            pp("wave_flux.pdf")
        #endif

    elif args[1] == "hcfluxden" or args[1] == "fd":

        if os.path.exists("WAVE.mhb"):
          
          import waveplot as w
          from waveplot import *
          
          for n in range(len(Nhead)):
            snam = Nhead[n][1]
            exec(snam + ' = nget("' + snam + '")')
          #endfor
          
          if Wif2p == 3.0:
            
            if Wispe:
              
              hcfluxden(clipe='no')
              optnstat()
              
              hs0=hget("h148000")
              
              s0max = hs0.y.max()
              legend("S0,     " + g3(s0max) + ", 1.000")
              
              #vwritexy(hs0.x,hs0.y,"s0.dat")
              
            else:
              print("*** No spectral data, check ISPEC in wave.in ***")
            #endif
            
            if Wiefo:
              hplot1d('h170000','same')
              hcfluxden('s0e','same')
              hs0e=hget("h170000")
              s0emax = hs0e.y.max()
              legend("S0_e,  " + g3(s0emax) + ", " + g3(s0emax/s0max))
            #endif
            
            if Wifol:
              if Wisto:
                #                hcfluxden('s0f','same')
                hs0f=hget("h160000")
              else:
                #                hcfluxden('fdf','same')
                hs0f=hget("h149000")
              #endif
              s0fmax = hs0f.y.max()
              legend("S0_f,  " + g3(s0fmax)  + ", " + g3(s0fmax/s0max))
            #endif
            
            if Wifol*Wiefo:
              # hcfluxden('s0ef','same')
              hs0ef=hget("h180000")
              s0efmax = hs0ef.y.max()
              legend("S0_ef, " + g3(s0efmax)  + ", " + g3(s0efmax/s0max))
            #endif
            
            #legend()
            #pp("wave_s0.pdf")
            
          elif Wif2p == 2.0:
            
            #lolo()
            
            npl(n3700,"ener:spec/1.e6","abs(y)<1.0e-6 and abs(z)<1.0e-6")
            
            xTit="photon energy [eV]"      
            yTit = 'N$_{\gamma}$' + '/s/' + str(Wbw) + ' %BW/' + str(int(Wcurr*1000.+0.5)) + "mA"
            txyz(tpinhole(),xTit,yTit)
              
          else:

            hcfluxden()
            
          #endif wif2p
          
        #endif WAVE.mhb
       
        ibck = 0
        icomp = 0
        if os.getcwd().split("/")[-1] == 'stage': icomp = 3

        if icomp == 4:
          nj=ncread("nj","e:fd","../job/waveplot_1.dat")
          print(n3700.spec.max()/1.e6/nj.fd.max())
          npllgs(nj)
        #endif

        icomp=0
        if icomp == 3:
          np00=ncread("np00","e:s0","pencil_00_100.dat")
          npllcs(np00,"e:s0")
          print("pencil_00: ",g5(n3700.spec.max()/1.e6/np00.s0.max()))
          text(0.85,0.9,"pencil_00: " + g5(n3700.spec.max()/1.e6/np00.s0.max()))
          np=ncread("np","e:s0","pencil_00005_100.dat")
          npllms(np,"e:s0")
          text(0.85,0.8,"pencil: " + g5(n3700.spec.max()/1.e6/np.s0.max()))
          print("pencil: ",g5(n3700.spec.max()/1.e6/np.s0.max()))
          text(0.85,0.8,"pencil: " + g5(n3700.spec.max()/1.e6/np.s0.max()))
          print("pencil_00: ",g5(n3700.spec.max()/1.e6/np.s0.max()))
          ne=ncread("ne","e:s0","emit_00_100.dat")
          npllbs(ne,"e:s0")
          text(0.85,0.7,"emit: " + g5(n3700.spec.max()/1.e6/ne.s0.max()))
          print("emit: ",g5(n3700.spec.max()/1.e6/ne.s0.max()))
          nb=ncread("nb","e:s0","beam_00_100.dat")
          npllgs(nb,"e:s0")
          text(0.85,0.6,"beam: " + g5(n3700.spec.max()/1.e6/nb.s0.max()))
          print("beam: ",g5(n3700.spec.max()/1.e6/nb.s0.max()))
        elif icomp > 0:

          if icomp == 2:
            #nbad=ncread("nbad","iel:yi:zi:ypi:zpi:i:t:dph:r0x:r0y:r0z:x:y:z:ox:oy:oz:dox:doy:doz:ie:rex:iex:rda:ida","bad.66")
            n66=ncread("n66","iel:yi:zi:ypi:zpi:i:t:dph:r0x:r0y:r0z:x:y:z:ox:oy:oz:dox:doy:doz:ie:rex:iex:rda:ida","fort.66")
            nsi=ncread("nsi","ie:t:x:z:reaz:imaz:spec:nx:nz","../job/n12.dat")
            ninfo(n66)
            ninfo(nsi)
          #endif

          nje=ncread("nje","e:s0","../job/emit_195-201.dat")
          #nje=ncread("nje","e:s0","../job/emit_195-201_1000m.dat")
          npllls(nje)
          text(0.8,0.8,g5(n3700.spec.max()/1.e6/nje.s0.max()))

          njb=ncread("njb","e:s0","../job/beam_195-201.dat")
          #njb=ncread("njb","e:s0","../job/beam_195-201_1000m.dat")
          npllgs(njb)
          text(0.8,0.7,g5(n3700.spec.max()/1.e6/njb.s0.max()))

          print("emit: ",g5(n3700.spec.max()/1.e6/nje.s0.max()))
          print("beam: ",g5(n3700.spec.max()/1.e6/njb.s0.max()))

          if ibck:
            eh = hs0.x
            s0h = hs0.y
            nbck=ncread("nbck","e:s0:s1:s2:s3","wave_stokes_selected.dat.bck")
            npllbs(nbck,"e:s0/1.e6")
            ebck = nbck.e
            s0bck = nbck.s0/1.e6
            vm = (s0bck+s0h)/2.0
            vplxy(eh,vm,"sameline",color='c')
          #endif bck

        #endif

    elif args[1] == "ray":

        if nargs > 2:  fray = args[2]
        else:          fray = "wave_ray.dat"

        waveray("hray",fray)
        Quit()

    elif args[1] == "stokesveraltet":

#        optconsole()
#        set_console_title("wavesDefault")

        if os.path.exists("WAVE.mhb"):
            import waveplot as w
            from waveplot import *
            n222=nget("n222")
            ncs0=ncread("ncs0","ener:s0:s1:s2:s3","wave_stokes__selected.dat",skiphead=2)
            ncs0f=ncread("ncs0f","ener:s0:s1:s2:s3","wave_stokesf__selected.dat",skiphead=2)
            ncs0e=ncread("ncs0e","ener:s0:s1:s2:s3","wave_stokese__selected.dat",skiphead=2)
            ncs0ef=ncread("ncs0ef","ener:s0:s1:s2:s3","wave_stokesef__selected.dat",skiphead=2)
            ns0=ncread("ns0","ener:s0:s1:s2:s3","wave_stokes_flux.dat",skiphead=2)
            ns0f=ncread("ns0f","ener:s0:s1:s2:s3","wave_stokesf_flux.dat",skiphead=2)
            ns0e=ncread("ns0e","ener:s0:s1:s2:s3","wave_stokese_flux.dat",skiphead=2)
            ns0ef=ncread("ns0ef","ener:s0:s1:s2:s3","wave_stokesef_flux.dat",skiphead=2)
            iflux=1
            if iflux == 1:
#                zone(1,2)
#                npl(ns0,"ener:s0")
                nnpl(ns0ef,"ener:s0","s0>1.")
                ene=(Wfhig+Wflow)/2.
                pp("wave_stokesef_flux_" + str(Wpinw*1000.) + "mm_x_" + str(Wpinh*1000) + "mm_" + str(ene) + "eV.pdf")
            #endif iflux == 1:
            Quit()
        #endif

    elif args[1] == "by":

#        optconsole()
#        set_console_title("wavesDefault")

        try:
          if os.path.exists("WAVE.mhb"):
            import waveplot as w
            from waveplot import *
            n222=nget("n222")
            nby()
            #endif
        except: Quit()

    elif args[1] == "b" or args[1] == "byz":

#        optconsole()
#        set_console_title("wavesDefault")
         
        import waveplot as w
        from waveplot import *
        
        try:
          if os.path.exists("WAVE.mhb"):
            n222=nget("n222")
            nbybz()
          #endif
        except:
          nlist()
        #endtry

    elif args[1] == "zy" or args[1] == "yz":

#        optconsole()
#        set_console_title("wavesDefault")

        try:
          if os.path.exists("WAVE.mhb"):
            import waveplot as w
            from waveplot import *
            n222=nget("n222")
            nzy()
            #endif
        except: Quit()

    elif args[1] == "zby" or args[1] == "zb":
      
      # optconsole()
      set_console_title("wavesDefault")
      
      import waveplot as w
      from waveplot import *
      n222=nget("n222")
      zone(1,2)
      nby()
      nextzone()
      nz()

    elif args[1] == "z":

#        optconsole()
#        set_console_title("wavesDefault")

        if os.path.exists("WAVE.mhb"):
            import waveplot as w
            from waveplot import *
            n222=nget("n222")
            nz()
        #endif

    elif args[1] == "y":

#        optconsole()
#        set_console_title("wavesDefault")

        if os.path.exists("WAVE.mhb"):
            import waveplot as w
            from waveplot import *
            n222=nget("n222")
            ny()
        #endif

    elif args[1] == "zp":

#        optconsole()
#        set_console_title("wavesDefault")

        if os.path.exists("WAVE.mhb"):
            import waveplot as w
            from waveplot import *
            n222=nget("n222")
            nzp()
        #endif

    elif args[1] == "messmodel":

        optconsole()
        set_console_title("wavesDefault")

        if os.path.exists("WAVE.mhb"):
            import waveplot as w
            from waveplot import *
            nhal=ncread("nhal","x:z:zp:by","model.dat",skiphead=1)
            nmess=ncread("nmess","x:z:zp:by","mess.dat",skiphead=1)
#            zone(1,2)
#            sel = 'abs(x+0.5)<0.1'
#            nplc(nhal,"x:by",sel)
#            nplcgs("n10","x:by",sel)
#            sel = 'abs(x+1.9)<0.1'
            zone(1,3)
            sel = 'abs(x+1.832)<0.01'
            nplc("n10","x:by",sel)
            nplcgs(nmess,"x:by",sel)
            sel = 'abs(x+1.218)<0.01'
            nnplc("n10","x:by",sel)
            nplcgs(nmess,"x:by",sel)
            sel = 'abs(x+0.535)<0.01'
            nnplc("n10","x:by",sel)
            nplcgs(nmess,"x:by",sel)
#            sel = 'abs(x+1.9)<0.1'
        #endif

    elif args[1] == "n9988":

        optconsole()
        set_console_title("wavesDefault")

        if os.path.exists("WAVE.mhb"):
            import waveplot as w
            from waveplot import *

            n99=ncread("n99","l:x:y:z:s:vx:vy:vz:b","fort.99")
            n88m=ncread("n88m","x:z:b:ang","mrad.88")
            n88w=ncread("n88w","x:z:b:ang","wave.88")
            n10=nget("n10")

            zone(1,3)

            sel="abs(x-0.55)<1."
            sell="l==9 and " + sel

            dot()
            npl(n10,"x:by",sel)
            #nplmls(n88m,"x:b","abs(x-0.6)<0.1")
            nplmls(n99,"x:b",sell)
            #nplmls(n88w,"x:b",sel)

            nextzone()
            #nplmls(n88m,"x:ang*0.017453292519943295e3","abs(x-0.6)<0.1")
            npl(n10,"x:atand(vz/vx)",sel)
            nplmgs(n99,"x:atand(vz/vx)",sell)
            #nplmgs(n88w,"x:ang*0.017453292519943295e3","abs(x-0.6)<0.1")
            #nplmgs(n99,"x:z*1000")            nextzone()

            nextzone()
            npl(n10,"x:z",sel)
            nplmgs(n99,"x:z",sell)

        else:
            print("*** WAVE.mhb not found, nothing to plot ***")
        #endif

        get_console("wavesDefault")
        wans('Hit Q or q to quit:')

    elif args[1] == "Byz":

        optconsole()
        set_console_title("wavesDefault")

        if os.path.exists("WAVE.mhb"):
            import waveplot as w
            from waveplot import *
            zone(1,2)
            nbybz()
            nextzone()
            nyz()
        else:
            print("*** WAVE.mhb not found, nothing to plot ***")
        #endif

        get_console("wavesDefault")
        wans('Hit Q or q to quit:')

    elif args[1] == "maus":
        set_console_title("wavesPython")
        nmd= ncread("nmd","x:z","maus_disp.dat")
        nms= ncread("nms","x:zp","maus_slope.dat")
    elif args[1] == "tribs_UE52":
        import waveplot as w
        from waveplot import *
        ntribs=ncread("ntribs","i:elem:s:z1:z2:z3:zp1:zp2:zp3","tribs_separation_displacement_angle.dat",skiphead=3)
        #nplc("n10","x:by","abs(x)<4.5")
#        zone(2,1)
        nplm("n10","x:z*1000.","abs(x)<4.5",color='black')
#        nplcs(ntribs,"s-75.:z3*1000.","abs(s-75.)<5") #IVUE32
        nplcrs(ntribs,"s-120.:z1*1000.","abs(s-120.)<5") #UE52
        nplcbs(ntribs,"s-120.:z2*1000.","abs(s-120.)<5") #UE52
        nplcgs(ntribs,"s-120.:z3*1000.","abs(s-120.)<5") #UE52
        txyz("UE52, TRIBS","x [m]","z [mm]")
        pp("ue52_tribs_z_vs_x.pdf")
#        nextzone()
#        ndistpow()
#        os.system("cp waveplot_3.dat ivue32_power_upstream_dipole_island-3.dat")
#        pp("ivue32_power_upstream_dipole_island-3.pdf")
    elif args[1] == "tribs_CPMU17":

        import waveplot as w
        from waveplot import *

        ntribs=ncread("ntribs","i:elem:s:z1:z2:z3:zp1:zp2:zp3","tribs_separation_displacement_angle.dat",skiphead=3)

        n10 = nget("n10")
        xmin = max(n10.x.min(),-4.6)
        xmax = n10.x.max()

        selx = "x >= " + str(xmin)
        sels = "s-165. >= " + str(xmin) + " and s-165. <= " + str(xmax)

        ms()
        nplm("n10","x:z*1000.",selx,color='black',legend='WAVE')

        nplcrs(ntribs,"s-165.:z1*1000.",sels,legend='Turn 1') #CPMU17
        nplcbs(ntribs,"s-165.:z2*1000.",sels,legend='Turn 2') #CPMU17
        nplcgs(ntribs,"s-165.:z3*1000.",sels,legend='Turn 3') #CPMU17

        legend()

        txyz("CPMU17, TRIBS","x [m]","z [mm]")
        pp("cpmu17_tribs_z_vs_x_"+ str(w.Wrun) + ".pdf")

    elif args[1] == "beta":
        import waveplot as w
        from waveplot import *
        optnstat()
        hbeta()
        #nbs=ncread("nbs","x:bx:by","~/spectra/betas.dat",skiphead=1)
        #setlinestyle("dashed")
        #nplcgs(nbs,"x:bx")
        #nplccs(nbs,"x:by")

    elif args[1] == "u41":
        nf81 = ncread("nf81","z:byi:bzi","sw2281_fit.dat")
        nw=ncread("nw","m:z:byi:bzi","fort.99")
        zone(2,1)
        nplc(nf81,"z:byi")
        nplms(nw,"z:byi","m==1")
        nextzone()
        nplc(nf81,"z:bzi")
        wans()
    elif args[1] == "u41fit":

        for i in range(9):
            com='n228'+str(i+1)+' = ncread("n228'+str(i+1)+'","z:byi:bzi","u41_field_integrals/Daten/sw228'+str(i+1)+'.yz1",skiphead=1)'
            exec(com)
        #endfor
        for i in range(2):
            com='n229'+str(i)+' = ncread("n229'+str(i)+'","z:byi:bzi","u41_field_integrals/Daten/sw229'+str(i)+'.yz1",skiphead=1)'
            exec(com)
        #endfor
        nlist()
        zone(2,1)

        nfitxy(n2282,"z/1000.:byi/1000.","",6)
        nf2282y=ncopn("Nfitxy","nf2282y")
        nf2282y.columns = ['z','byi','ey','fit']

        nfitxy(n2282,"z/1000.:bzi/1000.","",6)
        nf2282z=ncopn("Nfitxy","nf2282z")
        nf2282z.columns = ['z','bzi','ey','fit']

        n2282f = nmerge("nf2282y","nf2282z","n2282f","z:byi:fit","bzi:fit","z:byi:fy:bzi:fz")

        npl(n2282f,"z:byi")
        nplcs(n2282f,"z:fy")
        txyz("sw2282.yz1","z/m","ByInt1/Tm")

        nextzone()

        npl(n2282f,"z:bzi")
        nplcs(n2282f,"z:fz")
        txyz("sw2282.yz1","z/m","BzInt1/Tm")

        ndump(n2282f,"z:fy:fz","","sw2282_fit.dat")
        pp("sw2282_yz1.pdf")

        nextzone()
        npl(n2281,"z:byi")
        nplmbs(n2282,"z:byi")
        nplmgs(n2283,"z:byi")
        txyz("sw228-1/2/3.yz1","z/m","ByInt1/Tm")

        nextzone()
        npl(n2281,"z:bzi")
        nplmbs(n2282,"z:bzi")
        nplmgs(n2283,"z:bzi")
        txyz("sw228-1/2/3.yz1","z/m","BzInt1/Tm")
        pp("sw228_1_2_3_yz1.pdf")

    elif args[1] == "power_cpmu17":
        import waveplot as w
        from waveplot import *
        import cpmu17_power
        from cpmu17_power import *
        cpmu17_power()
    elif args[1] == "power":
        import waveplot as w
        from waveplot import *
#        null()
#        w,h=getplotsize()
#        print(w,h)
#        setplotsize(plt.gcf(),w,h)

        ndistpow()
    elif args[1] == "trib":
        #import waveplot as w
        #from waveplot import *
        set_console_title("wavesPython")
        nref = ncread('nref',"xr:yr:zr:ypr:zpr","traj_ref.dat")
        ntra = ncread('ntra',"x:y:z:yp:zp","traj_trib.dat")
        ntrib = nmerge(ntra,nref,"ntrib")
        nmd= ncread("nmd","x:z","./res_SepOrbitAngle_T6EMIL_1-turn_displacement_maus.dat")
        nms= ncread("nms","x:zp","res_SepOrbitAngle_T6EMIL_1-turn_slope_maus.dat")
        zone(1,2)
        null(nref.xr.min(),nref.xr.max(),-10.,10.)
        nplcbs(nmd,"x-165.:z")
        nplcs(ntrib,"x:z-zr")
        txyz("CPMU17 - TRIBS","x [m]","z [mm]")
        nextzone()
        null(nref.xr.min(),nref.xr.max(),-8.,8)
        nplcbs(nms,"x-165.:zp")
        nplcs(ntrib,"x:zp-zpr","x>-5.2 and x<=0")
        txyz("","x [m]","zp [mm]")
        pp("cpmu17_trib_vergleich.pdf")
        get_console()
        wans('Hit Q or q to quit:')
    elif args[1] == "track":
        import waveplot as w
        from waveplot import *
        set_console_title("wavesPython")
        optconsole()
        ndump("n20","x:by","x<0.0","traj.dat")
        Quit()
        get_console()
        wans('Hit Q or q to quit:')
    elif args[1] == "gui":
        if os.path.exists("WAVE.mhb"):
            import waveplot as w
            from waveplot import *
        else:
            print("*** WAVE.mhb not found, nothing to plot ***")
        #endif
#        set_console_title("wavesPython")
#        optconsole()
        #pp()
#endif nargs > 1

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

#try:
#  print(hmax(h148000)/3.856e12-1)
#except:
#  pass

#wans()
