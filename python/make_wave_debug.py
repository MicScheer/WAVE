
# +PATCH,//WAVE/PYTHON
# +DECK,wave_make_debug,T=PYTHON.

import os
import sys
import platform
import glob

def Quit(*args, delay=0):
  #reakpoint()
  nargs =  len(args)

  text = ''
  for i in range(nargs):
    text += str(args[i]) + " "
  #endif

  if delay > 0:

    if len(text):
      print("\n",text, "\nWaiting",delay," seconds before kill")
      #time.sleep(delay)
    else:
      print("\nWaiting",delay," seconds before kill")
      #time.sleep(delay)
    #endif len(text):

    if platform.system() == 'Windows':
      os.system("sleep " + str(delay) + " && taskkill /F /PID " + str(os.getpid()) + " &")
    else:
      os.system("sleep " + str(delay) + " && kill " + str(os.getpid()) + " &")
    #endif platform.system() == 'Windows'

  elif delay < 0:
    return
  else:
    print("\n",text)
    if platform.system() == 'Windows':
      os.system("taskkill /F /PID " + str(os.getpid()))
    else:
      os.system("kill " + str(os.getpid()))
    #endif platform.system() == 'Windows'

#enddef Quit(text = '', delay=0)

global Iverbose,Idry,Idebug,WI

args=sys.argv; nargs = len(args)

try:
  WI = os.environ['WAVE_INCL'] + "/"
except:
  WI = ''
  path = args[0].split("/")
  l = len(path)
  if l == 1:
    path = os.getcwd().split("/")
    path.append(args[0])
  elif l == 2:
    path = os.getcwd().split("/")
    pp = args[0].split("/")
    path.append(pp[0])
    path.append(pp[1])
  #endif
  for i in range(len(path)-2):
    WI += path[i] + "/"
  #endfor
  print("\n*** Warning: Shell variable WAVE_INCL not defined ***")
  print("*** Assuming: ",WI," ***")
  os.system('sleep 3')
#endtry


Iverbose = 0
Idebug = 0
Idry = 0

if nargs > 1:
  try:
    Iverbose = int(args[1])
  except:
    n = '\n'
    print(n)
    print("Usage: python3 " + WI + args[0] + " [verbose level]",n)
    print("To force total recompilation delete ",n,WI + "bin/wave_debug.exe",n)
    Quit()
  #end try
#endif

if nargs > 2: Idebug = int(args[2])

global Wave_tree,Scomp_all,Scomp_omp,Scomp,Texe,Tlib,Scomp_nowarn

Scomp = "gfortran -std=legacy -c -g -cpp -fbacktrace -ffpe-summary=invalid,zero,overflow -fdec -fd-lines-as-comments -Wno-align-commons -fno-automatic -ffixed-line-length-none -finit-local-zero -funroll-loops "
Scomp_nowarn = "gfortran -w -std=legacy -c -g -cpp -fbacktrace -ffpe-summary=invalid,zero,overflow -fdec -fd-lines-as-comments -Wno-align-commons -fno-automatic -ffixed-line-length-none -finit-local-zero -funroll-loops "
Scomp_all = "gfortran -std=legacy -c -g -cpp -fcheck=all -fbacktrace -ffpe-summary=invalid,zero,overflow -fdec -fd-lines-as-comments -Wno-align-commons -fno-automatic -ffixed-line-length-none -finit-local-zero -funroll-loops "
Scomp_omp = "gfortran -std=legacy -c -g -cpp -finit-local-zero -fcheck=all -fopenmp -fbacktrace -ffpe-summary=invalid,zero,overflow -fdec -fd-lines-as-comments -Wno-align-commons -ffixed-line-length-none -funroll-loops "

def get_wave_tree():

  global WI,Wave_tree,Iverbose,Idry,Idebug,Texe,Tlib

  try:
    Texe = os.stat(WI + '/bin/wave_debug.exe').st_mtime_ns
  except:
    Texe = 0
  #endtry

  top = glob.glob(WI+"/*")

  Wave_tree = []
  #reakpoint()

  for topd in top:

    dd = topd.split("/")[-1]

    if dd == 'cmz' or dd == 'doc' or dd == 'check_system' or dd == 'bin' \
    or dd == 'python' or dd == 'main' or dd == 'lib': continue

    t = os.stat(topd).st_mtime_ns

    modf = glob.glob(topd+"/mod/*.f")

    modfor = []
    for ff in modf:
      f = ff.split("/")[-1]
      tf = os.stat(ff).st_mtime_ns
      modfor.append([f,tf])
    #endfor

    modm = glob.glob(topd+"/*.mod")
    modmod = []
    for ff in modm:
      f = ff.split("/")[-1]
      tf = os.stat(ff).st_mtime_ns
      modmod.append([f,tf])
    #endfor

    cm = glob.glob(topd+"/*.cmn")
    cmn = []
    for ff in cm:
      f = ff.split("/")[-1]
      tf = os.stat(ff).st_mtime_ns
      cmn.append([f,tf])
    #endfor

    ff = glob.glob(topd+"/*.f")
    fort = []
    for fff in ff:
      f = fff.split("/")[-1]
      tf = os.stat(fff).st_mtime_ns
      fort.append([f,tf])
    #endfor

    Wave_tree.append([topd,t,modfor,modmod,cmn,fort])

  #endfor get_wave_tree

#enddef get_wave_tree

def wave_update():

  global WI,Wave_tree,Texe,Scomp_all,Scomp_omp,Scomp,Iverbose,Idry,Idebug,Scomp_nowarn

  kmain = 0

  get_wave_tree()

  for td in Wave_tree:

    dd = td[0]
    ds = dd + "/"
    dsm = dd + "/mod/"
    t = td[1]
    modfor = td[2]
    cmn = td[4]
    fort = td[5]

    scomp = Scomp

    lib = ''
    libm = ''
    ranl = 0
    ranlm = 0
    slibm = ''
    slib = ''

    ddd = dd.split("/")[-1]

    if Iverbose >= 0: print("\nProcessing",dd)
    #reakpoint()

    if ddd == 'mhbook':
      lib = WI + 'lib/libmhbook_debug.a'
      libm = WI + 'lib/libmhbook_modules_debug.a'
    elif ddd == 'mshcern':
      lib = WI + 'lib/libmshcern_debug.a'
      libm = WI + 'lib/libmshcern_modules_debug.a'
      scomp = Scomp_nowarn
      #breakpoint()
    elif ddd == 'mshplt':
      lib = WI + 'lib/libmshplt_debug.a'
      libm = WI + 'lib/libmshplt_modules.a'
    elif ddd == 'nomp':
      lib = WI + 'lib/libwave_debug.a'
      libm = WI + 'lib/libwave_modules_debug.a'
      scomp = Scomp_all
    elif ddd == 'omp':
      lib = WI + 'lib/libwave_omp_debug.a'
      libm = WI + 'lib/libwave_omp_modules_debug.a'
      scomp = Scomp_omp
    elif ddd == 'urad':
      lib = WI + 'lib/liburad_debug.a'
      libm = WI + 'lib/liburad_modules_debug.a'
      scomp = Scomp_all  # uradcfft does boundary tricks
    elif ddd == 'user':
      lib = WI + 'lib/libuser_debug.a'
      libm = WI + 'lib/libuser_modules_debug.a'
      scomp = Scomp_omp
    #endif

    klib = 0
    try:
      Tlib = os.stat(lib).st_mtime_ns
      if Tlib > Texe: klib = 1
    except:
      klib = 1
    #endtry

    klibm = 0
    try:
      Tlib = os.stat(libm).st_mtime_ns
      if Tlib > Texe: klibm = 1
    except:
      klibm = 1
    #endtry

    scompmod = "cd " + dd + "/mod && " + scomp
    scomp = "cd " + dd + " && " + scomp

    for f in modfor: # Compile modules

      ff = f[0]
      t = f[1]

      if t < Texe and klibm == 0: continue

      if Iverbose > 0: print(ff)

      fo = ff[:-1] + "o"
      fm = ff[:-1] + "mod"

      Flines = open(ds+"mod/"+ff,'r')

      while True:
        l = Flines.readline()
        if not l: break
        sl = l.split()
        if len(sl) == 0: continue
        key = sl[0].lower()
        if key== 'module':
          m = sl[1].lower()
          break
        #endif
      #end while
      Flines.close()

      if Iverbose > 0: print("\nModule:",m)

      scom = scompmod + "-o " + fo + " " + ff
      if Iverbose > 0: print("\n",scom,"\n")
      if Idry == 0: os.system(scom)

      scom = 'mv ' + dsm + m + ".mod " + dd
      if Iverbose > 0: print("\n",scom,"\n")
      if Idry == 0: os.system(scom)

      slibm += " " + dsm + fo
      ranlm = 1

      # Search use of module in *.cmn
      for ft in cmn:

        f = ft[0]
        t = ft[1]

        if t < Texe and klibm == 0: continue

        Flines = open(ds+f,'r')
        while True:
          l = Flines.readline()
          if Idebug > 1: print(l)
          if not l: break
          #if len(l) < 10: break
          sl = l.split()
          if len(sl) > 1:
            key = sl[0].lower()
            if key== 'use':
              if sl[1].lower() == m:
                scom = 'touch ' + ds+f
                if Iverbose > 0: print("\n",scom,"\n")
                if Idry == 0: os.system(scom)
                break
              #endif
            #endif
          #endif
        #end while
        Flines.close()

      #endfor

      # Search use of module in *.f
      for ft in fort:

        f = ft[0]
        t = ft[1]

        if Idebug > 1: print(f)

        Flines = open(ds+f,'r')
        while True:
          l = Flines.readline()
          if not l: break
          #if len(l) < 10: break
          sl = l.split()
          if len(sl) > 1:
            key = sl[0].lower()
            if key== 'implicit':
              if sl[1].lower() == 'none': break
            elif key== 'use':
              if sl[1].lower() == m:
                scom = 'touch ' + ds+f
                if Iverbose > 0: print("\n",scom,"\n")
                if Idry == 0: os.system(scom)
                break
              #endif
            #endif
          #endif
        #end while
        Flines.close()

      #endfor

    #endfor modfor

    if ranlm:
      scom = 'ar rc ' + libm + " " + slibm
      if Iverbose > 0: print("\n",scom,"\n")
      if Idry == 0: os.system(scom)
      scom = 'ranlib ' + libm
      if Iverbose > 0: print("\n",scom,"\n")
      if Idry == 0: os.system(scom)
      ranlm = 0
      slibm = ''
      kmain = 1
    #endif

    # Check *.cmn

    for ft in cmn:

      f = ft[0]
      t = os.stat(ds+f).st_mtime_ns

      if t < Texe and klib == 0: continue

      fcmn = f.split("/")[-1]

      for fft in fort:

        Flines = open(ds+fft[0],'r')
        while True:
          l = Flines.readline()
          if not l: break
          #if len(l) < 10: break
          sl = l.split()
          if len(sl) < 2: continue
          if sl[0][0] == '*' or sl[0][0] == '!' or len(sl[0]) < 7: continue
          key = sl[0].lower()
          if key== 'include':
            if sl[1].lower() == "'" + fcmn + "'" or sl[1].lower() == '"' + fcmn + '"':
              scom = 'touch ' + ds+fft[0]
              if Iverbose > 0: print("\n",scom,"\n")
              if Idry == 0: os.system(scom)
              break
            #endif
          #endif
        #end while
        Flines.close()
      #endfor fort

    #endfor cmn

    # Compile *.f if neccessary

    for ft in fort:

      f = ft[0]
      t = os.stat(ds+f).st_mtime_ns

      if t < Texe and klib == 0: continue

      fo = f[:-1] + "o"

      scom = scomp + "-o " + fo + " " + f
      if Iverbose > 0: print("\n",scom,"\n")
      if Idry == 0: os.system(scom)

      slib += " " + ds + fo
      ranl = 1

    #endfor

    if ranl:
      scom = 'ar rc ' + lib + " " + slib
      if Iverbose > 0: print("\n",scom,"\n")
      if Idry == 0: os.system(scom)
      scom = 'ranlib ' + lib
      if Iverbose > 0: print("\n",scom,"\n")
      if Idry == 0: os.system(scom)
      kmain = 1
      ranl = 0
      slib = ''
    #endif

  #endfor dir

  if kmain:
    scom = WI + "/shell/compile_wave_incl_debug.sh"
    if Iverbose > 0: print("\n",scom,"\n")
    if Idry == 0: os.system(scom)
    if Iverbose >=0: print("\n--- " + WI  + "bin/wave_debug.exe updated ---\n")
  else:
    if Iverbose >=0: print("\n--- No need to update " + WI  + "bin/wave_debug.exe ---\n")
  #endif


#enddef wave_update

wave_update()
