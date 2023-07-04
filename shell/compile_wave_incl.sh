#!/bin/sh

# +PATCH,//WAVE/SHELL
# +DECK,compile_wave_incl,T=SHELL.

OLDPWD=$PWD

cd $WAVE_INCL/main

#rm -f $WAVE_INCL/bin/wave.exe

echo
echo Making $WAVE_INCL/bin/wave.exe
echo

rm -f *.mod

cd mod

gfortran -c -O2 -cpp \
-fcheck=bounds \
-fbacktrace \
-ffpe-summary=invalid,zero,overflow \
-fdec -fd-lines-as-comments \
-Wno-align-commons -fno-automatic -ffixed-line-length-none \
-finit-local-zero \
-funroll-loops \
*.f

mv *.mod ..

cd ..

if test -e lib/lib_mhbook.a; then . shell/mklib_mhbook_incl.sh; fi
if test -e lib/lib_mshcern.a; then . shell/mklib_mshcern_incl.sh; fi
if test -e lib/lib_mshplt.a; then . shell/mklib_mshplt_incl.sh; fi
if test -e lib/lib_user.a; then . shell/mklib_user_incl.sh; fi
if test -e lib/libwave.a; then . shell/mklib_wave_incl.sh; fi
if test -e lib/libwave_omp.a; then . shell/mklib_wave_omp_incl.sh; fi

gfortran -O2 -cpp \
-fopenmp \
-fcheck=bounds \
-fbacktrace \
-ffpe-summary=invalid,zero,overflow \
-fdec -fd-lines-as-comments \
-Wno-align-commons \
-ffixed-line-length-none \
-finit-local-zero \
-funroll-loops \
-o $WAVE_INCL/bin/wave.exe \
wave_main.f \
$WAVE_INCL/lib/libwave.a \
$WAVE_INCL/lib/libuser.a \
$WAVE_INCL/lib/libuser_modules.a \
$WAVE_INCL/lib/libmhbook.a \
$WAVE_INCL/lib/libmhbook_modules.a \
$WAVE_INCL/lib/libmshplt.a \
$WAVE_INCL/lib/libmshcern.a \
$WAVE_INCL/lib/libwave_modules.a \
$WAVE_INCL/lib/libwave_omp.a \
$WAVE_INCL/lib/libwave_omp_modules.a \
$WAVE_INCL/lib/libwave.a \
$WAVE_INCL/lib/libwave_omp.a \

cd $OLDPWD
