#!/bin/sh

# +PATCH,//WAVE/SHELL
# +DECK,compile_wave_incl_debug,T=SHELL.

OLDPWD=$PWD

cd $WAVE_INCL/main

#rm -f $WAVE_INCL/bin/wave_debug.exe

echo
echo Making $WAVE_INCL/bin/wave_debug.exe
echo

rm -f *.mod

cd mod

gfortran -c -g -cpp \
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

if test -e lib/lib_mhbook_debug.a; then . shell/mklib_mhbook_incl_debug.sh; fi
if test -e lib/lib_mshcern_debug.a; then . shell/mklib_mshcern_incl_debug.sh; fi
if test -e lib/lib_mshplt_debug.a; then . shell/mklib_mshplt_incl_debug.sh; fi
if test -e lib/lib_user_debug.a; then . shell/mklib_user_incl_debug.sh; fi
if test -e lib/libwave_debug.a; then . shell/mklib_wave_incl_debug.sh; fi
if test -e lib/libwave_omp_debug.a; then . shell/mklib_wave_omp_incl_debug.sh; fi

gfortran -g -cpp \
-fopenmp \
-fcheck=bounds \
-fbacktrace \
-ffpe-summary=invalid,zero,overflow \
-fdec -fd-lines-as-comments \
-Wno-align-commons \
-ffixed-line-length-none \
-finit-local-zero \
-funroll-loops \
-o $WAVE_INCL/bin/wave_debug.exe \
wave_main.f \
$WAVE_INCL/lib/libwave_debug.a \
$WAVE_INCL/lib/libuser_debug.a \
$WAVE_INCL/lib/libuser_modules_debug.a \
$WAVE_INCL/lib/libmhbook_debug.a \
$WAVE_INCL/lib/libmhbook_modules_debug.a \
$WAVE_INCL/lib/libmshplt_debug.a \
$WAVE_INCL/lib/libmshcern_debug.a \
$WAVE_INCL/lib/libwave_modules_debug.a \
$WAVE_INCL/lib/libwave_omp_debug.a \
$WAVE_INCL/lib/libwave_omp_modules_debug.a \
$WAVE_INCL/lib/libwave_debug.a \
$WAVE_INCL/lib/libwave_omp_debug.a \

cd $OLDPWD
