#!/bin/sh

# +PATCH,//WAVE/SHELL
# +DECK,mklib_wave_incl_debug,T=SHELL.

OLDPWD=$PWD

cd $WAVE_INCL/nomp

echo
echo
echo Compiling modules in $WAVE_INCL/nomp/mod
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

echo
echo Making $WAVE_INCL/lib/libwave_modules_debug.a
echo

ar rc $WAVE_INCL/lib/libwave_modules_debug.a *.o
ranlib $WAVE_INCL/lib/libwave_modules_debug.a

cd ..

echo
echo Compiling "$WAVE_INCL/nomp/*.f"
echo

gfortran -c -g -cpp \
-fcheck=bounds \
-fbacktrace \
-ffpe-summary=invalid,zero,overflow \
-fdec -fd-lines-as-comments \
-Wno-align-commons -fno-automatic -ffixed-line-length-none \
-finit-local-zero \
-funroll-loops \
*.f

echo
echo Making $WAVE_INCL/lib/libwave_debug.a

ar rc $WAVE_INCL/lib/libwave_debug.a *.o
ranlib $WAVE_INCL/lib/libwave_debug.a

cd $OLDPWD
