#!/bin/sh

# +PATCH,//WAVE/SHELL
# +DECK,mklib_wave_omp_incl_debug,T=SHELL.

OLDPWD=$PWD

cd $WAVE_INCL/omp

echo
echo Compiling "$WAVE_INCL/omp/mod/*.f"

rm -f *.mod
cd mod

gfortran -c -g -cpp \
-finit-local-zero \
-fcheck=bounds \
-fopenmp \
-fbacktrace \
-ffpe-summary=invalid,zero,overflow \
-fdec -fd-lines-as-comments \
-Wno-align-commons \
-ffixed-line-length-none -funroll-loops \
*.f

mv *.mod ..

echo
echo Making $WAVE_INCL/lib/libwave_omp_modules_debug.a

ar rc $WAVE_INCL/lib/libwave_omp_modules_debug.a *.o
ranlib $WAVE_INCL/lib/libwave_omp_modules_debug.a

echo
echo Compiling "$WAVE_INCL/omp/mod/*.f"

cd ..

gfortran -c -g -cpp \
-finit-local-zero \
-fcheck=bounds \
-fopenmp \
-fbacktrace \
-ffpe-summary=invalid,zero,overflow \
-fdec -fd-lines-as-comments \
-Wno-align-commons \
-ffixed-line-length-none -funroll-loops \
*.f

echo
echo Making $WAVE_INCL/lib/libwave_omp_debug.a

ar rc $WAVE_INCL/lib/libwave_omp_debug.a *.o
ranlib $WAVE_INCL/lib/libwave_omp_debug.a


cd $OLDPWD
