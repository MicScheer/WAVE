#!/bin/sh

# +PATCH,//WAVE/SHELL
# +DECK,compile_wave_check_system_incl,T=SHELL.

OLDPWD=$PWD

cd $WAVE_INCL/check_system

echo
echo Compiling wave_check_system.f
echo

gfortran -static \
-ffpe-summary=invalid,zero,overflow \
-fdec -fd-lines-as-comments \
-Wno-align-commons -fno-automatic \
-ffixed-line-length-none -finit-local-zero -funroll-loops \
-o $WAVE_INCL/bin/wave_check_system_debug.exe \
*.f

cd $OLDPWD
