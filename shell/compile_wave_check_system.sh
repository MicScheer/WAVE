#!/bin/sh

# +PATCH,//WAVE/SHELL
# +DECK,compile_wave_check_system,T=SHELL.

OLDPWD=$PWD

cd $WAVE/for

echo
echo Compiling wave_check_system.f
echo

gfortran -static \
-ffpe-summary=invalid,zero,overflow \
-fdec -fd-lines-as-comments \
-Wno-align-commons -fno-automatic \
-ffixed-line-length-none -finit-local-zero -funroll-loops \
-o $WAVE/bin/wave_check_system_debug.exe \
wave_check_system.f

cd $OLDPWD
