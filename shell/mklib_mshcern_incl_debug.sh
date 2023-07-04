#!/bin/sh

# +PATCH,//WAVE/SHELL
# +DECK,mklib_mshcern_incl_debug,T=SHELL.

OLDPWD=$PWD

cd $WAVE_INCL/mshcern

echo
echo
echo Compiling "$WAVE_INCL/mshcern/*.f"
echo

gfortran -c -w -g -cpp \
-ffpe-summary=invalid,zero,overflow \
-fdec -fd-lines-as-comments \
-Wno-align-commons -fno-automatic -ffixed-line-length-none \
-finit-local-zero \
-funroll-loops \
*.f

echo
echo Making $WAVE_INCL/lib/libmshcern_debug.a
echo

ar rc $WAVE_INCL/lib/libmshcern_debug.a *.o
ranlib $WAVE_INCL/lib/libmshcern_debug.a

cd $OLDPWD
