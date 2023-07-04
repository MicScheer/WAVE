#!/bin/sh

# +PATCH,//WAVE/SHELL
# +DECK,mklib_mhbook_incl_debug,T=SHELL.

OLDPWD=$PWD

cd $WAVE_INCL/mhbook

echo
echo
echo Compiling $WAVE_INCL/mhbook/mod/mhbook_module.f
echo

rm -f *.mod

cd mod

gfortran -c -w -g -cpp \
-ffpe-summary=invalid,zero,overflow \
-fdec -fd-lines-as-comments \
-Wno-align-commons -fno-automatic -ffixed-line-length-none \
-finit-local-zero \
-funroll-loops \
*.f

mv *.mod ..

echo
echo Making $WAVE_INCL/lib/libmhbook_modules_debug.a
echo

ar rc $WAVE_INCL/lib/libmhbook_modules_debug.a *.o
ranlib $WAVE_INCL/lib/libmhbook_modules_debug.a

cd ..

echo
echo Compiling "$WAVE_INCL/mhbook/*.f"
echo

gfortran -c -w -g -cpp \
-ffpe-summary=invalid,zero,overflow \
-fdec -fd-lines-as-comments \
-Wno-align-commons -fno-automatic -ffixed-line-length-none \
-finit-local-zero \
-funroll-loops \
*.f

echo
echo Making $WAVE_INCL/lib/libmhbook_debug.a
echo

ar rc $WAVE_INCL/lib/libmhbook_debug.a *.o
ranlib $WAVE_INCL/lib/libmhbook_debug.a

cd $OLDPWD
