#!/bin/sh

# +PATCH,//WAVE/SHELL
# +DECK,mklib_user_incl_debug,T=SHELL.

OLDPWD=$PWD

cd $WAVE_INCL/user

rm -f *.mod 2>/dev/null
cd mod

echo

if test -e *.f; then

  echo
  echo Compiling "$WAVE_INCL/user/mod/*.f"
  echo

  gfortran -c -w -g -cpp \
    -finit-local-zero \
    -fcheck=bounds \
    -fno-automatic \
    -ffpe-summary=invalid,zero,overflow \
    -fdec -fd-lines-as-comments \
    -Wno-align-commons \
    -ffixed-line-length-none -funroll-loops \
  *.f

  mv *.mod ..

  echo
  echo Making $WAVE_INCL/lib/libuser_debug.a
  echo

  ar rc $WAVE_INCL/lib/libuser_modules_debug.a *.o
  ranlib $WAVE_INCL/lib/libuser_modules_debug.a

fi

echo
echo Compiling "$WAVE_INCL/user/*.f"
echo

cd ..

gfortran -c -w -g -cpp \
-finit-local-zero \
-fcheck=bounds \
-fno-automatic \
-ffpe-summary=invalid,zero,overflow \
-fdec -fd-lines-as-comments \
-Wno-align-commons \
-ffixed-line-length-none -funroll-loops \
*.f

echo
echo Making $WAVE_INCL/lib/libuser_debug.a
echo

ar rc $WAVE_INCL/lib/libuser_debug.a *.o
ranlib $WAVE_INCL/lib/libuser_debug.a

cd $OLDPWD
