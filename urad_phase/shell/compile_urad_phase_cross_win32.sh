# +PATCH,//WAVE/SHELL
# +DECK,compile_urad_phase_win32  ,T=SHELL.

cd for

echo
echo
echo
echo
echo
echo
echo
echo
echo
echo
echo
echo
echo '---------------------------------------------'
echo
echo 'Compiling mshcern.f'

# mshcern.f is created by cat $WAVE_INCL/mshcern/*.f > mshcern.f
x86_64-w64-mingw32-gfortran-win32 -c -O3 -cpp -w \
-ffpe-summary=invalid,zero,overflow \
-fdec -fd-lines-as-comments \
-Wno-align-commons \
-ffixed-line-length-none \
-finit-local-zero -funroll-loops \
mshcern.f

echo
echo 'Compiling urad_modules.f'
x86_64-w64-mingw32-gfortran-win32 -c -O3 -cpp \
-ffpe-summary=invalid,zero,overflow \
-fopenmp \
-fcheck=all \
-fdec -fd-lines-as-comments \
-Wno-align-commons \
-ffixed-line-length-none \
-finit-local-zero -funroll-loops \
urad_modules.f

x86_64-w64-mingw32-gfortran-win32 -c -O3 -cpp \
-ffpe-summary=invalid,zero,overflow \
-fopenmp \
-fcheck=all \
-fdec -fd-lines-as-comments \
-Wno-align-commons \
-ffixed-line-length-none \
-finit-local-zero -funroll-loops \
urad_util.f
echo
echo 'Compiling urad_util.f'

x86_64-w64-mingw32-gfortran-win32 -O3 -cpp \
-static \
-ffpe-summary=invalid,zero,overflow \
-fcheck=all \
-fopenmp \
-fdec -fd-lines-as-comments \
-Wno-align-commons \
-ffixed-line-length-none \
-finit-local-zero -funroll-loops \
-o ../bin/urad_phase_win32.exe \
urad_phase_main.f \
urad_modules.o urad_util.o \
mshcern.o \

echo
echo 'Compiling and link urad_phase_main.f'
echo
if test -e ../bin/urad_phase_win32.exe; then
  echo 'urad_phase_win32.exe succesfully created'
else
  echo '*** Failed to create urad_phase_win32.exe ***'
fi

cd ..
