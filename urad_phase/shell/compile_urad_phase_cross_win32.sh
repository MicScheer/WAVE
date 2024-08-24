# +PATCH,//WAVE/SHELL
# +DECK,compile_urad_phase  ,T=SHELL.
rm -f b/urad_phase_win32.exe

c f

echo
echo
echo

x86_64-w64-mingw32-gfortran-win32 -c -O3 -cpp \
-ffpe-summary=invalid,zero,overflow \
-fopenmp \
-fdec -fd-lines-as-comments \
-Wno-align-commons \
-ffixed-line-length-none \
-finit-local-zero -funroll-loops \
urad_modules.f

x86_64-w64-mingw32-gfortran-win32 -c -O3 -cpp \
-ffpe-summary=invalid,zero,overflow \
-fopenmp \
-fdec -fd-lines-as-comments \
-Wno-align-commons \
-ffixed-line-length-none \
-finit-local-zero -funroll-loops \
urad_util.f

x86_64-w64-mingw32-gfortran-win32 -static -O3 -cpp \
-ffpe-summary=invalid,zero,overflow \
-fopenmp \
-fdec -fd-lines-as-comments \
-Wno-align-commons \
-ffixed-line-length-none \
-finit-local-zero -funroll-loops \
urad_modules.o urad_util.o \
-o ../bin/urad_phase_win32.exe \
urad_phase_main.f

c ..
