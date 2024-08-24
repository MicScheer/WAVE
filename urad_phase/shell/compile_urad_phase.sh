# +PATCH,//WAVE/SHELL
# +DECK,compile_urad_phase  ,T=SHELL.
rm -f b/urad_phase.exe
rm -f f/*.o

c f

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

gfortran -c -O3 -cpp \
-ffpe-summary=invalid,zero,overflow \
-fopenmp \
-fcheck=all \
-fdec -fd-lines-as-comments \
-Wno-align-commons \
-ffixed-line-length-none \
-finit-local-zero -funroll-loops \
urad_modules.f

gfortran -c -O3 -cpp \
-ffpe-summary=invalid,zero,overflow \
-fopenmp \
-fcheck=all \
-fdec -fd-lines-as-comments \
-Wno-align-commons \
-ffixed-line-length-none \
-finit-local-zero -funroll-loops \
urad_util.f

gfortran -O3 -cpp \
-ffpe-summary=invalid,zero,overflow \
-fcheck=all \
-fopenmp \
-fdec -fd-lines-as-comments \
-Wno-align-commons \
-ffixed-line-length-none \
-finit-local-zero -funroll-loops \
urad_modules.o urad_util.o \
-o ../bin/urad_phase.exe \
urad_phase_main.f

c ..
