# +PATCH,//WAVE/SHELL
# +DECK,compile_urad_phase_debug  ,T=SHELL.
rm -f b/urad_phase_debug.exe
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

gfortran -c -g -cpp \
-ffpe-summary=invalid,zero,overflow \
-fcheck=all \
-fopenmp \
-fdec -fd-lines-as-comments \
-Wno-align-commons \
-ffixed-line-length-none \
-finit-local-zero -funroll-loops \
urad_modules.f

gfortran -c -g -cpp \
-ffpe-summary=invalid,zero,overflow \
-fopenmp \
-fdec -fd-lines-as-comments \
-Wno-align-commons \
-ffixed-line-length-none \
-finit-local-zero -funroll-loops \
-fcheck=all \
urad_util.f

gfortran -g -cpp \
-ffpe-summary=invalid,zero,overflow \
-fopenmp \
-fdec -fd-lines-as-comments \
-Wno-align-commons \
-ffixed-line-length-none \
-finit-local-zero -funroll-loops \
-fcheck=all \
urad_modules.o urad_util.o \
-o ../bin/urad_phase_debug.exe \
urad_phase_main.f

c ..
