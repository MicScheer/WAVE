Installation of WAVE:
---------------------

Make WAVE:
----------

Prerequisits: python3, gfortran (if you need to compile WAVE)

It is recommended to define the system Variable WAVE_INCL, i.e.
probably the parent directory of the file. If you have e.g. WAVE in
~/wave, this means

export WAVE_INCL=~/wave

To create it from scratch, delete all libraries in $WAVE_INCL/lib and
enter:

python3 $WAVE_INCL/python/make_wave.py

To run WAVE

export WAVE=~/wave
. $WAVE/shell/set_wave_environment.sh

