Installation of WAVE:
---------------------

Make WAVE:
----------

Prerequisits: python3, gfortran (if you need to compile WAVE)

It is recommended to define the system Variable WAVE_INCL, i.e.
probably the parent directory of the file. If you have e.g. WAVE in
~/wave, this means

export WAVE_INCL=~/wave/for_incl

To create it from scratch, use the script in $WAVE_INCL/shell or delete
$WAVE_INCL/bin/wave_debug.exe and use

python3 $WAVE_INCL/python/make_wave.py

To run WAVE, refor to wave_by_examples.pdf in $WAVE_INCL/doc.
