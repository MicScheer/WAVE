#!/bin/bash

# +PATCH,//WAVE/SHELL
# +DECK,mklib_wave_incl_all_debug,T=SHELL.

if test $WAVE_INCL; then

  cd $WAVE_INCL

  rm -f lib/*

  . $WAVE_INCL/shell/mklib_mshplt_incl_debug.sh
  . $WAVE_INCL/shell/mklib_mshcern_incl_debug.sh
  . $WAVE_INCL/shell/mklib_mhbook_incl_debug.sh
  . $WAVE_INCL/shell/mklib_user_incl_debug.sh
  . $WAVE_INCL/shell/mklib_wave_omp_incl_debug.sh
  . $WAVE_INCL/shell/mklib_urad_incl_debug.sh
  . $WAVE_INCL/shell/mklib_wave_incl_debug.sh

  cd $WAVE_INCL

fi
