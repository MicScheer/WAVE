# +PATCH,//WAVE/SHELL
# +DECK,set_wave_environment,T=SHELL.

if test -z $WAVE; then
  #echo
  #echo Environment variable WAVE is not defined. setting it to:
  cd ${0%/*}
  export WAVE=${PWD%/*}
  #echo $WAVE
  #echo
fi

if test -d $WAVE; then

  cd $WAVE

  alias es='$EDITOR $WAVE/shell/set_wave_environment.sh'
  alias ss='. $WAVE/shell/set_wave_environment.sh'
  alias cdwave='. $WAVE/shell/set_wave_environment.sh'

  alias wplot='cd $WAVE/stage; ipython3 -i $WAVE/python/waveplot.py'

  alias wave='$WAVE/stage/wave'
  alias undumag='$WAVE/stage/undumag'
  alias waves='cd $WAVE/stage; ipython3 -i $WAVE/python/waves.py'
  alias waveshop='cd $WAVE/stage; ipython3 -i $WAVE/python/waveshop.py'

  cd $WAVE/stage

else
  echo
  echo Error: $WAVE not found or it is not a directory!
  echo
fi
