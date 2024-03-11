# +PATCH,//WAVE/SHELL
# +DECK,set_wave_environment,T=SHELL.

echo
if test $SHELL = /bin/bash; then
  cd ${BASH_ARGV[0]%/*}
else
  cd ${0%/*}
fi

cd ..
export WAVE=${PWD}
echo Environment WAVE set to: $WAVE

if test -d $WAVE; then

  cd $WAVE

  alias es='$EDITOR $WAVE/shell/set_wave_environment.sh'
  alias ss='. $WAVE/shell/set_wave_environment.sh'
  alias cdwave='. $WAVE/shell/set_wave_environment.sh'

  alias wplot='cd $WAVE/stage; ipython3 -i $WAVE/python/waveplot.py'

  unalias r 2>/dev/null

  export UNDUMAG=$WAVE/undumag
  echo
  echo Environment variable UNDUMAG set to: $UNDUMAG
  echo
  alias wave='. $WAVE/stage/wave'
  alias undumag='. $WAVE/stage/undumag'
  alias waves='cd $WAVE/stage; bash $WAVE/shell/title.sh; ipython3 -i $WAVE/python/waves.py'
  alias wave
  alias undumag
  alias waves
  alias wplot
  echo

  cd $WAVE/stage

else
  echo
  echo Error: $WAVE not found or it is not a directory!
  echo
fi
