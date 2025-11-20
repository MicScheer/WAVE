echo " "
echo "----------------------------------------------------------------"
echo "This is set_wave_environment.sh"
echo "----------------------------------------------------------------"
echo " "

if test x$WAVE = x; then
   echo Shell variable WAVE not defined, trying to set it...

   cd ..

   no_wave='no'
   for d in 'stage' 'bin' 'python' 'shell'; do
     (ls -la | grep  -q $d) || no_wave='yes'
     if test x$no_wave = xyes; then
       echo '*** directory ' $d not found '***' && echo "This script must be run from a sub-directory of WAVE's home"
       return
     fi
   done
   export WAVE=`pwd`
fi

echo "Shell variable WAVE, i.e. WAVE's home is $WAVE"
echo " "

if test x$WAVE_INCL = x; then
  export WAVE_INCL=$WAVE
fi

echo " "
echo "Shell variable WAVE_INCL, i.e. WAVE's source directory is $WAVE_INCL"
echo " "

echo "Checking installation"
echo " "

   no_wave='no'
   for d in 'stage' 'bin' 'python'; do
     (ls -la $WAVE | grep  -q $d) || no_wave='yes'
     if test x$no_wave = xyes; then
       echo '*** directory ' $d not found '***' && echo "Check $WAVE"
       return
     fi
   done

   no_wave='no'
   #for f in 'bin/wave.exe' 'python/waveplot.py' 'python/waves.py' 'python/waveshop.py' 'stage/wave.in' 'stage/wave' 'stage/undumag'; do
   for f in 'bin/wave.exe' 'python/waveplot.py' 'python/waveshop.py' 'stage/wave.in' 'stage/wave' 'stage/undumag'; do
     ls -la $WAVE/$f > /dev/null || no_wave='yes'
     if test x$no_wave = xyes; then
       echo "'*** missing' $f  '***'"
       return
     else
       echo "--- found $f  ---"
     fi
   done

if ! echo $PATH | grep -q $WAVE/stage; then
  export PATH=$WAVE/bin:$WAVE/shell:$WAVE/python:$WAVE/stage:$PATH
  echo " "
  echo "Setting PATH variable:"
  echo $PATH
  echo " "
fi

echo " "
echo "Setting command aliases"
echo " "

  alias es='$EDITOR $WAVE/shell/set_wave_environment.sh'
  alias ss='. $WAVE/shell/set_wave_environment.sh'
  alias cdwave='. $WAVE/shell/set_wave_environment.sh'

  alias wplot='cd $WAVE/stage; ipython3 -i $WAVE/python/waveplot.py'
  alias waveplot='cd $WAVE/stage; ipython3 -i $WAVE/python/waveplot.py'

  unalias r 2>/dev/null

  export UNDUMAG=$WAVE/undumag
  echo
  echo Environment variable UNDUMAG set to: $UNDUMAG
  echo
  alias wave='. $WAVE/stage/wave'
  alias undumag='. $WAVE/stage/undumag'
  #alias waves='cd $WAVE/stage; python3 -i $WAVE/python/waves.py'
  alias waveshop='cd $WAVE/stage; ipython3 -i $WAVE/python/waveshop.py'
  alias wave
  alias undumag
  #alias waves
  alias wplot
  echo

  cd $WAVE/stage

echo " "
echo "To run WAVE, enter wave"
echo "For plotting results, enter wplot"
echo "To start GUI, try waveshop or waves"
echo " "
