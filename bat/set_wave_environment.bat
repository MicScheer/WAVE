@echo off

rem +PATCH,//WAVE/BAT
rem +DECK,set_wave_environment,T=BAT.

rem echo %0


if defined WAVE (cd %WAVE%\stage) else (

  set WAVE=%CD%

  cd bat

  extend_path.bat %WAVE%\bat
  extend_path.bat %WAVE%\bin

  cd ..\stage
)

