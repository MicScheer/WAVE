@echo off

rem +PATCH,//WAVE/BAT
rem +DECK,set_wave_environment,T=BAT.

if defined WAVE (
  cd %WAVE%\stage
) else (

  set WAVE=%CD%
  cd bat
  extend_path.bat %cd%\bat
  extend_path.bat %cd%\bin
  cd ..\stage
)

