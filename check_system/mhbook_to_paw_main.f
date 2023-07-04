*CMZ :  3.02/04 11/11/2014  16.45.38  by  Michael Scheer
*CMZ :  3.02/00 13/10/2014  09.15.08  by  Michael Scheer
*-- Author :    Michael Scheer   10/10/2014
      program mhbook_to_paw_main
c+seq,gplhint.

c--------------------------------------------------------------
c      Converts WAVE histogram file to HBOOK/PAW format of CERN
c      The name of the output file is converted to lower case!

      implicit none

      character(2048) filein,fileout,command
      integer narg,last,istat

      call get_command(command,last,istat)
      narg=command_argument_count()

      if (narg.eq.0) then
        filein='WAVE.mhb'
        fileout='wave_histo.his'
      else if (narg.eq.1) then
        call getarg(1,filein)
        fileout='wave_histo.his'
      else if (narg.eq.2) then
        call getarg(1,filein)
        call getarg(2,fileout)
      else
        stop '*** Expected input file [WAVE.mnb] and output file [wave_histo.his]!'
      endif

      call util_lower_case(fileout)
      call mhbook_to_paw(filein,fileout)

      stop
      end
      include 'mhbook_to_paw.f'
