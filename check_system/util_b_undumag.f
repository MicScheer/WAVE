*CMZ : 00.00/16 13/07/2015  11.27.12  by  Michael Scheer
*CMZ : 00.00/11 23/03/2011  16.28.17  by  Michael Scheer
*-- Author : Michael Scheer
      subroutine util_b_elleaume(perlen,gap,b,isilent)

      implicit none

      double precision perlen,gap,b,dk
      integer isilent

C Calulates magnetic field By for Hybrid-Undulator according to
c http://cas.web.cern.ch/cas/Belgium-2009/Lectures/PDFs/Bahrdt-3.pdf
c for M=4

      b=3.69d0*exp(-5.07d0*gap/perlen+1.52d0*(gap/perlen)**2)
      dk=0.934*b*perlen/10.

      if (isilent.eq.0) then
        print*,'B:',b
        print*,'K (for perlen in mm):',dk
      endif

      return
      end
