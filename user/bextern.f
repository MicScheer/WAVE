*CMZ :  4.01/00 12/03/2023  14.27.05  by  Michael Scheer
*CMZ :  4.00/11 06/07/2021  11.35.01  by  Michael Scheer
*CMZ :  3.03/02 15/12/2015  15.57.23  by  Michael Scheer
*CMZ :  3.01/02 24/01/2014  17.44.47  by  Michael Scheer
*CMZ :  3.01/00 04/07/2013  08.33.03  by  Michael Scheer
*CMZ :  3.00/01 02/04/2013  13.47.23  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.23  by  Michael Scheer
*CMZ :  2.63/05 12/05/2010  13.34.28  by  Michael Scheer
*CMZ :  2.61/06 12/04/2007  09.47.01  by  Michael Scheer
*CMZ :  2.58/01 15/03/2007  11.13.54  by  Michael Scheer
*CMZ :  2.57/05 24/08/2006  16.44.25  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.35.17  by  Michael Scheer
*CMZ : 00.01/02 04/11/94  14.47.19  by  Michael Scheer
*CMZ : 00.00/07 01/06/94  10.41.53  by  Michael Scheer
*-- Author : Michael Scheer
C-----------------------------------------------------------
      subroutine bextern(xin,yin,zin,bxout,byout,bzout,axout,ayout,azout)
C-----------------------------------------------------------

      implicit none

*KEEP,uservar.
      include 'uservar.cmn'
*KEND.

C--- USER ROUTINE TO CALCULATE MAGNETIC FIELD

C     INPUT: X, Y, Z
C     OUTPUT: MAGNETIC FIELD BXOUT,BYOUT,BZOUT
C             VECTOR POTENTIAL AXOUT,AYOUT,AZOUT

      double precision xtab(1000),bytab(1000),bztab(1000)
      integer ntab,luno

      double precision xin,yin,zin,bxout,byout,bzout,axout,ayout,azout



      return
      end
