*CMZ :  4.01/00 12/03/2023  14.28.40  by  Michael Scheer
*CMZ :  3.03/02 15/12/2015  16.08.16  by  Michael Scheer
*CMZ :  3.03/01 07/10/2015  14.45.17  by  Michael Scheer
*CMZ :  3.02/04 12/12/2014  15.55.26  by  Michael Scheer
*CMZ :  3.00/00 14/03/2013  11.59.45  by  Michael Scheer
*CMZ :  2.66/20 06/06/2012  15.22.24  by  Michael Scheer
*CMZ :  2.66/13 08/12/2010  09.44.42  by  Michael Scheer
*CMZ :  2.63/03 12/05/2010  13.34.28  by  Michael Scheer
*CMZ :  2.61/02 20/03/2007  16.42.23  by  Michael Scheer
*-- Author : Michael Scheer
      subroutine ustep(x,y,z,vx,vy,vz,vxp,vyp,vzp,dt,gamma,dgamma)

c user routine, called after each tracking step

      implicit none

*KEEP,ustep.
      include 'ustep.cmn'
*KEEP,uservar.
      include 'uservar.cmn'
*KEND.

      double precision x,y,z,vx,vy,vz,vxp,vyp,vzp,dt,gamma,dgamma

c     x,y,z is the current position of the electron in meter
c     vx,vy,vz is the current velocity of the electron in meter/sec
c     vxp,vyp,vzp the derivative of the velocity in meter/sec**2
c     dt the time step in sec
c     gamma the current gamma factor
c     dgamma the energy loss of step in terms of gamma

      return
      end
