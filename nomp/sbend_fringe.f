*CMZ :  4.01/07 24/11/2024  16.01.01  by  Michael Scheer
*-- Author :    Michael Scheer   23/11/2024
      subroutine sbend_fringe(cbmodel,fint,hgap,fringe,fa,fb,fc,istatus)

cdbuff dcbuff
      implicit none

      double precision fint,hgap,gap,fringe,fa,fb,fc,fringe1,fringe2,fringe3,fringe4,fringe5

      integer istatus
      character(32) cbmodel

*KEEP,phyconparam.
      include 'phyconparam.cmn'
*KEND.

      istatus=0
      fa=0.0d0
      fb=0.0d0
      fc=0.0d0
      gap=2.0d0*hgap

      if (cbmodel.eq."hard-edge") then
        fringe=0.0d0
      else if (cbmodel.eq."linear") then
        fringe=6.0d0*fint*gap
        fa=1.0d0/fringe
      else if (cbmodel.eq."cubic-spline") then
        fringe=70.0d0/9.0d0*fint*gap
        fringe2=fringe*fringe
        fringe3=fringe2*fringe
        fb=-2.0d0/fringe3
        fa=3.0d0/fringe2
      else if (cbmodel.eq."quintic-spline") then
        fringe=231.0d0*fint*gap/25.0d0
        fringe2=fringe*fringe
        fringe3=fringe2*fringe
        fringe4=fringe2*fringe2
        fringe5=fringe3*fringe2
        fa=10.0d0/fringe3
        fb=-15.0d0/fringe4
        fc=6.0d0/fringe5
      else
        istatus=-1
        fringe=0.0d0
      endif

      end
