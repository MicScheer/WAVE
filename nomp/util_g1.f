*CMZ :  4.00/04 17/05/2019  14.17.20  by  Michael Scheer
*CMZ :  3.03/02 19/11/2015  13.51.04  by  Michael Scheer
*CMZ :  3.02/05 27/03/2015  15.15.42  by  Michael Scheer
*CMZ :  2.70/05 02/01/2013  12.47.25  by  Michael Scheer
*CMZ :  2.67/04 15/05/2012  10.53.28  by  Michael Scheer
*-- Author :    Michael Scheer   10/05/2012
      subroutine util_g1(y,g1)
*KEEP,gplhint.
!******************************************************************************
!
!      Copyright 2013 Helmholtz-Zentrum Berlin (HZB)
!      Hahn-Meitner-Platz 1
!      D-14109 Berlin
!      Germany
!
!      Author Michael Scheer, Michael.Scheer@Helmholtz-Berlin.de
!
! -----------------------------------------------------------------------
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy (wave_gpl.txt) of the GNU General Public
!    License along with this program.
!    If not, see <http://www.gnu.org/licenses/>.
!
!    Dieses Programm ist Freie Software: Sie koennen es unter den Bedingungen
!    der GNU General Public License, wie von der Free Software Foundation,
!    Version 3 der Lizenz oder (nach Ihrer Option) jeder spaeteren
!    veroeffentlichten Version, weiterverbreiten und/oder modifizieren.
!
!    Dieses Programm wird in der Hoffnung, dass es nuetzlich sein wird, aber
!    OHNE JEDE GEWAEHRLEISTUNG, bereitgestellt; sogar ohne die implizite
!    Gewaehrleistung der MARKTFAEHIGKEIT oder EIGNUNG FueR EINEN BESTIMMTEN ZWECK.
!    Siehe die GNU General Public License fuer weitere Details.
!
!    Sie sollten eine Kopie (wave_gpl.txt) der GNU General Public License
!    zusammen mit diesem Programm erhalten haben. Wenn nicht,
!    siehe <http://www.gnu.org/licenses/>.
!
!******************************************************************************
*KEND.

c calculates G1(y) with an estimated precision of about 1.5e-3 for y<=30,
c and 1.5e-2 for y>30.


      implicit none

      integer ical,npoi,ieof,ipoi,npoilow,npoihigh

      double precision y,g1,g1_30,c_30,g1_5em5,c_5em5,y_30,y_5em5,ydum,g1dum,r1

      double precision, dimension(:), allocatable ::
     &  ywlow,ywhigh,coeflow,coefhigh,
     &  g1wlow,g1whigh,g1walow,g1wahigh,r1low,r1high,
     &  w1,w2,w3,w4

      save

      data ical/0/
      data y_30/30.0d0/
      data y_5em5/5.0d-5/
      data g1_30/6.580794488121591d-013/ !WAVE
      data g1_5em5/7.909860755922665E-002/ !WAVE

      if (ical.eq.0) then

        c_5em5=g1_5em5/y_5em5**(1.0d0/3.0d0)
        c_30=g1_30/(sqrt(y_30)*exp(-y_30))

        npoi=0
        npoilow=0
        npoihigh=0

        open(unit=99,file='wave-g1.dat',status='old',err=9999)

1       call util_skip_comment_end(99,ieof)
        if (ieof.ne.0) goto 9
        read(99,*)ydum
        npoi=npoi+1
        if (ydum.ge.y_5em5.and.ydum.le.4.0d0) then !zwei Abfragen wegen 4.0
          npoilow=npoilow+1
        endif
        if (ydum.ge.4.0d0.and.ydum.le.y_30) then
          npoihigh=npoihigh+1
        endif
        goto 1

9       rewind(99)

        allocate(ywlow(npoilow))
        allocate(g1wlow(npoilow))
        allocate(g1walow(npoilow))
        allocate(r1low(npoilow))
        allocate(ywhigh(npoihigh))
        allocate(g1whigh(npoihigh))
        allocate(g1wahigh(npoihigh))
        allocate(r1high(npoihigh))
        allocate(coeflow(npoilow))
        allocate(coefhigh(npoihigh))
        allocate(w1(max(npoilow,npoihigh)))
        allocate(w2(max(npoilow,npoihigh)))
        allocate(w3(max(npoilow,npoihigh)))
        allocate(w4(max(npoilow,npoihigh)))

        npoilow=0
        npoihigh=0
        do ipoi=1,npoi
          call util_skip_comment_end(99,ieof)
          read(99,*)ydum,g1dum
          if (ydum.ge.y_5em5.and.ydum.le.4.0d0) then
            npoilow=npoilow+1
            ywlow(npoilow)=ydum
            g1wlow(npoilow)=g1dum
            g1walow(npoilow)=
     &        391.8d0 * ydum**(1.0d0/3.0d0) * exp(-ydum*0.8307d0)
     &        -192.0d0 * sqrt(ydum) * exp(-ydum*0.7880d0)
          endif
          if (ydum.ge.4.0d0.and.ydum.le.y_30) then
            npoihigh=npoihigh+1
            ywhigh(npoihigh)=ydum
            g1whigh(npoihigh)=g1dum
            g1wahigh(npoihigh)=164.0d0*sqrt(ydum)* EXP(-ydum)
          endif
        enddo
        close(99)

        r1low(1:npoilow)=g1wlow(1:npoilow)/g1walow(1:npoilow)
        r1high(1:npoihigh)=g1whigh(1:npoihigh)/g1wahigh(1:npoihigh)

        call util_spline_coef(ywlow,r1low,npoilow,0.0d0,0.0d0,coeflow,
     &    w1,w2,w3,w4)
        call util_spline_coef(ywhigh,r1high,npoihigh,0.0d0,0.0d0,
     &    coefhigh,w1,w2,w3,w4)

        ical=1
      endif

      if (y.le.5.0d-5) then
        g1=c_5em5*y**(1.0d0/3.0d0)
      else if (y.ge.30.0d0) then
        g1=c_30*sqrt(y)*exp(-y)
      else

        if (y.ge.y_5em5.and.y.lt.4.0d0) then
          call util_spline_inter(ywlow,r1low,coeflow,npoilow,y,r1,-1)
          g1=r1*(
     &      391.8d0 * y**(1.0d0/3.0d0) * exp(-y*0.8307d0)
     &      -192.0d0 * sqrt(y) * exp(-y*0.7880d0))
        else if (y.ge.4.0d0.and.y.le.y_30) then
          call util_spline_inter(
     &      ywhigh,r1high,coefhigh,npoihigh,y,r1,-1)
          g1=r1*(164.0d0*sqrt(y)* EXP(-y))
        endif

      endif

      return
9999  stop '*** File wave-g1.dat not found ***'
      end
