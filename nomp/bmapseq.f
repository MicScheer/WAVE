*CMZ :          21/01/2025  16.53.27  by  Michael Scheer
*CMZ :  4.00/16 09/08/2022  09.07.08  by  Michael Scheer
*CMZ :  4.00/07 07/06/2020  15.15.28  by  Michael Scheer
*CMZ :  3.05/05 13/07/2018  11.51.31  by  Michael Scheer
*CMZ :  3.03/04 29/11/2017  10.21.49  by  Michael Scheer
*CMZ :  3.02/00 10/09/2014  14.10.09  by  Michael Scheer
*CMZ :  3.01/04 26/05/2014  16.15.21  by  Michael Scheer
*CMZ :  3.01/00 06/05/2013  09.13.42  by  Michael Scheer
*CMZ :  2.68/05 01/10/2012  14.11.38  by  Michael Scheer
*CMZ :  2.68/04 04/09/2012  09.22.40  by  Michael Scheer
*CMZ :  2.68/03 31/08/2012  09.01.42  by  Michael Scheer
*CMZ :  2.68/02 02/07/2012  13.51.47  by  Michael Scheer
*-- Author :    Michael Scheer   18/06/2012
C**********************************************************************
      subroutine bmapseq(imag,xin,yin,zin,bxout,byout,bzout)
C**********************************************************************

      use magseqf90m

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

      implicit none

*KEEP,contrl.
      include 'contrl.cmn'
*KEND.

      double precision xin,yin,zin,x,y,z,bx,by,bz,step,stepy,xold,yold,
     &  bxout,byout,bzout,xmin,xmax,ymin,ymax,zmin,zmax,
     &  x1,x2,a3(3),x3(3),b3(3),
     &  b11,b12,b13,b21,b22,b23,b31,b32,b33,bb1,bb2,bb3,dy,dz

      double precision b111(3),b211(3), b121(3),b221(3)
      double precision b112(3),b212(3), b122(3),b222(3)
      double precision b112111(3),b212211(3), b122121(3),b222221(3)
      double precision blow(3),bhig(3),b(3),dxx,dyy,dzz,xx,yy,zz,y1,z1

      real*8 :: eps=1.0d-6

      integer imag,last,ianf,
     &  kd,ix1,ix2,iy1,iy2,iz1,iz2,nyz,i,k,
     &  kx1,kx2,kx3,ky1,ky2,ky3,kz1,kz2,kz3,lunb0

      character(2048) cline

      seqmag(imag)%iwarnbmap=0

      if (seqmag(imag)%ical.eq.0) then

        write(lungfo,*)
        write(lungfo,*)'      Subroutine BMAPSEQ:'
        write(lungfo,*)
        write(lungfo,*)'      Reading file'
        write(lungfo,*)'      ',seqmag(imag)%cfilemap
        write(lungfo,*)

        seqmag(imag)%bmxmin= 1.0d99
        seqmag(imag)%bmxmax=-1.0d99
        seqmag(imag)%bmymin= 1.0d99
        seqmag(imag)%bmymax=-1.0d99
        seqmag(imag)%bmzmin= 1.0d99
        seqmag(imag)%bmzmax=-1.0d99

        seqmag(imag)%bmbxmin= 1.0d99
        seqmag(imag)%bmbxmax=-1.0d99
        seqmag(imag)%bmbymin= 1.0d99
        seqmag(imag)%bmbymax=-1.0d99
        seqmag(imag)%bmbzmin= 1.0d99
        seqmag(imag)%bmbzmax=-1.0d99

        seqmag(imag)%offsetx=0.0d0
        seqmag(imag)%offsety=0.0d0
        seqmag(imag)%offsetz=0.0d0

        open(unit=lunb0,file=seqmag(imag)%cfilemap,status='old')
        seqmag(imag)%ntot=0
        seqmag(imag)%nx=-1
        seqmag(imag)%ny=-1
        seqmag(imag)%nz=0
 1      read(lunb0,'(a)',end=9) cline
        last=len_trim(cline)

        if (last.le.1.or.
     &      cline(1:1).eq.'%'.or.
     &      cline(1:1).eq.'!'.or.
     &      cline(1:1).eq.'#'.or.
     &      cline(1:1).eq.'*'.or.
     &      cline(1:1).eq.'@'.or.
     &      cline(1:2).eq.' %'.or.
     &      cline(1:2).eq.' !'.or.
     &      cline(1:2).eq.' #'.or.
     &      cline(1:2).eq.' *'.or.
     &      cline(1:2).eq.' @'
     &      ) then

          write(lungfo,*) cline(1:last)

          ianf=index(cline,"scaling")
          if (ianf.gt.0) then
            ianf=index(cline,"=")
            read(cline(ianf+1:last),*)seqmag(imag)%scalex,seqmag(imag)%scaley,seqmag(imag)%scalez,
     &        seqmag(imag)%scalebx,seqmag(imag)%scaleby,seqmag(imag)%scalebz
          endif

          ianf=index(cline,"seqmag(imag)%offset")
          if (ianf.gt.0) then
            ianf=index(cline,"=")
            read(cline(ianf+1:last),*)seqmag(imag)%offsetx,seqmag(imag)%offsety,seqmag(imag)%offsetz,
     &        seqmag(imag)%offsetbx,seqmag(imag)%offsetby,seqmag(imag)%offsetbz
          endif

        else

          seqmag(imag)%ntot=seqmag(imag)%ntot+1
          read(cline(1:last),*)x,y,z,bx,by,bz

          x=x*seqmag(imag)%scalex+seqmag(imag)%offsetx
          y=y*seqmag(imag)%scaley+seqmag(imag)%offsety
          z=z*seqmag(imag)%scalez+seqmag(imag)%offsetz
          bx=bx*seqmag(imag)%scalebx+seqmag(imag)%offsetbx
          by=by*seqmag(imag)%scaleby+seqmag(imag)%offsetby
          bz=bz*seqmag(imag)%scalebz+seqmag(imag)%offsetbz

          if (x.lt.seqmag(imag)%bmxmin) seqmag(imag)%bmxmin=x
          if (x.gt.seqmag(imag)%bmxmax) seqmag(imag)%bmxmax=x
          if (y.lt.seqmag(imag)%bmymin) seqmag(imag)%bmymin=y
          if (y.gt.seqmag(imag)%bmymax) seqmag(imag)%bmymax=y
          if (z.lt.seqmag(imag)%bmzmin) seqmag(imag)%bmzmin=z
          if (z.gt.seqmag(imag)%bmzmax) seqmag(imag)%bmzmax=z
          if (bx.lt.seqmag(imag)%bmbxmin) seqmag(imag)%bmbxmin=bx
          if (bx.gt.seqmag(imag)%bmbxmax) seqmag(imag)%bmbxmax=bx
          if (by.lt.seqmag(imag)%bmbymin) seqmag(imag)%bmbymin=by
          if (by.gt.seqmag(imag)%bmbymax) seqmag(imag)%bmbymax=by
          if (bz.lt.seqmag(imag)%bmbzmin) seqmag(imag)%bmbzmin=bz
          if (bz.gt.seqmag(imag)%bmbzmax) seqmag(imag)%bmbzmax=bz

          if (seqmag(imag)%ntot.eq.1) then
            xold=x
            yold=y
          endif

          if (seqmag(imag)%ntot.eq.2.and.abs(x-xold).gt.eps) then
            write(6,*)'*** Error in BMAPSEQ: Bad file format! ***'
            write(6,*)'*** x must run latest! ***'
            write(lungfo,*)'*** Error in BMAPSEQ: Bad file format! ***'
            write(lungfo,*)'*** x must run latest! ***'
            stop '*** WAVE aborted ***'
          endif

          if (seqmag(imag)%nx.eq.-1) then
            if (seqmag(imag)%ny.eq.-1) then
              if (abs(y-yold).le.eps) then
                seqmag(imag)%nz=seqmag(imag)%nz+1
              else
                seqmag(imag)%ny=seqmag(imag)%nz+1
              endif
            else if (abs(x-xold).le.eps) then
              seqmag(imag)%ny=seqmag(imag)%ny+1
            else
              seqmag(imag)%ny=seqmag(imag)%ny/seqmag(imag)%nz
              seqmag(imag)%nx=0
            endif !seqmag(imag)%ny
          endif !seqmag(imag)%nx

        endif !line type

        goto 1
 9      rewind(lunb0)

        seqmag(imag)%nx=seqmag(imag)%ntot/(seqmag(imag)%ny*seqmag(imag)%nz)

        if (seqmag(imag)%nx.lt.2.and.irfilb0.eq.6.or.seqmag(imag)%nx.lt.3.and.irfilb0.eq.-6) then
          write(lungfo,*)'*** Error in BMAPSEQ: Too few data for field map on'
          write(lungfo,*)seqmag(imag)%cfilemap
          write(lungfo,*)'*** Program WAVE aborted ***'
          write(6,*)'*** Error in BMAPSEQ: Too few data for field map on'
          write(6,*)seqmag(imag)%cfilemap
          write(6,*)'*** Program WAVE aborted ***'
          stop
        endif

        allocate(seqmag(imag)%bmap(6,seqmag(imag)%ntot))

        seqmag(imag)%ntot=0
 11     read(lunb0,'(a)',end=99) cline

        last=len_trim(cline)

        if (last.le.1.or.
     &      cline(1:1).eq.'%'.or.
     &      cline(1:1).eq.'!'.or.
     &      cline(1:1).eq.'#'.or.
     &      cline(1:1).eq.'*'.or.
     &      cline(1:1).eq.'@'.or.
     &      cline(1:2).eq.' %'.or.
     &      cline(1:2).eq.' !'.or.
     &      cline(1:2).eq.' #'.or.
     &      cline(1:2).eq.' *'.or.
     &      cline(1:2).eq.' @'
     &      ) then
        else
          seqmag(imag)%ntot=seqmag(imag)%ntot+1
          read(cline(1:last),*) x,y,z,bx,by,bz
          x=x*seqmag(imag)%scalex+seqmag(imag)%offsetx
          y=y*seqmag(imag)%scaley+seqmag(imag)%offsety
          z=z*seqmag(imag)%scalez+seqmag(imag)%offsetz
          bx=bx*seqmag(imag)%scalebx+seqmag(imag)%offsetbx
          by=by*seqmag(imag)%scaleby+seqmag(imag)%offsetby
          bz=bz*seqmag(imag)%scalebz+seqmag(imag)%offsetbz
          seqmag(imag)%bmap(1,seqmag(imag)%ntot)=x
          seqmag(imag)%bmap(2,seqmag(imag)%ntot)=y
          seqmag(imag)%bmap(3,seqmag(imag)%ntot)=z
          seqmag(imag)%bmap(4,seqmag(imag)%ntot)=bx
          seqmag(imag)%bmap(5,seqmag(imag)%ntot)=by
          seqmag(imag)%bmap(6,seqmag(imag)%ntot)=bz
        endif
        goto 11
 99     close(lunb0)

        step=1.0d0/myinum
        stepy=step/10.0d0
        seqmag(imag)%bmapdy=1.0d0
        if (seqmag(imag)%ny.gt.1) seqmag(imag)%bmapdy=(seqmag(imag)%bmymax-seqmag(imag)%bmymin)/(seqmag(imag)%ny-1)
        seqmag(imag)%bmapdz=1.0d0
        if (seqmag(imag)%nz.gt.1) seqmag(imag)%bmapdz=(seqmag(imag)%bmzmax-seqmag(imag)%bmzmin)/(seqmag(imag)%nz-1)

        ix1=1
        ix2=seqmag(imag)%ntot
        seqmag(imag)%nyz=seqmag(imag)%ny*seqmag(imag)%nz

        xmin= seqmag(imag)%bmxmin+seqmag(imag)%xcen
        xmax= seqmag(imag)%bmxmax+seqmag(imag)%xcen
        ymin= seqmag(imag)%bmymin+seqmag(imag)%ycen
        ymax= seqmag(imag)%bmymax+seqmag(imag)%ycen
        zmin= seqmag(imag)%bmzmin+seqmag(imag)%zcen
        zmax= seqmag(imag)%bmzmax+seqmag(imag)%zcen

        seqmag(imag)%xmin=xmin
        seqmag(imag)%xmax=xmax
        seqmag(imag)%ymin=ymin
        seqmag(imag)%ymax=ymax
        seqmag(imag)%zmin=zmin
        seqmag(imag)%zmax=zmax

        write(lungfo,*)
        write(lungfo,*)' seqmag(imag)%nx, seqmag(imag)%ny, seqmag(imag)%nz, and number of data lines:',
     &    seqmag(imag)%nx,seqmag(imag)%ny,seqmag(imag)%nz,seqmag(imag)%ntot
        write(lungfo,*)' xmin, xmax                   :',sngl(xmin),sngl(xmax)
        write(lungfo,*)' ymin, ymax , step size of map:',sngl(ymin),sngl(ymax),sngl(seqmag(imag)%bmapdy)
        write(lungfo,*)' zmin, zmax , step size of map:',sngl(zmin),sngl(zmax),sngl(seqmag(imag)%bmapdz)
        write(lungfo,*)' Bxmin, Bxmax of map:',sngl(seqmag(imag)%bmbxmin),sngl(seqmag(imag)%bmbxmax)
        write(lungfo,*)' Bymin, Bymax of map:',sngl(seqmag(imag)%bmbymin),sngl(seqmag(imag)%bmbymax)
        write(lungfo,*)' Bzmin, Bzmax of map:',sngl(seqmag(imag)%bmbzmin),sngl(seqmag(imag)%bmbzmax)
        write(lungfo,*)

        seqmag(imag)%ical=1
      endif

      if (bxout.eq.-9999.0d0) then
        xin=seqmag(imag)%bmxmin+seqmag(imag)%xcen
        return
      else if (bxout.eq.9999.0d0) then
        xin=seqmag(imag)%bmxmax-seqmag(imag)%xcen
        return
      endif

      x=xin-seqmag(imag)%xcen
      y=yin-seqmag(imag)%ycen
      z=zin-seqmag(imag)%zcen

      IF (seqmag(imag)%iwarNX.EQ.0.AND.(X.LT.seqmag(imag)%bmXMIN-STEP.OR.X.GT.seqmag(imag)%bmXMAX+STEP)) THEN
c        WRITE(6,*)'*** WARNING: IN BMAPSEQ: X OUT OF RANGE'
c        WRITE(6,*)'X:',X
c        WRITE(6,*)'Y:',Y
c        WRITE(6,*)'Z:',Z
c        WRITE(LUNGFO,*)'*** WARNING IN BMAPSEQ: X OUT OF RANGE'
c        WRITE(LUNGFO,*)'X:',X
c        WRITE(LUNGFO,*)'Y:',Y
c        WRITE(LUNGFO,*)'Z:',Z
        seqmag(imag)%iwarNX=1
      ENDIF

      IF (seqmag(imag)%iwarNX.NE.0.AND.(X.LT.seqmag(imag)%bmXMIN-STEP.OR.X.GT.seqmag(imag)%bmXMAX+STEP)) THEN
        BXOUT=0.0D0
        BYOUT=0.0D0
        BZOUT=0.0D0
        seqmag(imag)%iwarnbmap=1
        RETURN
      ENDIF

      IF (seqmag(imag)%ny.gt.1.and.(Y.LT.seqmag(imag)%bmYMIN-STEPy.OR.Y.GT.seqmag(imag)%bmYMAX+STEPy)) THEN
        if (seqmag(imag)%iwarny.eq.0) then
          WRITE(6,*)'*** ERROR IN BMAPSEQ: Y OUT OF RANGE'
          WRITE(6,*)'X:',X
          WRITE(6,*)'Y:',Y
          WRITE(6,*)'Z:',Z
          WRITE(LUNGFO,*)'*** ERROR IN BMAPSEQ: Y OUT OF RANGE'
          WRITE(LUNGFO,*)'X:',X
          WRITE(LUNGFO,*)'Y:',Y
          WRITE(LUNGFO,*)'Z:',Z
          seqmag(imag)%iwarny=1
        endif
        BXOUT=0.0D0
        BYOUT=0.0D0
        BZOUT=0.0D0
        seqmag(imag)%iwarnbmap=1
        RETURN
c        STOP
      ENDIF

      IF (seqmag(imag)%nz.gt.1.and.(Z.LT.seqmag(imag)%bmZMIN-STEP.OR.Z.GT.seqmag(imag)%bmZMAX+STEP)) THEN
        if (seqmag(imag)%iwarnz.eq.0) then
          WRITE(6,*)'*** ERROR IN BMAPSEQ: Z OUT OF RANGE'
          WRITE(6,*)'X:',X
          WRITE(6,*)'Y:',Y
          WRITE(6,*)'Z:',Z
          WRITE(LUNGFO,*)'*** ERROR IN BMAPSEQ: Z OUT OF RANGE'
          WRITE(LUNGFO,*)'X:',X
          WRITE(LUNGFO,*)'Y:',Y
          WRITE(LUNGFO,*)'Z:',Z
          seqmag(imag)%iwarnz=1
        endif
        BXOUT=0.0D0
        BYOUT=0.0D0
        BZOUT=0.0D0
        seqmag(imag)%iwarnbmap=1
        RETURN
c        STOP
      ENDIF

      if (x.lt.seqmag(imag)%bmxmin) then
        ix1=1
        ix2=2
      else if (x.gt.seqmag(imag)%bmxmax) then
        ix1=seqmag(imag)%nx-1
        ix2=seqmag(imag)%nx
      else

        if (x.ge.seqmag(imag)%bmap(1,(ix1-1)*seqmag(imag)%nyz+1)) then
c hunt up
          kd=1
111       ix2=min(ix1+kd,seqmag(imag)%nx)
          if (x.gt.seqmag(imag)%bmap(1,seqmag(imag)%nyz*(ix2-1)+1)) then
            kd=2*kd
            ix1=ix2
            goto 111
          endif
        else    !(x.ge.seqmag(imag)%bmap(1,ix1))
c hunt down
          kd=1
          ix2=ix1
22        ix1=max(ix2-kd,1)
          if (x.lt.seqmag(imag)%bmap(1,(ix1-1)*seqmag(imag)%nyz+1)) then
            kd=2*kd
            ix2=ix1
            goto 22
          endif
        endif

1111    if (ix2-ix1.gt.1) then
          k=(ix2+ix1)/2
          if(seqmag(imag)%bmap(1,(k-1)*seqmag(imag)%nyz+1).gt.x)then
            ix2=k
          else
            ix1=k
          endif
          goto 1111
        endif
      endif

      x1=seqmag(imag)%bmap(1,(ix1-1)*seqmag(imag)%nyz+1)
      x2=seqmag(imag)%bmap(1,(ix2-1)*seqmag(imag)%nyz+1)
      dxx=(x-x1)/(x2-x1)

      if (seqmag(imag)%ny.gt.1) then
        if (y.lt.seqmag(imag)%bmymin) then
          iy1=1
        else
          iy1=int((y-seqmag(imag)%bmymin)/seqmag(imag)%bmapdy)+1
          iy1=max(1,iy1)
        endif
        iy2=iy1+1
        if (iy2.gt.seqmag(imag)%ny) then
          iy2=seqmag(imag)%ny
          iy1=iy2-1
        endif
      else
        iy1=1
        iy2=1
      endif

      if (seqmag(imag)%nz.gt.1) then
        if (z.lt.seqmag(imag)%bmzmin) then
          iz1=1
        else
          iz1=int((z-seqmag(imag)%bmzmin)/seqmag(imag)%bmapdz)+1
          iz1=max(1,iz1)
        endif
        iz2=iz1+1
        if (iz2.gt.seqmag(imag)%nz) then
          iz2=seqmag(imag)%nz
          iz1=iz2-1
        endif
      else
        iz1=1
        iz2=1
      endif

      if (irfilb0.eq.6) then

        dyy=(y-(seqmag(imag)%bmymin+(iy1-1)*seqmag(imag)%bmapdy))/seqmag(imag)%bmapdy
        dzz=(z-(seqmag(imag)%bmzmin+(iz1-1)*seqmag(imag)%bmapdz))/seqmag(imag)%bmapdz

        do i=1,3

          b111(i)=seqmag(imag)%bmap(3+i,iz1+(iy1-1)*seqmag(imag)%nz+(ix1-1)*seqmag(imag)%nyz)
          b211(i)=seqmag(imag)%bmap(3+i,iz2+(iy1-1)*seqmag(imag)%nz+(ix1-1)*seqmag(imag)%nyz)
          b121(i)=seqmag(imag)%bmap(3+i,iz1+(iy2-1)*seqmag(imag)%nz+(ix1-1)*seqmag(imag)%nyz)
          b221(i)=seqmag(imag)%bmap(3+i,iz2+(iy2-1)*seqmag(imag)%nz+(ix1-1)*seqmag(imag)%nyz)
          b112(i)=seqmag(imag)%bmap(3+i,iz1+(iy1-1)*seqmag(imag)%nz+(ix2-1)*seqmag(imag)%nyz)
          b212(i)=seqmag(imag)%bmap(3+i,iz2+(iy1-1)*seqmag(imag)%nz+(ix2-1)*seqmag(imag)%nyz)
          b122(i)=seqmag(imag)%bmap(3+i,iz1+(iy2-1)*seqmag(imag)%nz+(ix2-1)*seqmag(imag)%nyz)
          b222(i)=seqmag(imag)%bmap(3+i,iz2+(iy2-1)*seqmag(imag)%nz+(ix2-1)*seqmag(imag)%nyz)

          b112111(i)=b111(i)+(b112(i)-b111(i))*dxx
          b122121(i)=b121(i)+(b122(i)-b121(i))*dxx
          b212211(i)=b211(i)+(b212(i)-b211(i))*dxx
          b222221(i)=b221(i)+(b222(i)-b221(i))*dxx

          blow(i)=b112111(i)+(b212211(i)-b112111(i))*dzz
          bhig(i)=b122121(i)+(b222221(i)-b122121(i))*dzz

          b(i)=blow(i)+(bhig(i)-blow(i))*dyy

        enddo

      else !ifilb0

        if (ix1.gt.1) then
          kx1=ix1-1
          kx2=ix1
          kx3=ix1+1
        else
          kx1=ix1
          kx2=ix1+1
          kx3=ix1+2
        endif

        if (seqmag(imag)%ny.gt.1) then
          if (iy1.gt.1) then
            ky1=iy1-1
            ky2=iy1
            ky3=iy1+1
          else
            ky1=iy1
            ky2=iy1+1
            ky3=iy1+2
          endif
        else
          ky1=1
          ky2=1
          ky3=1
        endif

        if (seqmag(imag)%nz.gt.1) then
          if (iz1.gt.1) then
            kz1=iz1-1
            kz2=iz1
            kz3=iz1+1
          else
            kz1=iz1
            kz2=iz1+1
            kz3=iz1+2
          endif
        else
          kz1=1
          kz2=1
          kz3=1
        endif

        x3(1)=0.0d0

        x1=seqmag(imag)%bmap(1,kz1+(ky1-1)*seqmag(imag)%nz+(kx1-1)*seqmag(imag)%nyz)
        xx=x-x1

        y1=seqmag(imag)%bmap(2,kz1+(ky1-1)*seqmag(imag)%nz+(kx1-1)*seqmag(imag)%nyz)
        yy=y-y1
        dy=seqmag(imag)%bmap(2,kz1+(ky2-1)*seqmag(imag)%nz+(kx1-1)*seqmag(imag)%nyz)-y1

        z1=seqmag(imag)%bmap(3,kz1+(ky1-1)*seqmag(imag)%nz+(kx1-1)*seqmag(imag)%nyz)
        zz=z-z1
        dz=seqmag(imag)%bmap(3,kz2+(ky1-1)*seqmag(imag)%nz+(kx1-1)*seqmag(imag)%nyz)-z1

        do i=1,3

          x3(2)=seqmag(imag)%bmap(1,kz1+(ky1-1)*seqmag(imag)%nz+(kx2-1)*seqmag(imag)%nyz)-x1
          x3(3)=seqmag(imag)%bmap(1,kz1+(ky1-1)*seqmag(imag)%nz+(kx3-1)*seqmag(imag)%nyz)-x1

          b3(1)=seqmag(imag)%bmap(3+i,kz1+(ky1-1)*seqmag(imag)%nz+(kx1-1)*seqmag(imag)%nyz)
          b3(2)=seqmag(imag)%bmap(3+i,kz1+(ky1-1)*seqmag(imag)%nz+(kx2-1)*seqmag(imag)%nyz)
          b3(3)=seqmag(imag)%bmap(3+i,kz1+(ky1-1)*seqmag(imag)%nz+(kx3-1)*seqmag(imag)%nyz)
          call parabel_short(x3,b3,a3)
          b11=a3(1)+(a3(2)+a3(3)*xx)*xx

          b3(1)=seqmag(imag)%bmap(3+i,kz1+(ky2-1)*seqmag(imag)%nz+(kx1-1)*seqmag(imag)%nyz)
          b3(2)=seqmag(imag)%bmap(3+i,kz1+(ky2-1)*seqmag(imag)%nz+(kx2-1)*seqmag(imag)%nyz)
          b3(3)=seqmag(imag)%bmap(3+i,kz1+(ky2-1)*seqmag(imag)%nz+(kx3-1)*seqmag(imag)%nyz)
          call parabel_short(x3,b3,a3)
          b21=a3(1)+(a3(2)+a3(3)*xx)*xx

          b3(1)=seqmag(imag)%bmap(3+i,kz1+(ky3-1)*seqmag(imag)%nz+(kx1-1)*seqmag(imag)%nyz)
          b3(2)=seqmag(imag)%bmap(3+i,kz1+(ky3-1)*seqmag(imag)%nz+(kx2-1)*seqmag(imag)%nyz)
          b3(3)=seqmag(imag)%bmap(3+i,kz1+(ky3-1)*seqmag(imag)%nz+(kx3-1)*seqmag(imag)%nyz)
          call parabel_short(x3,b3,a3)
          b31=a3(1)+(a3(2)+a3(3)*xx)*xx

          b3(1)=seqmag(imag)%bmap(3+i,kz2+(ky1-1)*seqmag(imag)%nz+(kx1-1)*seqmag(imag)%nyz)
          b3(2)=seqmag(imag)%bmap(3+i,kz2+(ky1-1)*seqmag(imag)%nz+(kx2-1)*seqmag(imag)%nyz)
          b3(3)=seqmag(imag)%bmap(3+i,kz2+(ky1-1)*seqmag(imag)%nz+(kx3-1)*seqmag(imag)%nyz)
          call parabel_short(x3,b3,a3)
          b12=a3(1)+(a3(2)+a3(3)*xx)*xx

          b3(1)=seqmag(imag)%bmap(3+i,kz2+(ky2-1)*seqmag(imag)%nz+(kx1-1)*seqmag(imag)%nyz)
          b3(2)=seqmag(imag)%bmap(3+i,kz2+(ky2-1)*seqmag(imag)%nz+(kx2-1)*seqmag(imag)%nyz)
          b3(3)=seqmag(imag)%bmap(3+i,kz2+(ky2-1)*seqmag(imag)%nz+(kx3-1)*seqmag(imag)%nyz)
          call parabel_short(x3,b3,a3)
          b22=a3(1)+(a3(2)+a3(3)*xx)*xx

          b3(1)=seqmag(imag)%bmap(3+i,kz2+(ky3-1)*seqmag(imag)%nz+(kx1-1)*seqmag(imag)%nyz)
          b3(2)=seqmag(imag)%bmap(3+i,kz2+(ky3-1)*seqmag(imag)%nz+(kx2-1)*seqmag(imag)%nyz)
          b3(3)=seqmag(imag)%bmap(3+i,kz2+(ky3-1)*seqmag(imag)%nz+(kx3-1)*seqmag(imag)%nyz)
          call parabel_short(x3,b3,a3)
          b32=a3(1)+(a3(2)+a3(3)*xx)*xx

          b3(1)=seqmag(imag)%bmap(3+i,kz3+(ky1-1)*seqmag(imag)%nz+(kx1-1)*seqmag(imag)%nyz)
          b3(2)=seqmag(imag)%bmap(3+i,kz3+(ky1-1)*seqmag(imag)%nz+(kx2-1)*seqmag(imag)%nyz)
          b3(3)=seqmag(imag)%bmap(3+i,kz3+(ky1-1)*seqmag(imag)%nz+(kx3-1)*seqmag(imag)%nyz)
          call parabel_short(x3,b3,a3)
          b13=a3(1)+(a3(2)+a3(3)*xx)*xx

          b3(1)=seqmag(imag)%bmap(3+i,kz3+(ky2-1)*seqmag(imag)%nz+(kx1-1)*seqmag(imag)%nyz)
          b3(2)=seqmag(imag)%bmap(3+i,kz3+(ky2-1)*seqmag(imag)%nz+(kx2-1)*seqmag(imag)%nyz)
          b3(3)=seqmag(imag)%bmap(3+i,kz3+(ky2-1)*seqmag(imag)%nz+(kx3-1)*seqmag(imag)%nyz)
          call parabel_short(x3,b3,a3)
          b23=a3(1)+(a3(2)+a3(3)*xx)*xx

          b3(1)=seqmag(imag)%bmap(3+i,kz3+(ky3-1)*seqmag(imag)%nz+(kx1-1)*seqmag(imag)%nyz)
          b3(2)=seqmag(imag)%bmap(3+i,kz3+(ky3-1)*seqmag(imag)%nz+(kx2-1)*seqmag(imag)%nyz)
          b3(3)=seqmag(imag)%bmap(3+i,kz3+(ky3-1)*seqmag(imag)%nz+(kx3-1)*seqmag(imag)%nyz)
          call parabel_short(x3,b3,a3)
          b33=a3(1)+(a3(2)+a3(3)*xx)*xx

          x3(2)=dy
          x3(3)=x3(2)+dy

          if (seqmag(imag)%ny.gt.1) then

            b3(1)=b11
            b3(2)=b21
            b3(3)=b31
            call parabel_short(x3,b3,a3)
            bb1=a3(1)+(a3(2)+a3(3)*yy)*yy

            b3(1)=b12
            b3(2)=b22
            b3(3)=b32
            call parabel_short(x3,b3,a3)
            bb2=a3(1)+(a3(2)+a3(3)*yy)*yy

            b3(1)=b13
            b3(2)=b23
            b3(3)=b33
            call parabel_short(x3,b3,a3)
            bb3=a3(1)+(a3(2)+a3(3)*yy)*yy

          else
            bb3=b11
          endif

          if (seqmag(imag)%nz.gt.1) then
            x3(2)=dz
            x3(3)=x3(2)+dz
            b3(1)=bb1
            b3(2)=bb2
            b3(3)=bb3
            call parabel_short(x3,b3,a3)
            b(i)=a3(1)+(a3(2)+a3(3)*zz)*zz
          else
            b(i)=bb1
          endif

        enddo !i=1,3

      endif !irfilb0

      bxout=b(1)
      byout=b(2)
      bzout=b(3)

      return
      end
