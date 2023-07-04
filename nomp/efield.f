*CMZ :  2.66/20 21/11/2011  17.24.41  by  Michael Scheer
*-- Author :    Michael Scheer   18/11/2011
      subroutine efield(x,y,z,ex,ey,ez)
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

c calculates electrical field (ex,ey,ez) in V/m and corresponding
c dgamma for (x,y,z,dt)
c

c data format:
c x1 y1 z1 ex1 ey1 ez1
c       .
c x1 y1 znz ex111 ey111 ez11nx
c       .
c       .
c x1 y2 z1 ex121 ey121 ez121
c       .
c x1 yny znz ex1nynz ey1nynz ez1nynz
c       .
c       .
c xnx yny znz exnxnynz eynxnynz eznxnynz


*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,phycon.
      include 'phycon.cmn'
*KEND.

      integer nedimp
      parameter (nedimp=100000)

      double precision x,y,z,xx,yy,zz,
     &  ex,ey,ez,f(6,nedimp),
     &  xmax,xmin,ymax,ymin,zmax,zmin,
     &  exmax,exmin,eymax,eymin,ezmax,ezmin,
     &  x1,x2,y1,y2,z1,z2

      integer ical,lunein,nx,ny,nz,ix,iy,iz,ieof,nread,ixold,iyold,izold,k,
     &  klo,khi,kd,nzy

      data ical/0/
      data lunein/0/

      if (ical.eq.0) then
        nread=0
        nx=0
        x1=1.0d30
        y1=1.0d30
        z1=1.0d30
        x2=-1.0d30
        y2=-1.0d30
        z2=-1.0d30
        open(unit=lunein,file='efield.dat')
1       call util_skip_comment_end(lunein,ieof)
        if (ieof.ne.0) goto 9
        read(lunein,*)xx,yy,zz,ex,ey,ez
        nread=nread+1
        if (nread.gt.nedimp) stop '*** Error in efield: Dimension exceeded ***'
        f(1,nread)=xx
        f(2,nread)=yy
        f(3,nread)=zz
        f(4,nread)=ex
        f(5,nread)=ey
        f(6,nread)=ez
        if (xx.gt.xmax) xmax=xx
        if (yy.gt.ymax) ymax=yy
        if (zz.gt.zmax) zmax=zz
        if (xx.lt.xmin) xmin=xx
        if (yy.lt.ymin) ymin=yy
        if (zz.lt.zmin) zmin=zz
        if (ex.gt.exmax) exmax=ex
        if (ey.gt.eymax) eymax=ey
        if (ez.gt.ezmax) ezmax=ez
        if (ex.lt.exmin) exmin=ex
        if (ey.lt.eymin) eymin=ey
        if (ez.lt.ezmin) ezmin=ez
        if (f(1,nread).ne.x1) then
          nx=nx+1
          x1=xx
        endif
        goto 1
9       close(lunein)

        z1=f(3,1)
        do k=2,nread
          if (f(3,k).ne.z1) then
            nz=k
            goto 99
          endif
        enddo
99      ny=nread/nx/nz

        write(lungfo,*)' '
        write(lungfo,*)'      E-Field routine efield called first time:'
        write(lungfo,*)' '
        if (nx.lt.2) stop '*** Error in efield: Less then two x-values found ***'
        if (ny.lt.2) stop '*** Error in efield: Less then two y-values found ***'
        if (nz.lt.2) stop '*** Error in efield: Less then two z-values found ***'
        write(lungfo,*)'         Number of data points:',nread
        write(lungfo,*)' '
        write(lungfo,*)'         xmin, xmax:',xmin,xmax
        write(lungfo,*)'         ymin, ymax:',ymin,ymax
        write(lungfo,*)'         zmin, zmax:',zmin,zmax
        write(lungfo,*)' '
        write(lungfo,*)'         Exmin, Exmax:',Exmin,Exmax
        write(lungfo,*)'         Eymin, Eymax:',Eymin,Eymax
        write(lungfo,*)'         Ezmin, Ezmax:',Ezmin,Ezmax
        write(lungfo,*)' '
        if (nx*ny*nz.ne.nread) stop '*** Error in efield: Invalid file efield.dat ***'
      endif !ical

      klo=1
      khi=nz
      kd=nz
      if (ical.gt.0) then
        z1=f(3,izold)
        z2=f(3,izold+1)
      endif

      if (z.ge.z1.and.z.lt.z2) then
        iz=izold
      else !
        if (ical.ne.0) then
          if (z.ge.z1) then
            KD=1
121         KHI=MIN(KLO+KD,nz)
            IF (z.GT.f(3,KHI)) THEN
              KD=2*KD
              KLO=KHI
              GOTO 121
            ENDIF
          ELSE
            KD=1
            KHI=KLO
122         KLO=MAX(KHI-KD,1)
            IF (z.LT.f(3,KLO)) THEN
              KD=2*KD
              KHI=KLO
              GOTO 122
            ENDIF
          endif
        endif !ical
131     IF (KHI-KLO.GT.1) THEN
          K=(KHI+KLO)/2
          IF(f(3,K).GT.z)THEN
            KHI=K
          ELSE
            KLO=K
          ENDIF
          GOTO 131
        ENDIF
        iz=klo
      endif !z

      klo=1
      khi=ny
      kd=ny
      if (ical.ne.0) then
        y1=f(2,iyold)
        y2=f(2,iyold+nz)
      endif

      if (y.ge.y1.and.y.lt.y2) then
        iy=iyold
      else !
        if (ical.ne.0) then
          if (y.ge.y1) then
            KD=1
221         KHI=MIN(KLO+KD,ny)
            IF (y.GT.f(2,khi+nz)) THEN
              KD=2*KD
              KLO=KHI
              GOTO 221
            ENDIF
          ELSE
            KD=1
            KHI=KLO
222         KLO=MAX(KHI-KD,1)
            IF (y.LT.f(2,klo+nz)) THEN
              KD=2*KD
              KHI=KLO
              GOTO 222
            ENDIF
          endif
        endif !ical
231     IF (KHI-KLO.GT.1) THEN
          K=(KHI+KLO)/2
          IF(f(2,k+nz).GT.y)THEN
            KHI=K
          ELSE
            KLO=K
          ENDIF
          GOTO 231
        ENDIF
        iy=klo
      endif !y

      klo=1
      khi=nx
      kd=nx
      nzy=nz+ny
      if (ical.ne.0) then
        x1=f(1,ixold)
        x2=f(1,ixold+nzy)
      endif

      if (x.ge.x1.and.x.lt.x2) then
        ix=ixold
      else !
        if (ical.ne.0) then
          if (x.ge.x1) then
            KD=1
21          KHI=MIN(KLO+KD,nx)
            IF (X.GT.f(1,KHI+nzy)) THEN
              KD=2*KD
              KLO=KHI
              GOTO 21
            ENDIF
          ELSE
            KD=1
            KHI=KLO
22          KLO=MAX(KHI-KD,1)
            IF (X.LT.f(1,KLO+nzy)) THEN
              KD=2*KD
              KHI=KLO
              GOTO 22
            ENDIF
          endif
        endif !ical
31      IF (KHI-KLO.GT.1) THEN
          K=(KHI+KLO)/2
          IF(f(1,K+nzy).GT.X)THEN
            KHI=K
          ELSE
            KLO=K
          ENDIF
          GOTO 31
        ENDIF
        ix=klo
      endif !x

      ixold=ix
      iyold=iy
      izold=iz

      ical=1

      return
      end
