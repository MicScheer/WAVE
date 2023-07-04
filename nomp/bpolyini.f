*CMZ :  3.05/00 25/04/2018  13.09.51  by  Michael Scheer
*CMZ :  3.01/02 27/08/2013  09.19.16  by  Michael Scheer
*CMZ :  2.62/01 24/04/2007  12.00.02  by  Michael Scheer
*CMZ :  2.62/00 16/04/2007  16.00.14  by  Michael Scheer
*CMZ :  2.57/05 01/08/2006  19.29.50  by  Michael Scheer
*CMZ :  2.54/02 13/04/2005  12.35.45  by  Michael Scheer
*CMZ :  2.53/05 23/02/2005  18.54.13  by  Michael Scheer
*CMZ :  2.52/05 16/08/2004  14.49.57  by  Michael Scheer
*CMZ :  2.48/00 01/03/2004  17.55.46  by  Michael Scheer
*CMZ :  2.47/23 17/02/2004  13.47.20  by  Michael Scheer
*CMZ :  0.99/07 16/02/2004  17.22.16  by  Michael Scheer
*CMZ :  0.99/03 12/02/2004  13.55.05  by  Michael Scheer
*CMZ :  0.99/00 26/01/2004  17.03.49  by  Michael Scheer
*CMZ :  0.00/08 23/01/2004  12.52.23  by  Michael Scheer
*CMZ :  0.00/07 16/01/2004  11.05.44  by  Michael Scheer
*CMZ :  0.00/06 09/01/2004  15.55.17  by  Michael Scheer
*CMZ :  0.00/05 23/12/2003  14.52.54  by  Michael Scheer
*CMZ :  0.00/04 23/12/2003  10.15.07  by  Michael Scheer
*CMZ :  0.00/02 15/12/2003  12.43.34  by  Michael Scheer
*CMZ :  0.00/01 10/12/2003  17.56.52  by  Michael Scheer
*-- Author :    Michael Scheer   02/12/2003
      subroutine bpolyini(xstart,xstop,lungfo)
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

*KEEP,bwpolyederf90u.
      include 'bwpolyederf90u.cmn'
*KEND.

      implicit none

*KEEP,bforce.
      include 'bforce.cmn'
*KEND.

      double precision, dimension (:), allocatable :: xsort,ysort,bpexpos0

      integer, dimension (:), allocatable :: ibpeplan0
      integer, dimension (:,:), allocatable :: ibpecorn0

      double precision, dimension (:,:), allocatable :: bpebc0
      double precision, dimension (:,:,:), allocatable :: bpemat0
      double precision, dimension (:,:,:,:), allocatable :: bpemag0,shuffle,
     &  bperot0,bpetm0

      integer, dimension (:), allocatable :: ibpecol

      double precision x0,y0,z0,bc,xm,ym,zm,rmag(3),vmaglab(3),
     &  vnormlab(3),r1(3),bdum,dum,r1lab(3),
     &  xstart,xstop,xmn,xmx,
     &  rausch,offrausch

      double precision ts(3,3),tsinv(3,3),p1(3),p2(3),p3(3)
      double precision v1x,v1y,v1z,v2x,v2y,v2z,vsx,vsy,vsz
      double precision xmin,xmax,ymin,ymax,zmin,zmax
      double precision space,vspace(3),vbsym(3),rotmod(3,3),
     &  xmod,ymod,zmod,q,qsign,a,b,r2(3),det,
     &  xlen,ylen,zlen,vx,vy,vn,ca,sa,tz(3,3),ws(3,3),
     &  ts1(3,3),ts1inv(3,3),
     &  bxi,byi,bzi,xint,yint,zint,
     &  vxint,vyint,vzint,
     &  bxint,byint,bzint,ddum,rmin,rmax,y0max,y0min,ylenmax,ylenmin,ycen,
     &  shll,shlr,shul,shur,hgap

      real xplmin,xplmax,yplmin,yplmax,zplmin,zplmax,theta,phi

      integer iwarn,ifilemodus,ifail,lungfo,iplot,idum,
     &  iallowin,ishfield,ibextern,iseed

      integer ntupp,nx,ny,nz,nmagdiv,ndiv,idiv,n2div
      parameter(ntupp=10)

      real htup(ntupp)

      integer ical
      integer nplan,ncol,ncorn,nmodule,ncopy
      integer nmagmax,nplanmax,ncornmax,kmag1,kmag2
      integer imag,jmag,kmag,iplan,icorn,imodul,icopy,i,
     &      j,ip1,ip2,iflange,nflange,nmago

      character(256) usercom,cfile

      data ical/0/

      if (ical.eq.0) then

        tiny=1.0d-6

        lunbpe=lunpm
        filebpe=filepm

        open(unit=lunbpe,file=filebpe,form='formatted',status='old')

        call util_skip_comment(lunbpe)
        read(lunbpe,*)ifilemodus,iplot,
     &    xplmin,xplmax,yplmin,yplmax,zplmin,zplmax,theta,phi

        call util_skip_comment(lunbpe)
        read(lunbpe,'(a)')usercom

        call util_skip_comment(lunbpe)
        read(lunbpe,*)iaxint,x0int,y0int,z0int,nstepint,ranginti,ranginte

        call util_skip_comment(lunbpe)
        read(lunbpe,*)idum,idum,idum,iforcol

        call util_skip_comment(lunbpe)
        read(lunbpe,*)ddum,ddum,ddum

        call util_skip_comment(lunbpe)
        read(lunbpe,*)ddum,ddum,ddum

        call util_skip_comment(lunbpe)
        read(lunbpe,*)ddum,ddum,ddum

        call util_skip_comment(lunbpe)
        read(lunbpe,*)idum,idum,ddum

        call util_skip_comment(lunbpe)
        read(lunbpe,*)ddum,ddum,ddum,ddum

        call util_skip_comment(lunbpe)
        read(lunbpe,*)ddum,ddum,ddum

        call util_skip_comment(lunbpe)
        read(lunbpe,*)nx,xmin,xmax

        call util_skip_comment(lunbpe)
        read(lunbpe,*)ny,ymin,ymax

        call util_skip_comment(lunbpe)
        read(lunbpe,*)nz,zmin,zmax

        call util_skip_comment(lunbpe)
        read(lunbpe,*)ddum

        nmagmax=0
        nplanmax=0
        ncornmax=0

11      call util_skip_comment(lunbpe)
        read(lunbpe,*)nmag

        nmago=nmag

        if (nmag.eq.0) then
          call util_skip_comment(lunbpe)
          read(lunbpe,'(a)')cfile
          lunbpe=197
          open(unit=lunbpe,file=cfile,form='formatted',status='old')
          call util_skip_comment(lunbpe)
          read(lunbpe,*)nmag
        endif !nmago

        if (nmag.eq.9999) goto 99

c get dimensions

        do imag=1,nmag

          call util_skip_comment(lunbpe)
          read(lunbpe,*)x0,y0,z0

          call util_skip_comment(lunbpe)
          read(lunbpe,*)bc,x0,y0,z0

          call util_skip_comment(lunbpe)
          read(lunbpe,*)nplan,ncol

          if ((nplan.eq.-1.or.nplan.eq.-6).and.nplanmax.lt.6) then
            nplanmax=6
          else if (abs(nplan).gt.nplanmax) then
            nplanmax=abs(nplan)
          endif

          if (nplan.gt.0) then

            do iplan=1,nplan

              call util_skip_comment(lunbpe)
              read(lunbpe,*)ncorn

              if (ncorn.lt.3) then
                print *
                print *,'*** ERROR IN BPOLYINI: Too few points'
                print *,'Magnet, plane: ',imag,iplan
                print *
                stop
              endif

              if (ncorn+1.gt.ncornmax) ncornmax=ncorn+1

              do icorn=1,ncorn
                call util_skip_comment(lunbpe)
                read(lunbpe,*)x0,y0,z0
              enddo

              ncorn=ncorn+1

            enddo !nplan

          else if (nplan.eq.-6) then

              ncorn=4
              if (ncorn+1.gt.ncornmax) ncornmax=ncorn+1
              ncorn=ncorn+1

              call util_skip_comment(lunbpe)
              read(lunbpe,*)xlen,ylen,zlen

          else if (nplan.eq.-1) then

            ncorn=4
            if (ncorn+1.gt.ncornmax) ncornmax=ncorn+1
            ncorn=ncorn+1

            call util_skip_comment(lunbpe)
            read(lunbpe,*)xlen,ylen,zlen,ndiv

            if (xlen.gt.0.0d0) then
              nmagdiv=nmagdiv+2*ndiv-1
              n2div=2*ndiv
            else
              nmagdiv=nmagdiv+ndiv-1
              n2div=ndiv
            endif

          else !(nplan.eq.-6)
            print*, '*** Error in BPOLYINI: Bad value of number of planes'
            print*, '*** Magnet:', imag
            stop '*** Program aborted ***'
          endif !nplan.gt.0

        enddo !nmag

        nmag=nmag+nmagdiv

        lunbpe=lunpm

        call util_skip_comment(lunbpe)
        read(lunbpe,*)nmodule

        do imodul=1,nmodule

          call util_skip_comment(lunbpe)
          read(lunbpe,*)x0,y0,z0

          call util_skip_comment(lunbpe)
          read(lunbpe,*)x0,y0,z0

          call util_skip_comment(lunbpe)
          read(lunbpe,*)x0,y0,z0

          call util_skip_comment(lunbpe)
          read(lunbpe,*)x0,y0,z0

          call util_skip_comment(lunbpe)
          read(lunbpe,*)ncopy

          call util_skip_comment(lunbpe)
          read(lunbpe,*)space,vspace(1),vspace(2),vspace(3)

          call util_skip_comment(lunbpe)
          read(lunbpe,*)vbsym(1),vbsym(2),vbsym(3)

          nmagmax=nmagmax+ncopy*nmag

        enddo !imodul=1,nomodul

        goto 11

99      continue

        call util_skip_comment(lunbpe)
        rausch=0.0
        offrausch=0.0

        read(lunbpe,*,end=991)rausch,offrausch
        offrausch=offrausch/1000.

        read(lunbpe,*)iseed

        call util_skip_comment(lunbpe)
        read(lunbpe,*,end=991)iallowin

        call util_skip_comment(lunbpe)
        read(lunbpe,*,end=991)ishfield

        call util_skip_comment(lunbpe)
        read(lunbpe,*,end=991)ibextern

        call util_skip_comment(lunbpe)
        read(lunbpe,*,end=991)shll
        call util_skip_comment(lunbpe)
        read(lunbpe,*,end=991)shlr
        call util_skip_comment(lunbpe)
        read(lunbpe,*,end=991)shul
        call util_skip_comment(lunbpe)
        read(lunbpe,*,end=991)shur
        call util_skip_comment(lunbpe)
c        read(lunbpe,*,end=991)gappm
        read(lunbpe,*,end=991)dum

c        shiftll=shiftll+shll
c        shiftlr=shiftlr+shlr
c        shiftul=shiftul+shul
c        shiftur=shiftur+shur

991     continue

        allocate(bpexpos(nmagmax))
        allocate(bpebc(8,nmagmax))
        allocate(bpemat(3,3,nmagmax))
        allocate(bpemag(6,ncornmax,nplanmax,nmagmax))
        allocate(shuffle(7,ncornmax,nplanmax,nmagmax))
        allocate(bperot(3,ncornmax,nplanmax,nmagmax))
        allocate(bpetm(3,8,nplanmax,nmagmax))
        allocate(ibpeplan(nmagmax))
        allocate(ibpecol(nmagmax))
        allocate(ibpecorn(nplanmax,nmagmax))
        allocate(bflange(7,ncornmax*nplanmax))

        allocate(bpexpos0(nmagmax))
        allocate(bpebc0(8,nmagmax))
        allocate(bpemat0(3,3,nmagmax))
        allocate(bpemag0(6,ncornmax,nplanmax,nmagmax))
        allocate(bperot0(3,ncornmax,nplanmax,nmagmax))
        allocate(bpetm0(3,8,nplanmax,nmagmax))
        allocate(ibpeplan0(nmagmax))
        allocate(ibpecorn0(nplanmax,nmagmax))

        rewind(lunbpe)

c read and store data

        nmagmax=0
        nplanmax=0
        ncornmax=0

        call util_skip_comment(lunbpe)
        read(lunbpe,*)ifilemodus,iplot,
     &    xplmin,xplmax,yplmin,yplmax,zplmin,zplmax,theta,phi

        call util_skip_comment(lunbpe)
        read(lunbpe,'(a)')usercom

        call util_skip_comment(lunbpe)
        read(lunbpe,*)iaxint,x0int,y0int,z0int,nstepint,ranginti,ranginte

        call util_skip_comment(lunbpe)
        read(lunbpe,*)idum,idum,idum,iforcol

        call util_skip_comment(lunbpe)
        read(lunbpe,*)ddum,ddum,ddum

        call util_skip_comment(lunbpe)
        read(lunbpe,*)ddum,ddum,ddum

        call util_skip_comment(lunbpe)
        read(lunbpe,*)ddum,ddum,ddum

        forcol=iforcol

        bflenxmm=bflenx*1000.0d0
        bflenymm=bfleny*1000.0d0
        bflenzmm=bflenz*1000.0d0

        bfcenxmm=bfcenx*1000.0d0
        bfcenymm=bfceny*1000.0d0
        bfcenzmm=bfcenz*1000.0d0

        torqcenxmm=torqcenx*1000.0d0
        torqcenymm=torqceny*1000.0d0
        torqcenzmm=torqcenz*1000.0d0

        forxpl(1)=bfcenxmm-bflenxmm/2.
        forxpl(2)=bfcenxmm+bflenxmm/2.
        forypl(1)=bfcenymm-bflenymm/2.
        forypl(2)=bfcenymm+bflenymm/2.
        forzpl(1)=bfcenzmm-bflenzmm/2.
        forzpl(2)=bfcenzmm+bflenzmm/2.

        call util_skip_comment(lunbpe)
        read(lunbpe,*)idum,idum,ddum

        call util_skip_comment(lunbpe)
        read(lunbpe,*)ddum,ddum,ddum,ddum

        call util_skip_comment(lunbpe)
        read(lunbpe,*)ddum,ddum,ddum

        call util_skip_comment(lunbpe)
        read(lunbpe,*)nx,xmin,xmax

        call util_skip_comment(lunbpe)
        read(lunbpe,*)ny,ymin,ymax

        call util_skip_comment(lunbpe)
        read(lunbpe,*)nz,zmin,zmax

        call util_skip_comment(lunbpe)
        read(lunbpe,*)ddum

1       call util_skip_comment(lunbpe)
        read(lunbpe,*)nmag

        if (nmag.eq.9999) goto 9

        if (nmag.eq.0) then
          call util_skip_comment(lunbpe)
          read(lunbpe,'(a)')cfile
          lunbpe=197
          rewind(lunbpe)
          call util_skip_comment(lunbpe)
          read(lunbpe,*)nmag
        endif !nmago

        kmag1=nmagmax+1
        kmag2=nmagmax+nmag

        do imag=1,nmag

          nmagmax=nmagmax+1

          call util_skip_comment(lunbpe)
          read(lunbpe,*)bpebc0(1,nmagmax),bpebc0(2,nmagmax),bpebc0(3,nmagmax)

          call util_skip_comment(lunbpe)
          read(lunbpe,*)bc,xm,ym,zm !magnetization vector M

          bdum=sqrt(xm*xm+ym*ym+zm*zm)
          if (bdum.eq.0.d0) then
            print *
            print *,'*** ERROR IN BPOLYEDER: Bad magnetization vector'
            print *,'Magnet ',imag
            print *
            stop
          endif
          bdum=bc/bdum

          bpebc0(4,nmagmax)=xm*bdum
          bpebc0(5,nmagmax)=ym*bdum
          bpebc0(6,nmagmax)=zm*bdum

          call util_skip_comment(lunbpe)
          read(lunbpe,*)nplan,ncol

          ibpecol(nmagmax)=ncol

          if (nplan.gt.0) then

            bpebc0(7,nmagmax)=1

            do iplan=1,nplan

              call util_skip_comment(lunbpe)
              read(lunbpe,*)ncorn

              ncorn=ncorn+1

              ibpecorn(iplan,nmagmax)=ncorn

              do icorn=1,ncorn

                if (icorn.lt.ncorn) then

                  call util_skip_comment(lunbpe)
                  read(lunbpe,*)x0,y0,z0

                  bpemag0(1,icorn,iplan,nmagmax)=x0
                  bpemag0(2,icorn,iplan,nmagmax)=y0
                  bpemag0(3,icorn,iplan,nmagmax)=z0

                else ! icorn.lt.ncorn

                  bpemag0(1,icorn,iplan,nmagmax)=bpemag0(1,1,iplan,nmagmax)
                  bpemag0(2,icorn,iplan,nmagmax)=bpemag0(2,1,iplan,nmagmax)
                  bpemag0(3,icorn,iplan,nmagmax)=bpemag0(3,1,iplan,nmagmax)

                endif ! icorn.lt.ncorn

              enddo !icorn

            enddo !iplan

            ibpeplan(nmagmax)=nplan

          else if (nplan.eq.-6) then

            bpebc0(7,nmagmax)=-6

            nplan=6

            call util_skip_comment(lunbpe)
            read(lunbpe,*)xlen,ylen,zlen

            ibpeplan(nmagmax)=6

            iplan=1
            ibpecorn(iplan,nmagmax)=5
            icorn=1
            bpemag0(1,icorn,iplan,nmagmax)=-xlen/2.0d0
            bpemag0(2,icorn,iplan,nmagmax)=-ylen/2.0d0
            bpemag0(3,icorn,iplan,nmagmax)=-zlen/2.0d0
            icorn=2
            bpemag0(1,icorn,iplan,nmagmax)=-xlen/2.0d0
            bpemag0(2,icorn,iplan,nmagmax)=-ylen/2.0d0
            bpemag0(3,icorn,iplan,nmagmax)=+zlen/2.0d0
            icorn=3
            bpemag0(1,icorn,iplan,nmagmax)=-xlen/2.0d0
            bpemag0(2,icorn,iplan,nmagmax)=+ylen/2.0d0
            bpemag0(3,icorn,iplan,nmagmax)=+zlen/2.0d0
            icorn=4
            bpemag0(1,icorn,iplan,nmagmax)=-xlen/2.0d0
            bpemag0(2,icorn,iplan,nmagmax)=+ylen/2.0d0
            bpemag0(3,icorn,iplan,nmagmax)=-zlen/2.0d0

            iplan=2
            ibpecorn(iplan,nmagmax)=5
            icorn=4
            bpemag0(1,icorn,iplan,nmagmax)=+xlen/2.0d0
            bpemag0(2,icorn,iplan,nmagmax)=-ylen/2.0d0
            bpemag0(3,icorn,iplan,nmagmax)=-zlen/2.0d0
            icorn=3
            bpemag0(1,icorn,iplan,nmagmax)=+xlen/2.0d0
            bpemag0(2,icorn,iplan,nmagmax)=-ylen/2.0d0
            bpemag0(3,icorn,iplan,nmagmax)=+zlen/2.0d0
            icorn=2
            bpemag0(1,icorn,iplan,nmagmax)=+xlen/2.0d0
            bpemag0(2,icorn,iplan,nmagmax)=+ylen/2.0d0
            bpemag0(3,icorn,iplan,nmagmax)=+zlen/2.0d0
            icorn=1
            bpemag0(1,icorn,iplan,nmagmax)=+xlen/2.0d0
            bpemag0(2,icorn,iplan,nmagmax)=+ylen/2.0d0
            bpemag0(3,icorn,iplan,nmagmax)=-zlen/2.0d0

            iplan=3
            ibpecorn(iplan,nmagmax)=5
            icorn=4
            bpemag0(1,icorn,iplan,nmagmax)=-xlen/2.0d0
            bpemag0(2,icorn,iplan,nmagmax)=-ylen/2.0d0
            bpemag0(3,icorn,iplan,nmagmax)=-zlen/2.0d0
            icorn=3
            bpemag0(1,icorn,iplan,nmagmax)=-xlen/2.0d0
            bpemag0(2,icorn,iplan,nmagmax)=-ylen/2.0d0
            bpemag0(3,icorn,iplan,nmagmax)=+zlen/2.0d0
            icorn=2
            bpemag0(1,icorn,iplan,nmagmax)=+xlen/2.0d0
            bpemag0(2,icorn,iplan,nmagmax)=-ylen/2.0d0
            bpemag0(3,icorn,iplan,nmagmax)=+zlen/2.0d0
            icorn=1
            bpemag0(1,icorn,iplan,nmagmax)=+xlen/2.0d0
            bpemag0(2,icorn,iplan,nmagmax)=-ylen/2.0d0
            bpemag0(3,icorn,iplan,nmagmax)=-zlen/2.0d0

            iplan=4
            ibpecorn(iplan,nmagmax)=5
            icorn=1
            bpemag0(1,icorn,iplan,nmagmax)=-xlen/2.0d0
            bpemag0(2,icorn,iplan,nmagmax)=+ylen/2.0d0
            bpemag0(3,icorn,iplan,nmagmax)=-zlen/2.0d0
            icorn=2
            bpemag0(1,icorn,iplan,nmagmax)=-xlen/2.0d0
            bpemag0(2,icorn,iplan,nmagmax)=+ylen/2.0d0
            bpemag0(3,icorn,iplan,nmagmax)=+zlen/2.0d0
            icorn=3
            bpemag0(1,icorn,iplan,nmagmax)=+xlen/2.0d0
            bpemag0(2,icorn,iplan,nmagmax)=+ylen/2.0d0
            bpemag0(3,icorn,iplan,nmagmax)=+zlen/2.0d0
            icorn=4
            bpemag0(1,icorn,iplan,nmagmax)=+xlen/2.0d0
            bpemag0(2,icorn,iplan,nmagmax)=+ylen/2.0d0
            bpemag0(3,icorn,iplan,nmagmax)=-zlen/2.0d0

            iplan=5
            ibpecorn(iplan,nmagmax)=5
            icorn=1
            bpemag0(1,icorn,iplan,nmagmax)=-xlen/2.0d0
            bpemag0(2,icorn,iplan,nmagmax)=-ylen/2.0d0
            bpemag0(3,icorn,iplan,nmagmax)=-zlen/2.0d0
            icorn=2
            bpemag0(1,icorn,iplan,nmagmax)=-xlen/2.0d0
            bpemag0(2,icorn,iplan,nmagmax)=+ylen/2.0d0
            bpemag0(3,icorn,iplan,nmagmax)=-zlen/2.0d0
            icorn=3
            bpemag0(1,icorn,iplan,nmagmax)=+xlen/2.0d0
            bpemag0(2,icorn,iplan,nmagmax)=+ylen/2.0d0
            bpemag0(3,icorn,iplan,nmagmax)=-zlen/2.0d0
            icorn=4
            bpemag0(1,icorn,iplan,nmagmax)=+xlen/2.0d0
            bpemag0(2,icorn,iplan,nmagmax)=-ylen/2.0d0
            bpemag0(3,icorn,iplan,nmagmax)=-zlen/2.0d0

            iplan=6
            ibpecorn(iplan,nmagmax)=5
            icorn=4
            bpemag0(1,icorn,iplan,nmagmax)=-xlen/2.0d0
            bpemag0(2,icorn,iplan,nmagmax)=-ylen/2.0d0
            bpemag0(3,icorn,iplan,nmagmax)=+zlen/2.0d0
            icorn=3
            bpemag0(1,icorn,iplan,nmagmax)=-xlen/2.0d0
            bpemag0(2,icorn,iplan,nmagmax)=+ylen/2.0d0
            bpemag0(3,icorn,iplan,nmagmax)=+zlen/2.0d0
            icorn=2
            bpemag0(1,icorn,iplan,nmagmax)=+xlen/2.0d0
            bpemag0(2,icorn,iplan,nmagmax)=+ylen/2.0d0
            bpemag0(3,icorn,iplan,nmagmax)=+zlen/2.0d0
            icorn=1
            bpemag0(1,icorn,iplan,nmagmax)=+xlen/2.0d0
            bpemag0(2,icorn,iplan,nmagmax)=-ylen/2.0d0
            bpemag0(3,icorn,iplan,nmagmax)=+zlen/2.0d0

            do iplan=1,6
              bpemag0(1,5,iplan,nmagmax)=bpemag0(1,1,iplan,nmagmax)
              bpemag0(2,5,iplan,nmagmax)=bpemag0(2,1,iplan,nmagmax)
              bpemag0(3,5,iplan,nmagmax)=bpemag0(3,1,iplan,nmagmax)
            enddo !iplan

          else if (nplan.eq.-1) then

            call util_skip_comment(lunbpe)
            read(lunbpe,*)xlen,ylen,zlen,ndiv

            nmag=nmag+n2div-1
            kmag2=kmag2+n2div-1

            rmin=xlen
            rmax=ylen

            ylenmax=2.0d0*rmax/ndiv
            ylenmin=2.0d0*rmin/ndiv

            xm=bpebc0(4,nmagmax)
            ym=bpebc0(5,nmagmax)
            zm=bpebc0(6,nmagmax)

            x0=bpebc0(1,nmagmax)
            ycen=bpebc0(2,nmagmax)
            y0max=ycen-rmax-ylenmax/2.0d0
            y0min=ycen-rmin-ylenmin/2.0d0
            z0=bpebc0(3,nmagmax)

            nmagmax=nmagmax-1

            do idiv=1,n2div

c divide cylindrical magnet into rectangular blocks

              nmagmax=nmagmax+1

              if ((idiv/2)*2.ne.idiv.or.rmin.le.0.0d0)  then
                y0max=y0max+ylenmax
                xlen=2.0d0*
     &            sqrt((rmax+(y0max-ycen))*(rmax-(y0max-ycen)))
                ylen=ylenmax
            bpebc0(2,nmagmax)=y0max
                bpebc0(4,nmagmax)=xm
                bpebc0(5,nmagmax)=ym
                bpebc0(6,nmagmax)=zm
              else
                y0MIN=y0MIN+ylenMIN
                xlen=2.0d0*
     &           sqrt((rMIN+(y0MIN-ycen))*(rMIN-(y0MIN-ycen)))
                ylen=ylenmin
            bpebc0(2,nmagmax)=y0min
                bpebc0(4,nmagmax)=-xm
                bpebc0(5,nmagmax)=-ym
                bpebc0(6,nmagmax)=-zm
              endif

              bpebc0(1,nmagmax)=x0
              bpebc0(3,nmagmax)=z0

              bpebc0(7,nmagmax)=-6
              nplan=6

              ibpeplan(nmagmax)=6

              iplan=1
              ibpecorn(iplan,nmagmax)=5
              icorn=1
              bpemag0(1,icorn,iplan,nmagmax)=-xlen/2.0d0
              bpemag0(2,icorn,iplan,nmagmax)=-ylen/2.0d0
              bpemag0(3,icorn,iplan,nmagmax)=-zlen/2.0d0
              icorn=2
              bpemag0(1,icorn,iplan,nmagmax)=-xlen/2.0d0
              bpemag0(2,icorn,iplan,nmagmax)=-ylen/2.0d0
              bpemag0(3,icorn,iplan,nmagmax)=+zlen/2.0d0
              icorn=3
              bpemag0(1,icorn,iplan,nmagmax)=-xlen/2.0d0
              bpemag0(2,icorn,iplan,nmagmax)=+ylen/2.0d0
              bpemag0(3,icorn,iplan,nmagmax)=+zlen/2.0d0
              icorn=4
              bpemag0(1,icorn,iplan,nmagmax)=-xlen/2.0d0
              bpemag0(2,icorn,iplan,nmagmax)=+ylen/2.0d0
              bpemag0(3,icorn,iplan,nmagmax)=-zlen/2.0d0

              iplan=2
              ibpecorn(iplan,nmagmax)=5
              icorn=4
              bpemag0(1,icorn,iplan,nmagmax)=+xlen/2.0d0
              bpemag0(2,icorn,iplan,nmagmax)=-ylen/2.0d0
              bpemag0(3,icorn,iplan,nmagmax)=-zlen/2.0d0
              icorn=3
              bpemag0(1,icorn,iplan,nmagmax)=+xlen/2.0d0
              bpemag0(2,icorn,iplan,nmagmax)=-ylen/2.0d0
              bpemag0(3,icorn,iplan,nmagmax)=+zlen/2.0d0
              icorn=2
              bpemag0(1,icorn,iplan,nmagmax)=+xlen/2.0d0
              bpemag0(2,icorn,iplan,nmagmax)=+ylen/2.0d0
              bpemag0(3,icorn,iplan,nmagmax)=+zlen/2.0d0
              icorn=1
              bpemag0(1,icorn,iplan,nmagmax)=+xlen/2.0d0
              bpemag0(2,icorn,iplan,nmagmax)=+ylen/2.0d0
              bpemag0(3,icorn,iplan,nmagmax)=-zlen/2.0d0

              iplan=3
              ibpecorn(iplan,nmagmax)=5
              icorn=4
              bpemag0(1,icorn,iplan,nmagmax)=-xlen/2.0d0
              bpemag0(2,icorn,iplan,nmagmax)=-ylen/2.0d0
              bpemag0(3,icorn,iplan,nmagmax)=-zlen/2.0d0
              icorn=3
              bpemag0(1,icorn,iplan,nmagmax)=-xlen/2.0d0
              bpemag0(2,icorn,iplan,nmagmax)=-ylen/2.0d0
              bpemag0(3,icorn,iplan,nmagmax)=+zlen/2.0d0
              icorn=2
              bpemag0(1,icorn,iplan,nmagmax)=+xlen/2.0d0
              bpemag0(2,icorn,iplan,nmagmax)=-ylen/2.0d0
              bpemag0(3,icorn,iplan,nmagmax)=+zlen/2.0d0
              icorn=1
              bpemag0(1,icorn,iplan,nmagmax)=+xlen/2.0d0
              bpemag0(2,icorn,iplan,nmagmax)=-ylen/2.0d0
              bpemag0(3,icorn,iplan,nmagmax)=-zlen/2.0d0

              iplan=4
              ibpecorn(iplan,nmagmax)=5
              icorn=1
              bpemag0(1,icorn,iplan,nmagmax)=-xlen/2.0d0
              bpemag0(2,icorn,iplan,nmagmax)=+ylen/2.0d0
              bpemag0(3,icorn,iplan,nmagmax)=-zlen/2.0d0
              icorn=2
              bpemag0(1,icorn,iplan,nmagmax)=-xlen/2.0d0
              bpemag0(2,icorn,iplan,nmagmax)=+ylen/2.0d0
              bpemag0(3,icorn,iplan,nmagmax)=+zlen/2.0d0
              icorn=3
              bpemag0(1,icorn,iplan,nmagmax)=+xlen/2.0d0
              bpemag0(2,icorn,iplan,nmagmax)=+ylen/2.0d0
              bpemag0(3,icorn,iplan,nmagmax)=+zlen/2.0d0
              icorn=4
              bpemag0(1,icorn,iplan,nmagmax)=+xlen/2.0d0
              bpemag0(2,icorn,iplan,nmagmax)=+ylen/2.0d0
              bpemag0(3,icorn,iplan,nmagmax)=-zlen/2.0d0

              iplan=5
              ibpecorn(iplan,nmagmax)=5
              icorn=1
              bpemag0(1,icorn,iplan,nmagmax)=-xlen/2.0d0
              bpemag0(2,icorn,iplan,nmagmax)=-ylen/2.0d0
              bpemag0(3,icorn,iplan,nmagmax)=-zlen/2.0d0
              icorn=2
              bpemag0(1,icorn,iplan,nmagmax)=-xlen/2.0d0
              bpemag0(2,icorn,iplan,nmagmax)=+ylen/2.0d0
              bpemag0(3,icorn,iplan,nmagmax)=-zlen/2.0d0
              icorn=3
              bpemag0(1,icorn,iplan,nmagmax)=+xlen/2.0d0
              bpemag0(2,icorn,iplan,nmagmax)=+ylen/2.0d0
              bpemag0(3,icorn,iplan,nmagmax)=-zlen/2.0d0
              icorn=4
              bpemag0(1,icorn,iplan,nmagmax)=+xlen/2.0d0
              bpemag0(2,icorn,iplan,nmagmax)=-ylen/2.0d0
              bpemag0(3,icorn,iplan,nmagmax)=-zlen/2.0d0

              iplan=6
              ibpecorn(iplan,nmagmax)=5
              icorn=4
              bpemag0(1,icorn,iplan,nmagmax)=-xlen/2.0d0
              bpemag0(2,icorn,iplan,nmagmax)=-ylen/2.0d0
              bpemag0(3,icorn,iplan,nmagmax)=+zlen/2.0d0
              icorn=3
              bpemag0(1,icorn,iplan,nmagmax)=-xlen/2.0d0
              bpemag0(2,icorn,iplan,nmagmax)=+ylen/2.0d0
              bpemag0(3,icorn,iplan,nmagmax)=+zlen/2.0d0
              icorn=2
              bpemag0(1,icorn,iplan,nmagmax)=+xlen/2.0d0
              bpemag0(2,icorn,iplan,nmagmax)=+ylen/2.0d0
              bpemag0(3,icorn,iplan,nmagmax)=+zlen/2.0d0
              icorn=1
              bpemag0(1,icorn,iplan,nmagmax)=+xlen/2.0d0
              bpemag0(2,icorn,iplan,nmagmax)=-ylen/2.0d0
              bpemag0(3,icorn,iplan,nmagmax)=+zlen/2.0d0

              do iplan=1,6
                bpemag0(1,5,iplan,nmagmax)=bpemag0(1,1,iplan,nmagmax)
                bpemag0(2,5,iplan,nmagmax)=bpemag0(2,1,iplan,nmagmax)
                bpemag0(3,5,iplan,nmagmax)=bpemag0(3,1,iplan,nmagmax)
              enddo !iplan

            enddo !ndiv

          endif !nplan.gt.0

        enddo !imag

        lunbpe=lunpm

        call util_skip_comment(lunbpe)
        read(lunbpe,*)nmodule

        nmagmax=nmagmax-nmag

        do imodul=1,nmodule

          call util_skip_comment(lunbpe)
          read(lunbpe,*)xmod,ymod,zmod

          call util_skip_comment(lunbpe)
          read(lunbpe,*)rotmod(1,1),rotmod(1,2),rotmod(1,3)

          call util_skip_comment(lunbpe)
          read(lunbpe,*)rotmod(2,1),rotmod(2,2),rotmod(2,3)

          call util_skip_comment(lunbpe)
          read(lunbpe,*)rotmod(3,1),rotmod(3,2),rotmod(3,3)

          call util_determinante(3,rotmod,det,ifail)
          if (ifail.ne.0.or.abs(abs(det)-1.0d0).gt.tiny) then
            print *,
     &        '*** Error in bpolyini: Bad rotation matrix'
            print *,'magnet, plane: ',imag,iplan
            stop
          endif

          call util_skip_comment(lunbpe)
          read(lunbpe,*)ncopy

          call util_skip_comment(lunbpe)
          read(lunbpe,*)space,vspace(1),vspace(2),vspace(3)

          call util_vnorm(3,vspace,vspace)

          call util_skip_comment(lunbpe)
          read(lunbpe,*)vbsym(1),vbsym(2),vbsym(3)

          do icopy=1,ncopy

            do imag=kmag1,kmag2

              nmagmax=nmagmax+1

              ibpecol(nmagmax)=ibpecol(imag)
              ibpeplan(nmagmax)=ibpeplan(imag)

              x0=xmod+bpebc0(1,imag)+(icopy-1)*space*vspace(1)
              y0=ymod+bpebc0(2,imag)+(icopy-1)*space*vspace(2)
              z0=zmod+bpebc0(3,imag)+(icopy-1)*space*vspace(3)

              bpebc(1,nmagmax)=x0
              bpebc(2,nmagmax)=y0
              bpebc(3,nmagmax)=z0

              bpebc(4,nmagmax)=bpebc0(4,imag)*vbsym(1)
              bpebc(5,nmagmax)=bpebc0(5,imag)*vbsym(2)
              bpebc(6,nmagmax)=bpebc0(6,imag)*vbsym(3)

              bpebc(8,nmagmax)=bpebc0(7,imag)

              do iplan=1,ibpeplan(nmagmax)

                ibpecorn(iplan,nmagmax)=ibpecorn(iplan,imag)

                do icorn=1,ibpecorn(iplan,nmagmax)

                  bpemag(1,icorn,iplan,nmagmax)=
     &              bpemag0(1,icorn,iplan,imag)
                  bpemag(2,icorn,iplan,nmagmax)=
     &              bpemag0(2,icorn,iplan,imag)
                  bpemag(3,icorn,iplan,nmagmax)=
     &              bpemag0(3,icorn,iplan,imag)

c rotation of magnet
                  x0=
     &              rotmod(1,1)*bpemag(1,icorn,iplan,nmagmax)+
     &              rotmod(1,2)*bpemag(2,icorn,iplan,nmagmax)+
     &              rotmod(1,3)*bpemag(3,icorn,iplan,nmagmax)

                  y0=
     &              rotmod(2,1)*bpemag(1,icorn,iplan,nmagmax)+
     &              rotmod(2,2)*bpemag(2,icorn,iplan,nmagmax)+
     &              rotmod(2,3)*bpemag(3,icorn,iplan,nmagmax)

                  z0=
     &              rotmod(3,1)*bpemag(1,icorn,iplan,nmagmax)+
     &              rotmod(3,2)*bpemag(2,icorn,iplan,nmagmax)+
     &              rotmod(3,3)*bpemag(3,icorn,iplan,nmagmax)

                  bpemag(1,icorn,iplan,nmagmax)=x0+bpebc(1,nmagmax)
                  bpemag(2,icorn,iplan,nmagmax)=y0+bpebc(2,nmagmax)
                  bpemag(3,icorn,iplan,nmagmax)=z0+bpebc(3,nmagmax)

                enddo !icorn

                if (det.lt.0.0d0) then

                  do icorn=1,ibpecorn(iplan,nmagmax)
                    ip1=ibpecorn(iplan,nmagmax)-icorn+1
                    shuffle(1,ip1,iplan,nmagmax)=bpemag(1,icorn,iplan,nmagmax)
                    shuffle(2,ip1,iplan,nmagmax)=bpemag(2,icorn,iplan,nmagmax)
                    shuffle(3,ip1,iplan,nmagmax)=bpemag(3,icorn,iplan,nmagmax)
                  enddo !icorn

                  do icorn=1,ibpecorn(iplan,nmagmax)
                    bpemag(1,icorn,iplan,nmagmax)=
     &                shuffle(1,icorn,iplan,nmagmax)
                    bpemag(2,icorn,iplan,nmagmax)=
     &                shuffle(2,icorn,iplan,nmagmax)
                    bpemag(3,icorn,iplan,nmagmax)=
     &                shuffle(3,icorn,iplan,nmagmax)
                  enddo !icorn

                endif !det

              enddo !iplan

            enddo !imag

          enddo !icopy

        enddo !imodul=1,nmodul

        goto 1

9     continue !end of input file

      nmag=nmagmax

      close(lunbpe)

      do imag=1,nmag
c center of magnet in lab

c Anders als in POLYMAG, werden hier alle Magnete geschoben

        hgap=gappm/2.0d0

        if (bpebc(2,imag).lt.0.0d0.and.bpebc(3,imag).lt.0.0d0) then
          bpebc(1,imag)=bpebc(1,imag)+shiftll
          bpebc(2,imag)=bpebc(2,imag)-hgap
          do iplan=1,ibpeplan(imag)
            do icorn=1,ibpecorn(iplan,imag)
              bpemag(1,icorn,iplan,imag)=bpemag(1,icorn,iplan,imag)+shiftll
              bpemag(2,icorn,iplan,imag)=bpemag(2,icorn,iplan,imag)-hgap
            enddo
          enddo
        else if (bpebc(2,imag).lt.0.0d0.and.bpebc(3,imag).gt.0.0d0) then
          bpebc(1,imag)=bpebc(1,imag)+shiftlr
          bpebc(2,imag)=bpebc(2,imag)-hgap
          do iplan=1,6
            do icorn=1,5
              bpemag(1,icorn,iplan,imag)=bpemag(1,icorn,iplan,imag)+shiftlr
              bpemag(2,icorn,iplan,imag)=bpemag(2,icorn,iplan,imag)-hgap
            enddo
          enddo
        else if (bpebc(2,imag).gt.0.0d0.and.bpebc(3,imag).gt.0.0d0) then
          bpebc(1,imag)=bpebc(1,imag)+shiftur
          bpebc(2,imag)=bpebc(2,imag)+hgap
          do iplan=1,6
            do icorn=1,5
              bpemag(1,icorn,iplan,imag)=bpemag(1,icorn,iplan,imag)+shiftur
              bpemag(2,icorn,iplan,imag)=bpemag(2,icorn,iplan,imag)+hgap
            enddo
          enddo
        else if (bpebc(2,imag).gt.0.0d0.and.bpebc(3,imag).lt.0.0d0) then
          bpebc(1,imag)=bpebc(1,imag)+shiftul
          bpebc(2,imag)=bpebc(2,imag)+hgap
          do iplan=1,6
            do icorn=1,5
              bpemag(1,icorn,iplan,imag)=bpemag(1,icorn,iplan,imag)+shiftul
              bpemag(2,icorn,iplan,imag)=bpemag(2,icorn,iplan,imag)+hgap
            enddo
          enddo
        endif

        rmag(1)=bpebc(1,imag)
        rmag(2)=bpebc(2,imag)
        rmag(3)=bpebc(3,imag)

c magnetization vector in lab

        vmaglab(1)=bpebc(4,imag)
        vmaglab(2)=bpebc(5,imag)
        vmaglab(3)=bpebc(6,imag)

        bc=sqrt(vmaglab(1)**2+vmaglab(2)**2+vmaglab(3)**2)

        do iplan=1,ibpeplan(imag)

c three points defining plane (lab.-system)

          p1(1)=bpemag(1,1,iplan,imag)
          p1(2)=bpemag(2,1,iplan,imag)
          p1(3)=bpemag(3,1,iplan,imag)

          p2(1)=bpemag(1,2,iplan,imag)
          p2(2)=bpemag(2,2,iplan,imag)
          p2(3)=bpemag(3,2,iplan,imag)

          p3(1)=bpemag(1,3,iplan,imag)
          p3(2)=bpemag(2,3,iplan,imag)
          p3(3)=bpemag(3,3,iplan,imag)

          call bpen(imag,iplan,p1,p2,p3,vnormlab)

c check if normal vector is perpendicular to magnetization vector
c if mag. vector is parallel, skip plane

          if (bc.ne.0.d0) then
            dum=abs(
     &        (vnormlab(1)*vmaglab(1)+vnormlab(2)*vmaglab(2)+
     &        vnormlab(3)*vmaglab(3))
     &        /bc
     &        )
          else
            dum=0.0d0
          endif

          if (dum.lt.1.d-20.and.bpebc(8,imag).ne.-6) then

            ibpecorn(iplan,imag)=-ibpecorn(iplan,imag)
            bpetm(1,8,iplan,imag)=vnormlab(1)
            bpetm(2,8,iplan,imag)=vnormlab(2)
            bpetm(3,8,iplan,imag)=vnormlab(3)

          else

            bpetm(1,7,iplan,imag)=
     &        vmaglab(1)*vnormlab(1)+
     &        vmaglab(2)*vnormlab(2)+
     &        vmaglab(3)*vnormlab(3)
            bpetm(1,8,iplan,imag)=vnormlab(1)
            bpetm(2,8,iplan,imag)=vnormlab(2)
            bpetm(3,8,iplan,imag)=vnormlab(3)

c get matrices ts and tsinv

            call bpet(vnormlab,ts,tsinv)

            if (bpebc(8,imag).eq.-6) then

c for rectangular magnets, we rotate plan such, that the flanges coinside with
c the axis of the coord.-system.


c All planes are rotated to the system of the
c first plane

              if (iplan.eq.1) then
                ts1=ts
                ts1inv=tsinv
              else
                ts=ts1
                tsinv=ts1inv
              endif !(iplan.eq.1)

              do icorn=1,2

                r1lab(1)=bpemag(1,icorn,iplan,imag)
                r1lab(2)=bpemag(2,icorn,iplan,imag)
                r1lab(3)=bpemag(3,icorn,iplan,imag)

                r1(1)=ts(1,1)*r1lab(1)+ts(1,2)*r1lab(2)+ts(1,3)*r1lab(3)
                r1(2)=ts(2,1)*r1lab(1)+ts(2,2)*r1lab(2)+ts(2,3)*r1lab(3)
                r1(3)=ts(3,1)*r1lab(1)+ts(3,2)*r1lab(2)+ts(3,3)*r1lab(3)

                bperot(1,icorn,iplan,imag)=r1(1)
                bperot(2,icorn,iplan,imag)=r1(2)
                bperot(3,icorn,iplan,imag)=r1(3)

              enddo !icorn=1,ncorn

              vx=bperot(1,2,iplan,imag)-bperot(1,1,iplan,imag)
              vy=bperot(2,2,iplan,imag)-bperot(2,1,iplan,imag)
              vn=sqrt(vx*vx+vy*vy)

              sa=vy/vn
              ca=vx/vn

              tz=ts

              ts(1,1)=ca
              ts(1,2)=sa
              ts(1,3)=0.0d0

              ts(2,1)=-sa
              ts(2,2)=ca
              ts(2,3)=0.0d0

              ts(3,1)=0.0d0
              ts(3,2)=0.0d0
              ts(3,3)=1.0d0

              call util_matrix_multiplication(3,3,3,ts,tz,ts,ws)

              do i=1,3
                do j=1,3
                  tsinv(i,j)=ts(j,i)
                enddo
              enddo

            endif ! if (bpebc(8,imag).eq.-6)

            do i=1,3
              do j=1,3
                bpetm(i,j,iplan,imag)=ts(i,j)
                bpetm(i,j+3,iplan,imag)=tsinv(i,j)
              enddo
            enddo

          endif !check if normal vector is perpendicular to magnetization vector

        enddo !iplan

      enddo !imag

      do imag=1,nmag

c check, if all flanges appear twice, i.e. volume is closed

        iwarn=0

        nflange=0
        do iplan=1,ibpeplan(imag)
          do icorn=1,iabs(ibpecorn(iplan,imag))-1
            nflange=nflange+1
            ip1=icorn
            ip2=ip1+1
            bflange(1,nflange)=bpemag(1,ip1,iplan,imag)
            bflange(2,nflange)=bpemag(2,ip1,iplan,imag)
            bflange(3,nflange)=bpemag(3,ip1,iplan,imag)
            bflange(4,nflange)=bpemag(1,ip2,iplan,imag)
            bflange(5,nflange)=bpemag(2,ip2,iplan,imag)
            bflange(6,nflange)=bpemag(3,ip2,iplan,imag)
          enddo ! icorn
        enddo !iplan

        do iflange=1,nflange
          bflange(7,iflange)=1.d0
        enddo

        do iflange=1,nflange

          do i=iflange+1,nflange

            if (
     &          bflange(1,i).eq.bflange(1,iflange) .and.
     &          bflange(2,i).eq.bflange(2,iflange) .and.
     &          bflange(3,i).eq.bflange(3,iflange) .and.
     &          bflange(4,i).eq.bflange(4,iflange) .and.
     &          bflange(5,i).eq.bflange(5,iflange) .and.
     &          bflange(6,i).eq.bflange(6,iflange)
     &          .or.
     &          bflange(4,i).eq.bflange(1,iflange) .and.
     &          bflange(5,i).eq.bflange(2,iflange) .and.
     &          bflange(6,i).eq.bflange(3,iflange) .and.
     &          bflange(1,i).eq.bflange(4,iflange) .and.
     &          bflange(2,i).eq.bflange(5,iflange) .and.
     &          bflange(3,i).eq.bflange(6,iflange)
     &          ) then

              bflange(7,iflange)=bflange(7,iflange)+1.d0
              bflange(7,i)=bflange(7,i)+1.d0

            endif !hit

          enddo !i

          if (bflange(7,iflange).ne.2.d0.and.iwarn.eq.0) then
            iwarn=1
            print *,
c     &        '*** Error in bpolyini: Magnet is not a closed volume'
     &        '*** Warning in bpolyini: Magnet is not a closed volume'
            print *,'magnet: ',imag
            print *
c            stop
          endif

        enddo !iflange

c center of gravity is a point inside the magnet since shape is convex

        x0=0.d0
        y0=0.d0
        z0=0.d0

        i=0
        do iplan=1,ibpeplan(imag)
          do icorn=1,iabs(ibpecorn(iplan,imag))-1
            i=i+1
            x0=x0+bpemag(1,icorn,iplan,imag)
            y0=y0+bpemag(2,icorn,iplan,imag)
            z0=z0+bpemag(3,icorn,iplan,imag)
          enddo ! icorn
        enddo !iplan

        x0=x0/i
        y0=y0/i
        z0=z0/i

        bpexpos(imag)=x0

        do iplan=1,ibpeplan(imag)

          vnormlab(1)=bpetm(1,8,iplan,imag)
          vnormlab(2)=bpetm(2,8,iplan,imag)
          vnormlab(3)=bpetm(3,8,iplan,imag)

c does normal vector point outside?

          vsx=bpemag(1,1,iplan,imag)-x0
          vsy=bpemag(2,1,iplan,imag)-y0
          vsz=bpemag(3,1,iplan,imag)-z0

          if ( vsx*vnormlab(1) + vsy*vnormlab(2) + vsz*vnormlab(3)
     &        .lt. 0.d0 ) then
            print *
            print *,
     &        '*** Error in bpolyini: Normal vector is not pointing outside'
            print *,'magnet, plane: ',imag,iplan
            print *
            stop
          endif

          do icorn=3,iabs(ibpecorn(iplan,imag))-1

            ip1=icorn-2
            ip2=icorn-1

            v1x=bpemag(1,ip2,iplan,imag)-bpemag(1,ip1,iplan,imag)
            v1y=bpemag(2,ip2,iplan,imag)-bpemag(2,ip1,iplan,imag)
            v1z=bpemag(3,ip2,iplan,imag)-bpemag(3,ip1,iplan,imag)

            v2x=bpemag(1,icorn,iplan,imag)-bpemag(1,ip2,iplan,imag)
            v2y=bpemag(2,icorn,iplan,imag)-bpemag(2,ip2,iplan,imag)
            v2z=bpemag(3,icorn,iplan,imag)-bpemag(3,ip2,iplan,imag)

            vsx=v1y*v2z-v1z*v2y
            vsy=v1z*v2x-v1x*v2z
            vsz=v1x*v2y-v1y*v2x

            if ( abs(v2x*vnormlab(1)+ v2y*vnormlab(2)+ v2z*vnormlab(3))
     &          .gt.tiny ) then
              print *
              print *,'*** Error in bpolyini: Points not in a plane'
              print *,'magnet, plane, point: ',imag,iplan,icorn
              print *
              stop
            endif

            if ( vsx*vnormlab(1) + vsy*vnormlab(2) + vsz*vnormlab(3)
     &          .lt. 0.d0 ) then
              print *
              print *,'*** Error in bpolyini: Direction of rotation not unique'
              print *,'magnet, plane, point ',imag,iplan,icorn
              print *
              stop
            endif

          enddo !icorn=1,ncorn

        enddo ! iplan=1,nplan

      enddo ! imag=1,nmag

c transform everything to the nz=(0,0,1) system

      do imag=1,nmag

        qsign=0.d0

        do iplan=1,ibpeplan(imag)

          if (ibpecorn(iplan,imag).gt.0) then

            do i=1,3
              do j=1,3
                ts(i,j)=bpetm(i,j,iplan,imag)
              enddo
            enddo

            vnormlab(1)=bpetm(1,8,iplan,imag)
            vnormlab(2)=bpetm(2,8,iplan,imag)
            vnormlab(3)=bpetm(3,8,iplan,imag)

            do icorn=1,ibpecorn(iplan,imag)

              r1lab(1)=bpemag(1,icorn,iplan,imag)
              r1lab(2)=bpemag(2,icorn,iplan,imag)
              r1lab(3)=bpemag(3,icorn,iplan,imag)

              r1(1)=ts(1,1)*r1lab(1)+ts(1,2)*r1lab(2)+ts(1,3)*r1lab(3)
              r1(2)=ts(2,1)*r1lab(1)+ts(2,2)*r1lab(2)+ts(2,3)*r1lab(3)
              r1(3)=ts(3,1)*r1lab(1)+ts(3,2)*r1lab(2)+ts(3,3)*r1lab(3)

              bperot(1,icorn,iplan,imag)=r1(1)
              bperot(2,icorn,iplan,imag)=r1(2)
              bperot(3,icorn,iplan,imag)=r1(3)

            enddo !icorn=1,ncorn

            do icorn=1,iabs(ibpecorn(iplan,imag))-1

              ip2=icorn+1

              r1(1)=bperot(1,icorn,iplan,imag)
              r1(2)=bperot(2,icorn,iplan,imag)
              r1(3)=bperot(3,icorn,iplan,imag)

              r2(1)=bperot(1,ip2,iplan,imag)
              r2(2)=bperot(2,ip2,iplan,imag)
              r2(3)=bperot(3,ip2,iplan,imag)

              if (abs(r1(1)-r2(1)).gt.tiny) then

                a=(r2(2)-r1(2))/(r2(1)-r1(1))
                b=r1(2)-a*r1(1)

              else

                a=0.0d0
                b=r1(2)

              endif !(abs(r1(1)-r2(1)).gt.tiny)

              q=-((a*r1(1)+ a*r2(1) + 2*b)*(r1(1) - r2(1)))/2.d0

              qsign=qsign+q*(
     &           vnormlab(1)*bpebc(4,imag)
     &          +vnormlab(2)*bpebc(5,imag)
     &          +vnormlab(3)*bpebc(6,imag))

            enddo ! icorn

          endif !(ibpecorn(iplan,imag).gt.0) then

        enddo ! iplan=1,nplan

        if (abs(qsign).gt.1.d-10.and.bpebc(8,imag).ne.-6) then
          print *
          print *,
     &      '*** Warning in BPOLYINI: Sum of magnetic charge not zero ***'
          print *,'magnet: ',imag
          print *
        endif

      enddo ! imag=1,nmag

      open(unit=99,file='polymag.mag',form='formatted',status='unknown')

      xmn=1.0d30
      xmx=-1.0d30

      do imag=1,nmag
        do iplan=1,ibpeplan(imag)
          do icorn=1,abs(ibpecorn(iplan,imag))

            if (xmn.gt.bpemag(1,icorn,iplan,imag))
     &          xmn=bpemag(1,icorn,iplan,imag)
            if (xmx.lt.bpemag(1,icorn,iplan,imag))
     &          xmx=bpemag(1,icorn,iplan,imag)

            htup(1)=imag
            htup(2)=ibpecol(imag)
            htup(3)=iplan
            htup(4)=icorn*sign(1,ibpecorn(iplan,imag))
            htup(5)=bpemag(1,icorn,iplan,imag)
            htup(6)=bpemag(2,icorn,iplan,imag)
            htup(7)=bpemag(3,icorn,iplan,imag)
            htup(8)=bpebc(4,imag)
            htup(9)=bpebc(5,imag)
            htup(10)=bpebc(6,imag)

            vmaglab(1)=bpebc(4,imag)
            vmaglab(2)=bpebc(5,imag)
            vmaglab(3)=bpebc(6,imag)

            bc=sqrt(vmaglab(1)**2+vmaglab(2)**2+vmaglab(3)**2)
            bpebc(7,imag)=bc

            if (bc.ne.0.d0) then
              if (xmn.gt.bpemag(1,icorn,iplan,imag))
     &          xmn=bpemag(1,icorn,iplan,imag)
              if (xmx.lt.bpemag(1,icorn,iplan,imag))
     &          xmx=bpemag(1,icorn,iplan,imag)
              write(99,'(4f7.0,6e15.5)')htup
            endif

      enddo !ncorn
      enddo !nplan
      enddo !nmag

      close(99)

      if (xstart.eq.9999.d0) xstart=xmn/1000.0d0-rangpm
      if (xstop .eq.9999.d0) xstop =xmx/1000.0d0+rangpm

      open(unit=34,file='polymag.err',status='unknown',form='formatted')

      write(lungfo,*)
      write(lungfo,*)'      SUBROUTINE BPOLYINI:'
      write(lungfo,*)
      write(lungfo,*)'      Comment on file polymag.in:'
      write(lungfo,*)'      ',usercom
      write(lungfo,*)'      WINPM, RANGPM: ',WINPM,RANGPM
      write(lungfo,*)'      BSCALEPM: ',BSCALEPM
      write(lungfo,*)
      write(lungfo,*)'      SHIFTLL,SHIFTLR:',SHIFTLL,SHIFTLR
      write(lungfo,*)'      SHIFTUL,SHIFTUR:',SHIFTUL,SHIFTUR
      write(lungfo,*)'      GAPPM:          ',GAPPM
      write(lungfo,*)

c sort by x-position and weed zero-magnets {

c store values
      do imag=1,nmagmax

        bpexpos0(imag)=bpexpos(imag)

        do i=1,8
          bpebc0(i,imag)=bpebc(i,imag)
        enddo

        ibpeplan0(imag)=ibpeplan(imag)

        do i=1,3
          do j=1,3
            bpemat0(j,i,imag)=bpemat(j,i,imag)
          enddo
        enddo

        do iplan=1,ibpeplan(imag)

          ibpecorn0(iplan,imag)=ibpecorn(iplan,imag)

          do icorn=1,abs(ibpecorn(iplan,imag))

            do i=1,6
              bpemag0(i,icorn,iplan,imag)=bpemag(i,icorn,iplan,imag)
            enddo

            do i=1,3
              bperot0(i,icorn,iplan,imag)=bperot(i,icorn,iplan,imag)
              do j=1,8
                bpetm0(i,j,iplan,imag)=bpetm(i,j,iplan,imag)
              enddo
            enddo

          enddo !icorn
        enddo !iplan
      enddo !imag

c weed

      imag=0
      do kmag=1,nmagmax

        if (bpebc0(7,kmag).eq.0.0d0) goto 999

      imag=imag+1

        bpexpos(imag)=bpexpos0(kmag)

        do i=1,8
          bpebc(i,imag)=bpebc0(i,kmag)
        enddo

        ibpeplan(imag)=ibpeplan0(kmag)

        do i=1,3
          do j=1,3
            bpemat(j,i,imag)=bpemat0(j,i,kmag)
          enddo
        enddo

        do iplan=1,ibpeplan0(kmag)

          ibpecorn(iplan,imag)=ibpecorn0(iplan,kmag)

          do icorn=1,abs(ibpecorn0(iplan,kmag))

            do i=1,6
              bpemag(i,icorn,iplan,imag)=bpemag0(i,icorn,iplan,kmag)
            enddo

            do i=1,3
              bperot(i,icorn,iplan,imag)=bperot0(i,icorn,iplan,kmag)
              do j=1,8
                bpetm(i,j,iplan,imag)=bpetm0(i,j,iplan,kmag)
              enddo
            enddo

          enddo
        enddo
999   continue
      enddo !imag

      nmag=imag
      nmagmax=nmag

      allocate(xsort(nmagmax))
      allocate(ysort(nmagmax))

c store values
      do imag=1,nmagmax

        do i=1,8
          bpebc0(i,imag)=bpebc(i,imag)
        enddo

        ibpeplan0(imag)=ibpeplan(imag)

        do i=1,3
          do j=1,3
            bpemat0(j,i,imag)=bpemat(j,i,imag)
          enddo
        enddo

        do iplan=1,ibpeplan(imag)

          ibpecorn0(iplan,imag)=ibpecorn(iplan,imag)

          do icorn=1,abs(ibpecorn(iplan,imag))

            do i=1,6
              bpemag0(i,icorn,iplan,imag)=bpemag(i,icorn,iplan,imag)
            enddo

            do i=1,3
              bperot0(i,icorn,iplan,imag)=bperot(i,icorn,iplan,imag)
              do j=1,8
                bpetm0(i,j,iplan,imag)=bpetm(i,j,iplan,imag)
              enddo
            enddo

          enddo !icorn
        enddo !iplan
      enddo !imag

      do imag=1,nmagmax
        bpexpos0(imag)=bpexpos(imag)
        xsort(imag)=bpexpos(imag)
        ysort(imag)=imag
      enddo

      call util_sort_func(nmagmax,xsort,ysort)

      imag=0
      do jmag=1,nmagmax

      imag=imag+1

        kmag=ysort(jmag)

        bpexpos(imag)=bpexpos0(kmag)

        do i=1,8
          bpebc(i,imag)=bpebc0(i,kmag)
        enddo

        ibpeplan(imag)=ibpeplan0(kmag)

        do i=1,3
          do j=1,3
            bpemat(j,i,imag)=bpemat0(j,i,kmag)
          enddo
        enddo

        do iplan=1,ibpeplan0(kmag)

          ibpecorn(iplan,imag)=ibpecorn0(iplan,kmag)

          do icorn=1,abs(ibpecorn0(iplan,kmag))

            do i=1,6
              bpemag(i,icorn,iplan,imag)=bpemag0(i,icorn,iplan,kmag)
            enddo

            do i=1,3
              bperot(i,icorn,iplan,imag)=bperot0(i,icorn,iplan,kmag)
              do j=1,8
                bpetm(i,j,iplan,imag)=bpetm0(i,j,iplan,kmag)
              enddo
            enddo

          enddo
        enddo
      enddo !imag

      deallocate(bpebc0)
      deallocate(bpemat0)
      deallocate(bpemag0)
      deallocate(bperot0)
      deallocate(bpetm0)
      deallocate(ibpeplan0)
      deallocate(ibpecorn0)
      deallocate(xsort)
      deallocate(ysort)
      deallocate(bpexpos0)

c sort by x-position }

      if (iplot.ne.0) then
        call bpolyplot(iplot,xplmin,xplmax,yplmin,yplmax,zplmin,zplmax,
     &    theta,phi,usercom)
      endif !iplot

      if (iaxint.eq.1.or.iaxint.eq.2.or.iaxint.eq.3) then

        xint=x0int*0.001d0
        yint=y0int*0.001d0
        zint=z0int*0.001d0

        vxint=0.0d0
        vyint=0.0d0
        vzint=0.0d0

        if (iaxint.eq.1) then
          vxint=1.0d0
        else if (iaxint.eq.2) then
          vyint=1.0d0
        else if (iaxint.eq.3) then
          vzint=1.0d0
        endif

        if (nstepint.ge.0) then

          bxi=0.0d0
          byi=0.0d0
          bzi=0.0d0

          do imag=1,nmag

            call bpolyint(imag,xint,yint,zint,
     &        vxint,vyint,vzint,
     &        bxint,byint,bzint)

            bxi=bxi+bxint
            byi=byi+byint
            bzi=bzi+bzint

          enddo

          write(lungfo,*)
          write(lungfo,*)'First field integrals [Tm]'
          write(lungfo,*)
     &      '(analytically,only correct for ensemble of rectangular magnets)'
          write(lungfo,*)

          if (iaxint.eq.1) then
            write(lungfo,*)'integration parallel to x-axis, through [mm]:'
            write(lungfo,*)sngl(x0int),sngl(y0int),sngl(z0int)
            write(lungfo,*)
          else if (iaxint.eq.2) then
            write(lungfo,*)'integration parallel to y-axis, through [mm]:'
            write(lungfo,*)sngl(x0int),sngl(y0int),sngl(z0int)
            write(lungfo,*)
          else if (iaxint.eq.3) then
            write(lungfo,*)'integration parallel to z-axis, through [mm]:'
            write(lungfo,*)sngl(x0int),sngl(y0int),sngl(z0int)
            write(lungfo,*)
          endif

          write(lungfo,*)bxi
          write(lungfo,*)byi
          write(lungfo,*)bzi

        endif !(nstepint.ge.0)

        if (nstepint.ne.0) then

          nstepint=abs(nstepint)

          bxi=0.0d0
          byi=0.0d0
          bzi=0.0d0

          do imag=1,nmag

            call bpolyintnum(imag,xint,yint,zint,
     &        vxint,vyint,vzint,
     &        bxint,byint,bzint)

            bxi=bxi+bxint
            byi=byi+byint
            bzi=bzi+bzint

          enddo

          write(lungfo,*)
          write(lungfo,*)'First field integrals [Tm]'
          write(lungfo,*)'(numerically)'
          write(lungfo,*)
          write(lungfo,*)
          write(lungfo,*)bxi
          write(lungfo,*)byi
          write(lungfo,*)bzi
          write(lungfo,*)
          write(lungfo,*)

        endif !(nstepint.ne.0)

      endif ! (iaxint.eq.1.or.iaxint.eq.2.or.iaxint.eq.3)

      ical=1

      else

        return

      endif !ical

      return
      end
