*CMZ :  3.05/00 25/04/2018  13.09.51  by  Michael Scheer
*CMZ :  3.01/02 18/09/2013  12.46.29  by  Michael Scheer
*CMZ :  2.63/05 22/07/2009  08.34.06  by  Michael Scheer
*CMZ :  2.54/00 28/02/2005  17.30.57  by  Michael Scheer
*CMZ :  2.53/05 24/02/2005  14.01.19  by  Michael Scheer
*CMZ :  2.52/05 17/08/2004  12.50.56  by  Michael Scheer
*CMZ :  1.02/01 09/08/2004  14.48.51  by  Michael Scheer
*CMZ :  1.02/00 27/07/2004  13.36.14  by  Michael Scheer
*CMZ :  1.01/01 22/07/2004  13.27.23  by  Michael Scheer
*CMZ :  1.00/00 26/02/2004  17.30.03  by  Michael Scheer
*CMZ :  0.99/07 16/02/2004  16.00.03  by  Michael Scheer
*CMZ :  0.99/03 12/02/2004  10.51.58  by  Michael Scheer
*CMZ :  0.99/01 11/02/2004  13.42.37  by  Michael Scheer
*CMZ :  0.99/00 29/01/2004  12.50.00  by  Michael Scheer
*CMZ :  0.00/08 23/01/2004  15.23.23  by  Michael Scheer
*CMZ :  0.00/07 20/01/2004  16.45.51  by  Michael Scheer
*CMZ :  0.00/06 16/01/2004  10.43.52  by  Michael Scheer
*CMZ :  0.00/05 23/12/2003  16.07.01  by  Michael Scheer
*CMZ :  0.00/04 23/12/2003  10.43.15  by  Michael Scheer
*CMZ :  0.00/02 15/12/2003  12.43.34  by  Michael Scheer
*CMZ :  0.00/01 10/12/2003  17.56.52  by  Michael Scheer
*-- Author :    Michael Scheer   02/12/2003
      subroutine bpolyeder(xin,yin,zin,bxout,byout,bzout)
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
c
c       Calculation of magnetic field of polyhedron according to
c      Oleb Chubar, Pascal Elleaume and Joel Chavanne
c       J. Synchrotron Rad. (1998) 5, 481-484
c
c Paper contains an error: The rotation matrix is wrong, since for nz=-1 the
c                          determinant is -1, which yields to errors??.


*KEEP,bwpolyederf90u.
      include 'bwpolyederf90u.cmn'
*KEND.

      implicit none

      double precision xstart,xstop,win2,xino

      double precision xin,yin,zin,bxout,byout,bzout,
     &  bxm,bym,bzm,bxp,byp,bzp,xxm,yym,zzm,xxp,yyp,zzp,xx0,yy0,zz0,
     &  rr0,rrm

      double precision r1(3),r2(3),dlab(3),blab(3)
      double precision ts(3,3),tsinv(3,3),bplan(3),bcvn,vnormlab(3)
      double precision xx,yy,zz,xxrot,yyrot,zzrot
      double precision a,b,z,qx,qy,qz,pi4inv,reverse,tiny2
      double precision q(3,3),vmagrot(3),vmaglab(3),h(3),
     &  xr(2),yr(2),zr(2),dum,dume

      parameter (pi4inv=0.0795774715459477d0)

      integer ical
      integer iout,kmag1,iwarni,itiny,iwtiny,jtiny
      integer imag,iplan,icorn,i,j,k,ip2,kwarn,kwarntot,lwarn,lwarnmax,mwarn

      data ical/0/
      data iwarni/0/
      data lwarnmax/20/
      data lwarn/0/
      data xino/1.0d30/
      data kmag1/1/

      kwarntot=0
      mwarn=0

      if (ical.eq.0) then
      win2=winpm*500.d0
        tiny2=tiny*tiny
        call bpolyini(xstart,xstop,16)
        ical=1
      endif !ical


c calculate field at (xin,yin,zin)

      xx=xin*1000.0d0
      yy=yin*1000.0d0
      zz=zin*1000.0d0

      xxp=xx
      yyp=yy
      zzp=zz

      xx0=xx
      yy0=yy
      zz0=zz

      xxm=xx
      yym=yy
      zzm=zz

      itiny=0
      jtiny=0
      iwtiny=0

1     continue

c      if (kwarn.gt.1) then
      if (kwarn.gt.0.or.jtiny.ne.0) then !24feb05
        if (kwarn.gt.0) then
          lwarn=lwarn+1
          mwarn=mwarn+1
          if (mwarn.gt.1000) then
            print *,'*** Error in BPOLYEDER: Endless loop?'
            print *,' ical:',ical
            print *,'x,y,z (mm): '
            print *,xin*1000.d0,yin*1000.d0,zin*1000.d0
            print *,'Bx,By,Bz: '
            print *,bxout,byout,bzout
            print *,' '
            stop '***Program WAVE aborted ***'
          else if (mwarn.gt.0.and.ibpnowarn.eq.0) then
            print *,'*** Warning in BPOLYEDER'
            print *,'Counter for this warning:',mwarn
            print *,'x,y,z (mm): '
            print *,xx,yy,zz
            print *,'Bx,By,Bz: '
            print *,bxout,byout,bzout
            print *,'x,y,z will be changed by',1000.*tiny
            print *,' '
          endif
        endif !kwarn
2       if (itiny.eq.0) then
          if (kwarn.eq.0.and.jtiny.eq.1) then
            itiny=1
            iwtiny=0
            goto 2
          else
            xxp=xxp+1000.*tiny
            yyp=yyp+1000.*tiny
            zzp=zzp+1000.*tiny
            xx=xxp
            yy=yyp
            zz=zzp
            jtiny=1
            iwtiny=1
          endif
        else if (itiny.eq.1) then
          if (kwarn.eq.0.and.jtiny.eq.2) then
            itiny=2
            iwtiny=0
            goto 2
          else
            xxm=xxm-1000.*tiny
            yym=yym-1000.*tiny
            zzm=zzm-1000.*tiny
            xx=xxm
            yy=yym
            zz=zzm
            jtiny=2
            iwtiny=2
          endif
        else if (itiny.eq.2) then
          rr0=sqrt((xxp-xxm)**2+(yyp-yym)**2+(zzp-zzm)**2)
          rrm=sqrt((xxm-xx0)**2+(yym-yy0)**2+(zzm-zz0)**2)
          bxout=bxm+(bxp-bxm)/rr0*rrm
          byout=bym+(byp-bym)/rr0*rrm
          bzout=bzm+(bzp-bzm)/rr0*rrm
          return
        endif
      endif !kwarn

      bxout=0.d0
      byout=0.d0
      bzout=0.d0

      if (xin.lt.xino) kmag1=1

      do imag=kmag1,nmag
        if (abs(xx-bpexpos(imag)).le.win2) then
          kmag1=imag
          goto 901
        endif
      enddo

901   xino=xin

      do imag=kmag1,nmag

        if (abs(xx-bpexpos(imag)).gt.win2) goto 999

        if (bpebc(7,imag).ne.0.0d0) then !non-zero magnetization

          if(bpebc(8,imag).eq.1) then !not rectangular magnet

c check, if we are inside of magnet; we assume convex shape

            iout=-1

            do iplan=1,ibpeplan(imag)

              dlab(1)=xx-bpemag(1,1,iplan,imag)
              dlab(2)=yy-bpemag(2,1,iplan,imag)
              dlab(3)=zz-bpemag(3,1,iplan,imag)

              vnormlab(1)=bpetm(1,8,iplan,imag)
              vnormlab(2)=bpetm(2,8,iplan,imag)
              vnormlab(3)=bpetm(3,8,iplan,imag)

              if( dlab(1)*vnormlab(1)+dlab(2)*vnormlab(2)+
     &            dlab(3)*vnormlab(3).gt.0.d0) then
                iout=1
                goto 9
              endif

            enddo !iplan

9           continue

            do iplan=1,ibpeplan(imag)

              bcvn=-bpetm(1,7,iplan,imag)*pi4inv

              bplan(1)=0.d0
              bplan(2)=0.d0
              bplan(3)=0.d0

c transform everything to the nz=(0,0,1) system

              if (ibpecorn(iplan,imag).gt.0) then

                do i=1,3
                  do j=1,3
                    ts(i,j)=bpetm(i,j,iplan,imag)
                    tsinv(i,j)=bpetm(i,j+3,iplan,imag)
                  enddo
                enddo

                xxrot=ts(1,1)*xx+ts(1,2)*yy+ts(1,3)*zz
                yyrot=ts(2,1)*xx+ts(2,2)*yy+ts(2,3)*zz
                zzrot=ts(3,1)*xx+ts(3,2)*yy+ts(3,3)*zz

                do icorn=1,ibpecorn(iplan,imag)-1

                  ip2=icorn+1

                  r1(1)=bperot(1,icorn,iplan,imag)-xxrot
                  r1(2)=bperot(2,icorn,iplan,imag)-yyrot
                  r1(3)=bperot(3,icorn,iplan,imag)-zzrot

                  r2(1)=bperot(1,ip2,iplan,imag)-xxrot
                  r2(2)=bperot(2,ip2,iplan,imag)-yyrot
                  r2(3)=bperot(3,ip2,iplan,imag)-zzrot

                  if (abs(r1(1)-r2(1)).gt.tiny) then

                    a=(r2(2)-r1(2))/(r2(1)-r1(1))
                    b=r1(2)-a*r1(1)

                    if (abs(a).lt.tiny2) then
                      a=0.0d0
                      b=r1(2)
                    endif

                    z=r1(3)

                    call bpeq(r1(1),r2(1),a,b,z,qx,qy,qz,
     &                tiny,reverse,kwarn)

                    if (kwarn.eq.1) then
                      kwarntot=1 !close to magnet boundary, position changed
                      goto 1 !24feb05
                    else if (kwarn.gt.1) then
                      kwarntot=2
                      goto 1
                    endif

                    bplan(1)=bplan(1)-qx*bcvn
                    bplan(2)=bplan(2)-qy*bcvn
                    bplan(3)=bplan(3)-qz*bcvn

                  endif !(abs(r1(1)-r2(1)).gt.tiny)

                enddo !icorn=1,ncorn

                blab(1)=tsinv(1,1)*bplan(1)+tsinv(1,2)*bplan(2)+tsinv(1,3)*bplan(3)
                blab(2)=tsinv(2,1)*bplan(1)+tsinv(2,2)*bplan(2)+tsinv(2,3)*bplan(3)
                blab(3)=tsinv(3,1)*bplan(1)+tsinv(3,2)*bplan(2)+tsinv(3,3)*bplan(3)

                bxout=bxout+blab(1)
                byout=byout+blab(2)
                bzout=bzout+blab(3)

              endif !magnetization parallel to normal vector

            enddo ! iplan=1,nplan

            if (iout.eq.-1) then
              bxout=bxout+bpebc(4,imag)
              byout=byout+bpebc(5,imag)
              bzout=bzout+bpebc(6,imag)
              if (iwarni.eq.0) then
                iwarni=1
                print *
                print *,'*** Error in BPOLYEDER: Inside Magnet! ***'
                print *,'imag,x,y,z (mm): ',imag,sngl(xx),sngl(yy),sngl(zz)
                print*
                stop '*** program WAVE aborted ***'
              endif
            endif !iout

          else !bpebc(8,imag) .eq. 1

c rectangular magnet
c check, if we are inside of magnet; we assume convex shape

            iout=-1

            do iplan=1,ibpeplan(imag)

              dlab(1)=xx-bpemag(1,1,iplan,imag)
              dlab(2)=yy-bpemag(2,1,iplan,imag)
              dlab(3)=zz-bpemag(3,1,iplan,imag)

              vnormlab(1)=bpetm(1,8,iplan,imag)
              vnormlab(2)=bpetm(2,8,iplan,imag)
              vnormlab(3)=bpetm(3,8,iplan,imag)

              if( dlab(1)*vnormlab(1)+dlab(2)*vnormlab(2)+
     &            dlab(3)*vnormlab(3).gt.0.d0) then
                iout=1
                goto 91
              endif

            enddo !iplan

91          continue


            vmaglab(1)=bpebc(4,imag)
            vmaglab(2)=bpebc(5,imag)
            vmaglab(3)=bpebc(6,imag)

c transform everything to the nz=(0,0,1) system and rotate it parallel to x-axis

            do i=1,3
              do j=1,3
                ts(i,j)=bpetm(i,j,1,imag)
                tsinv(i,j)=bpetm(i,j+3,1,imag)
              enddo
            enddo

            xxrot=ts(1,1)*xx+ts(1,2)*yy+ts(1,3)*zz
            yyrot=ts(2,1)*xx+ts(2,2)*yy+ts(2,3)*zz
            zzrot=ts(3,1)*xx+ts(3,2)*yy+ts(3,3)*zz

            vmagrot(1)=
     &        ts(1,1)*vmaglab(1)+ts(1,2)*vmaglab(2)+ts(1,3)*vmaglab(3)
            vmagrot(2)=
     &        ts(2,1)*vmaglab(1)+ts(2,2)*vmaglab(2)+ts(2,3)*vmaglab(3)
            vmagrot(3)=
     &        ts(3,1)*vmaglab(1)+ts(3,2)*vmaglab(2)+ts(3,3)*vmaglab(3)

            xr(1)=bperot(1,1,1,imag)-xxrot
            xr(2)=bperot(1,2,1,imag)-xxrot
            yr(1)=bperot(2,1,1,imag)-yyrot
            yr(2)=bperot(2,3,1,imag)-yyrot

            zr(1)=bperot(3,1,1,imag)-zzrot
            zr(2)=bperot(3,1,3,imag)-zzrot

            if (abs(xr(2)-xr(1)).lt.tiny) then
              print*,'abs(xr(2)-xr(1)).lt.tiny'
            endif

            if (abs(yr(2)-yr(1)).lt.tiny) then
              print*,'abs(yr(2)-yr(1)).lt.tiny'
            endif

            if (abs(zr(2)-zr(1)).lt.tiny) then
              print*,'abs(zr(2)-zr(1)).lt.tiny'
            endif

            kwarn=0

            if (
     &          abs(xr(1)).lt.tiny.or.
     &          abs(xr(2)).lt.tiny.or.
     &          abs(yr(1)).lt.tiny.or.
     &          abs(yr(2)).lt.tiny.or.
     &          abs(zr(1)).lt.tiny.or.
     &          abs(zr(2)).lt.tiny) then
              kwarn=10
              kwarntot=3
              goto 1
            endif

            q(1,1)=0.0d0
            q(2,2)=0.0d0
            q(3,3)=0.0d0

            q(1,2)=1.0d0
            q(1,3)=1.0d0
            q(2,3)=1.0d0

            do i=1,2
              do j=1,2
                do k=1,2

                  q(1,1)=q(1,1)+
     &              (-1)**(i+j+k+1)*
     &              atan(
     &              yr(j)/xr(i)*zr(k)/
     &              sqrt(xr(i)**2+yr(j)**2+zr(k)**2)
     &              )

                  q(2,2)=q(2,2)+
     &              (-1)**(i+j+k+1)*
     &              atan(
     &              xr(j)/yr(i)*zr(k)/
     &              sqrt(yr(i)**2+xr(j)**2+zr(k)**2)
     &              )

                  q(3,3)=q(3,3)+
     &              (-1)**(i+j+k+1)*
     &              atan(
     &              yr(j)/zr(i)*xr(k)/
     &              sqrt(zr(i)**2+yr(j)**2+xr(k)**2)
     &              )

                  dum=zr(k)+sqrt(xr(i)**2+yr(j)**2+zr(k)**2)
                  dume=(-1.0d0)**(i+j+k)

                  if (dum.ne.0.0d0) then
                    if (dume.gt.0.0d0) then
                      q(1,2)=q(1,2)*dum
                    else
                      q(1,2)=q(1,2)/dum
                    endif
                  endif

                  dum=yr(k)+sqrt(xr(i)**2+zr(j)**2+yr(k)**2)
                  if (dum.ne.0.0d0) then
                    if (dume.gt.0.0d0) then
                      q(1,3)=q(1,3)*dum
                    else
                      q(1,3)=q(1,3)/dum
                    endif
                  endif

                  dum=xr(k)+sqrt(zr(i)**2+yr(j)**2+xr(k)**2)
                  if (dum.ne.0.0d0) then
                    if (dume.gt.0.0d0) then
                      q(2,3)=q(2,3)*dum
                    else
                      q(2,3)=q(2,3)/dum
                    endif
                  endif

                enddo !k
              enddo !j
            enddo !i

            if (
     &          q(1,2).lt.tiny2.or.
     &          q(1,3).lt.tiny2.or.
     &          q(2,3).lt.tiny2) then
              kwarn=10
              kwarntot=10
              goto 1
            endif

            q(1,2)=log(q(1,2))
            q(1,3)=log(q(1,3))
            q(2,3)=log(q(2,3))

            q(2,1)=q(1,2)
            q(3,1)=q(1,3)
            q(3,2)=q(2,3)

            h(1)=
     &        -(q(1,1)*vmagrot(1)+q(1,2)*vmagrot(2)+q(1,3)*vmagrot(3))*
     &        pi4inv
            h(2)=
     &        -(q(2,1)*vmagrot(1)+q(2,2)*vmagrot(2)+q(2,3)*vmagrot(3))*
     &        pi4inv
            h(3)=
     &        -(q(3,1)*vmagrot(1)+q(3,2)*vmagrot(2)+q(3,3)*vmagrot(3))*
     &        pi4inv

            bplan(1)=h(1)
            bplan(2)=h(2)
            bplan(3)=h(3)

            blab(1)=tsinv(1,1)*bplan(1)+tsinv(1,2)*bplan(2)+tsinv(1,3)*bplan(3)
            blab(2)=tsinv(2,1)*bplan(1)+tsinv(2,2)*bplan(2)+tsinv(2,3)*bplan(3)
            blab(3)=tsinv(3,1)*bplan(1)+tsinv(3,2)*bplan(2)+tsinv(3,3)*bplan(3)

            bxout=bxout+blab(1)
            byout=byout+blab(2)
            bzout=bzout+blab(3)


            if (iout.eq.-1) then
              bxout=bxout+bpebc(4,imag)
              byout=byout+bpebc(5,imag)
              bzout=bzout+bpebc(6,imag)
              if (iwarni.eq.0) then
                iwarni=1
                print *
                print *,'*** Error in BPOLYEDER: Inside Magnet! ***'
                print *,'imag,x,y,z (mm): ',imag,sngl(xx),sngl(yy),sngl(zz)
                print*
                stop '*** program WAVE aborted ***'
              endif
            endif !iout

          endif !(bpebc(8,imag).eq.1)

        endif !non-zero magnetization

      enddo !imag=1,nmag

999   continue

      bxout=bxout*bscalepm
      byout=byout*bscalepm
      bzout=bzout*bscalepm

      if (jtiny.eq.1) then
        bxp=bxout
        byp=byout
        bzp=bzout
        goto 1
      else if (jtiny.eq.2) then
        bxm=bxout
        bym=byout
        bzm=bzout
        goto 1
      endif

      if (kwarntot.ne.0.and.lwarn.lt.lwarnmax) then
        print *
        print *,'x,y,z (mm): '
        print *,xx,yy,zz
        print *,'Bx,By,Bz: '
        print *,bxout,byout,bzout
        print *,'*** Questionable data written to polymag.err'
        print *,'------------------------------------------'
      endif

      if (lwarn.eq.lwarnmax) then
        lwarn=lwarn+1
        print *
        print *,'*** Maximum number of printed warning reached!!'
        print *,'*** Questionable data written to polymag.err'
        print *,'------------------------------------------'
      endif

      if (kwarntot.ne.0) then
        write(34,'(6e15.5)')
     &    sngl(xin),sngl(yin),sngl(zin),
     &    sngl(bxout),sngl(byout),sngl(bzout)
      endif

      return
      end
