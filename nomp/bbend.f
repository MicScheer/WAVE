*CMZ :  4.00/11 27/07/2021  13.26.18  by  Michael Scheer
*-- Author :    Michael Scheer   01/07/2021
      subroutine bbend(xin,yin,zin,bxout,byout,bzout,axout,ayout,azout,
     &  fint,gap,center,b,pin,vnin,pout,vnout,
     &  modus,istatus,inowarn)

      implicit none

      double precision xin,yin,zin,bxout,byout,bzout,axout,ayout,azout,
     &  fringe,gap,Center(3),B(3),Pin(3),vnin(3),Pout(3),vnout(3),bout(3),bn(3),
     &  bnn,dum(3),r(3),rn(3),rnn,dumn,rin(3),rino(3),
     &  qin(3),qout(3),distin,distout,x,y,z,vnormorb(3),fint

      integer :: modus, istatus, inowarn, ibad=0,istat,iover, ical=0

c The orbit plane contains the points Center,Pin, and Pout. The field vector
c B is normal to the orbit plane.
c The entrace plane is perpendicular the orbit plane, phiin is the
c angle around B, and Pin is a point in the plane;  the exit plane is defind
c accordingly

c The field is calculated with respect to the entrace or exit plane closer to
c xin,yin,zin. Yin is the distance to the orbit plane.

      save rino

      istatus=0
      !ical=ical+1

      axout=0.0d0
      ayout=0.0d0
      azout=0.0d0

      bxout=0.0d0
      byout=0.0d0
      bzout=0.0d0

      bnn=norm2(b)
      if (bnn.eq.0.0d0) goto 9999

      bn=b/bnn

      rin(1)=xin
      rin(2)=yin
      rin(3)=zin

      call util_plane(pout,center,pin,rin,vnormorb,y,iover,istat)

      if (istat.ne.0) then
        if (inowarn.eq.0) then
          print*,"*** Error in bbend: Could not calculate distance to orbit plane  ***"
        endif
        istatus=-2
        return
      endif

      if (norm2(abs(vnormorb)-abs(bn)).gt.1.0d-9) then
        if (inowarn.eq.0) then
          print*,"*** Error in bbend: Normal vector of  orbit plane not parallel to B ***"
        endif
        istatus=-1
        return
      endif

      dumn=norm2(vnin)
      if (dumn.eq.0.0d0) then
        if (inowarn.eq.0) then
          print*,"*** Error in bbend: Zero normal vector of  entrance plane ***"
        endif
        istatus=-1
        return
      else
        vnin=vnin/dumn
      endif

      dumn=norm2(vnout)
      if (dumn.eq.0.0d0) then
        if (inowarn.eq.0) then
          print*,"*** Error in bbend: Zero normal vector of exit plane ***"
        endif
        istatus=-1
        return
      else
        vnout=vnout/dumn
      endif

      if (abs(dot_product(vnin,bn)).gt.1.0d-9) then
        if (inowarn.eq.0) then
          print*,"*** Error in bbend: Normal vector of  entrance plane not perpendicular to B ***"
        endif
        istatus=-1
        return
      endif

      if (abs(dot_product(vnout,bn)).gt.1.0d-9) then
        if (inowarn.eq.0) then
          print*,"*** Error in bbend: Normal vector of  exit plane not perpendicular to B ***"
        endif
        istatus=-1
        return
      endif

      call util_vcross(vnin,bn,qin)
      qin=pin+qin

      call util_plane(pin,pin+bn,qin,rin,vnin,distin,iover,istat)

      if (istat.ne.0) then
        if (inowarn.eq.0) then
          print*,"*** Error in bbend: Could not calculate distance to entrance plane  ***"
        endif
        istatus=-3
        return
      endif

      if (distin.lt.0.0d0) then
        goto 9999 ! outside planes
      endif

      call util_vcross(vnout,bn,qout)
      qout=pout+qout

      call util_plane(pout,pout+bn,qout,rin,vnout,distout,iover,istat)

c      if (ical.ge.4729) then
c        print*, ical,rin,distin,distout
c        stop
c      endif

      if (distout.gt.0.0d0) then
        goto 9999 ! outside planes
      endif

      if (istat.ne.0) then
        if (inowarn.eq.0) then
          print*,"*** Error in bbend: Could not calculate distance to exit plane  ***"
        endif
        istatus=-4
        return
      endif

      if (fint.le.0.0d0.or.gap.le.0.0d0) then

        byout=bnn

      else

        if (abs(distin).le.abs(distout)) then
          x=distin
        else
          x=-distout
        endif !distin

        if (modus.eq.3) then
          call mrad_fringe_cubic_spline(x,y,zin,
     &      bout(1),bout(2),bout(3),axout,ayout,azout,fint,gap,fringe,
     &      istatus)
        else if (modus.eq.5) then
          call mrad_fringe_quintic_spline(x,y,zin,
     &      bout(1),bout(2),bout(3),axout,ayout,azout,fint,gap,fringe,
     &      istatus)
        else
          bout=b
          goto 9999
        endif !modus

        if (istat.ne.0) then
          if (inowarn.eq.0) then
            print*,"*** Error in bbend: Bad return from mrad_fringe_cubic_spline  ***"
          endif
        endif

        if (norm2(bout).eq.0.0d0) goto 9999

        if (abs(distin).le.abs(distout)) then
          bout=(-bout(1)*vnin+bout(2)*bn)*bnn
        else
          bout=(bout(1)*vnout+bout(2)*bn)*bnn
        endif !distin

        bxout=bout(1)
        byout=bout(2)
        bzout=bout(3)

      endif !fringe

c      if (ical.gt.0) then
c        write(88,*) xin,zin,byout,atan2(rin(3)-rino(3),rin(1)-rino(1))*57.29577951308232
c      endif
c      rino=rin
c      ical=ical+1

9999  continue

      return
      end
