*CMZ :          19/01/2025  09.27.21  by  Michael Scheer
*-- Author :    Michael Scheer   23/11/2024
      subroutine sbend(nmag,imag,cbmodel,rho,dbrho,angi,fint,hgap,
     &  cposmodel,xpos,zpos,angex,ebeam,b0,bmovecut,ds,icharge,fringe,fa,fb,fc,istatus)

cdbuff dcbuff
      implicit none

      integer, parameter :: lp=1000

      real*8 r,rho,angle,fint,hgap,posi(7,2),posdum(7,5),ebeam,de,bmovecut,ds,
     &  fringe,fa,fb,fc,b0,a,sina,cosa,d,e1,e2,vx,vz,angi,phi,zenith(3,3),apara(3),zp(3),xopt,zopt,
     &  X2,Y2,Z2,VX2,VY2,VZ2,DTIM,GAMMA,beta,dbrho,da,agoal,v0,anga(2,lp),xe,ze,xi,zi,
     &  y2a(lp),cen(3),ango,arot,vin(3),vout(3),xexit,zexit,dx,dz,xpos,zpos,angex,vrot(3),adjust

      integer nmag,imag,icharge,istatus,i,l
      character(32) cbmodel,cposmodel

*KEEP,mgsqc.
      include 'mgsqc.cmn'
*KEEP,phyconparam.
      include 'phyconparam.cmn'
*KEND.

      call sbend_fringe(cbmodel,fint,hgap,fringe,fa,fb,fc,istatus)
      if (istatus.ne.0) then
        print*,"*** Warning in SBEND: Bad return from SBEND_FRINGE, be careful ***"
      endif

      call sbend_int_fringe(cbmodel,fringe,fa,fb,fc,ebeam,b0,bmovecut,ds,posdum,icharge,istatus)
      if (istatus.ne.0) then
        print*,"*** Warning in SBEND_INT_FRINGE: Bad return from SBEND_FRINGE, be careful ***"
      endif

      angle=angi*grarad1
      gamma=ebeam/emassg1
      beta=dsqrt((1.d0-1.d0/gamma)*(1.d0+1.d0/gamma))
      dtim=ds/(clight1*beta)
      v0=beta*clight1

      r=icharge*rho

      if (rho.lt.0.0d0) then
        agoal=angi*grarad1/2.0d0
      else
        agoal=-angi*grarad1/2.0d0
      endif

      da=abs(atan(posdum(6,2)/posdum(4,2)))

      vin=[cos(agoal),0.0d0,-sin(agoal)]
      adjust=0.0d0

      posi(4,1)=vin(1)
      posi(6,1)=vin(3)
      posi(4,2)=vin(1)
      posi(6,2)=-vin(3)

      xe=abs(r*vin(3))
      ze=r*(1.0d0-vin(1))

      if (nmag.gt.0) then
        pmag(3,imag)=posi(4,1)
        pmag(4,imag)=posi(6,1)
        pmag(7,imag)=posi(4,2)
        pmag(8,imag)=posi(6,2)
      else
        pmag(3,imag)=1.0d0
        pmag(4,imag)=0.0d0
        pmag(7,imag)=1.0d0
        pmag(8,imag)=0.0d0
      endif

      if (cbmodel.eq.'hard-edge') then
        pmag(9,imag)=0.0d0
      else if (cbmodel.eq.'linear') then
        pmag(9,imag)=1.0d0
      else if (cbmodel.eq.'cubic-spline') then
        pmag(9,imag)=3.0d0
      else if (cbmodel.eq.'quintic-spline') then
        pmag(9,imag)=5.0d0
      else
        pmag(9,imag)=0.0d0
        cbmodel="hard-edge"
      endif

      pmag(10,imag)=dbrho/r
c        pmag(11,imag)=edge(1)
c        pmag(12,imag)=edge(2)
      pmag(13,imag)=fint
      pmag(14,imag)=2.0d0*hgap

      pmag(15,imag)=fringe
      pmag(16,imag)=fa
      pmag(17,imag)=fb
      pmag(18,imag)=fc

      do l=1,lp

        if (l.eq.1) then
          angle=angi*grarad1-da
          adjust=-fringe/2.0d0
        else if (l.eq.2) then
          angle=angi*grarad1+da
          adjust=fringe/2.0d0
        endif

        if (rho.lt.0.0d0) then
          a=abs(angle)/2.0d0
        else
          a=-abs(angle)/2.0d0
        endif

        sina=sin(a)
        cosa=cos(a)

        d=adjust

        posi(1,2)=xe+d*cosa
        posi(3,2)=ze+d*sina

        posi(1,1)=-posi(1,2)
        posi(3,1)= posi(3,2)

        pmag(1,imag)=posi(1,1)
        pmag(2,imag)=posi(3,1)

        pmag(5,imag)=posi(1,2)
        pmag(6,imag)=posi(3,2)

        call TRACKBEND(posi(1,1),0.0d0,posi(3,1),v0*vin(1),0.0d0,v0*vin(3),
     &    posi(1,2),posi(2,2),posi(3,2),
     &    vin(1),0.0d0,-vin(3),
     &    X2,Y2,Z2,VX2,VY2,VZ2,zenith,DTIM,GAMMA,iabs(nmag),imag,
     &    bmovecut,0,icharge,0)
        vout=[vx2/v0,0.0d0,vz2/v0]

        if (l.ge.3) then
          anga(1:2,l-1)=anga(1:2,l)
        endif

        anga(1,l)=adjust
        anga(2,l)=abs(acos(dot_product(vin,vout)))
        angle=anga(2,l)

        if (l.ge.2) then
          call util_adjust(angi*grarad1/1.0d0,anga(1,l-1),anga(2,l-1),anga(1,l),anga(2,l))
          adjust=anga(1,l)
        endif

        if (abs((anga(2,l)*radgra1-angi/1.0d0)/(anga(2,l)*radgra1+angi/2.0d0)).le.1.0d-9) then
          exit
        endif

      enddo !Loop to adjust angle

      call util_parabel(
     &  [zenith(1,1),zenith(1,2),zenith(1,3)],
     &  [zenith(3,1),zenith(3,2),zenith(3,3)],
     &  apara,zp,xopt,zopt,istatus)

      if (l.ge.lp) then
        print*,"*** Warning in SBEND: Adjustment of angle failed, be careful ***"
      else
        istatus=0
      endif

      if (angex.ne.0.0d0) then

        cen=0.0d0
        phi=icharge*angex*grarad1
        vrot=[0.0d0,1.0d0,0.0d0]

        vin=[pmag(1,imag),0.0d0,pmag(2,imag)]
        call util_rotate(cen,vrot,phi,vin,vout,istatus)
        pmag(1,imag)=vout(1)
        pmag(2,imag)=vout(3)

        vin=[xopt,0.0d0,zopt]
        call util_rotate(cen,vrot,phi,vin,vout,istatus)
        xopt=vout(1)
        zopt=vout(3)

        vin=[pmag(3,imag),0.0d0,pmag(4,imag)]
        call util_rotate(cen,vrot,phi,vin,vout,istatus)
        pmag(3,imag)=vout(1)
        pmag(4,imag)=vout(3)

        vin=[pmag(5,imag),0.0d0,pmag(6,imag)]
        call util_rotate(cen,vrot,phi,vin,vout,istatus)
        pmag(5,imag)=vout(1)
        pmag(6,imag)=vout(3)

        vin=[pmag(7,imag),0.0d0,pmag(8,imag)]
        call util_rotate(cen,vrot,phi,vin,vout,istatus)
        pmag(7,imag)=vout(1)
        pmag(8,imag)=vout(3)

      endif

      if (cposmodel.eq.'entrance') then
        dx=xpos-pmag(1,imag)
        dz=zpos-pmag(2,imag)
      else if (cposmodel.eq.'zenith') then
        dx=xpos-xopt
        dz=zpos-zopt
      else if (cposmodel.eq.'exit') then
        dx=xpos-pmag(5,imag)
        dz=zpos-pmag(6,imag)
      endif

      pmag(1,imag)=pmag(1,imag)+dx
      pmag(2,imag)=pmag(2,imag)+dz
      pmag(5,imag)=pmag(5,imag)+dx
      pmag(6,imag)=pmag(6,imag)+dz

      return
      end
