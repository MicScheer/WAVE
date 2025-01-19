*CMZ :  4.01/07 13/01/2025  22.54.56  by  Michael Scheer
*-- Author :    Michael Scheer   23/11/2024
      subroutine sbend_int_fringe(cbmodel,fringe,fa,fb,fc,ebeam,b0,bmovecut,ds,posi,icharge,
     &  istatus)

      implicit none

cdbuff dcbuff
      real*8 angle,fint,hgap,posi(7,5),ebeam,de,bmovecut,ds,v(3),dgsum,t,
     &  fringe,fa,fb,fc,beta,gamma,dtim,dtim0,angi,s,x2,x1,y2,y1,z2,z1,vx1,vx2,vy2,vy1,
     &  vz2,vz1,vn,hit(3),dist1,distn,ang,b0,ds2,xb,yb,zb,bx,by,bz,dgamma,vxp,vyp,vzp,v0

      integer :: iustep=0,ieneloss=0,icharge,istatus
      character(32) cbmodel

*KEEP,phyconparam.
      include 'phyconparam.cmn'
*KEND.

      gamma=ebeam/emassg1 !GeV
      beta=dsqrt((1.0d0-1.0d0/gamma)*(1.0d0+1.0d0/gamma))
      v0=beta*clight1

      posi=0.0d0

      posi(1:7,1)=[0.0d0,0.0d0,0.0d0,1.0d0,0.0d0,0.0d0,0.0d0]
      posi(1:7,2)=[fringe,0.0d0,0.0d0,1.0d0,0.0d0,0.0d0,0.0d0]

      if (cbmodel.eq.'hard-edge') then
        return
      endif

      t=0.0d0

      x2=0.0d0
      y2=0.0d0
      z2=0.0d0

      vx2=v0
      vy2=0.0d0
      vz2=0.0d0

      s=0.0d0

      de=0.0d0
      dtim0=ds/clight1/beta
      dtim=dtim0

      ang=atan2(vz2,vx2)
      if (ang.lt.0.0d0) ang=ang+twopi1

      !track throught fringe

      do while (.true.)

        x1=x2
        y1=y2
        z1=z2

        vx1=vx2
        vy1=vy2
        vz1=vz2

        vn=sqrt(vx1**2+(vy1**2+vz1**2))

        call util_plane_hit_hesse(
     &    posi(1:3,2),posi(4:6,2),
     &    [x1,y1,z1],[vx1,vy1,vz1],
     &    hit,dist1,distn,istatus)

        if (dist1.le.ds) exit

        xb=x1+vx1*dtim/2.0d0
        yb=y1+vy1*dtim/2.0d0
        zb=z1+vz1*dtim/2.0d0

        call bfringe(cbmodel,xb,yb,zb,b0,bx,by,bz,fringe,fa,fb,fc,icharge,istatus)

        call bmovetayl(x1,y1,z1,vx1,vy1,vz1,bx,by,bz,dtim,
     &    x2,y2,z2,vx2,vy2,vz2,vxp,vyp,vzp,gamma,icharge,bmovecut,
     &    iustep,ieneloss,dgamma)

        t=t+dtim
        s=s+sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)

        dgsum=dgsum+dgamma

      enddo !while

      dtim=dtim0*dist1/ds

      xb=x1+vx1*dtim/2.0d0
      yb=y1+vy1*dtim/2.0d0
      zb=z1+vz1*dtim/2.0d0

      call bfringe(cbmodel,xb,yb,zb,b0,bx,by,bz,fringe,fa,fb,fc,icharge,istatus)

      call bmovetayl(x1,y1,z1,vx1,vy1,vz1,bx,by,bz,dtim,
     &  x2,y2,z2,vx2,vy2,vz2,vxp,vyp,vzp,gamma,icharge,bmovecut,
     &  iustep,ieneloss,dgamma)

      t=t+dtim
      dtim=dtim0*(dist1)/ds

      s=s+sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
      vn=sqrt(vx1**2+(vy1**2+vz1**2))

      posi(1:7,2)=[x2,y2,z2,vx2/vn,vy2/vn,vz2/vn,s]

      end
