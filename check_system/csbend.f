*CMZ :  4.00/13 28/10/2021  11.24.17  by  Michael Scheer
*CMZ :  4.00/11 02/07/2021  10.12.35  by  Michael Scheer
*CMZ :  3.06/00 11/02/2019  13.10.51  by  Michael Scheer
*CMZ :  3.05/28 07/01/2019  14.51.12  by  Michael Scheer
*-- Author :    Michael Scheer   22/12/2018
      subroutine csbend(cbmodel,strength,angle,dlength,edge,seclen,
     &  posi,fint,hgap,de,ebeam,bmovecut,ds,istatus)

      ! angle and edge in degree

      implicit none

      double precision strength,angle,dlength,seclen,posi(4,3),fringe,hgap,
     &  gamma,elmom,brho,rho,ang,fint,gap,edge(2),beta,ds,dtim,dtim0,dt,
     &  x2b,y2b,z2b,x2,y2,z2,x1,y1,z1,vxp,vyp,vzp,vx1,vy1,vz1,vx2,vy2,vz2,
     &  s1,s2,dgamma,edg1,edg2,efx2,efy2,efz2,bx2,by2,bz2,ang1,ang2,slen,dang,
     &  fringe2,fringe3,x3,fa,fb,x4,x5,y3,fringe4,fringe5,fc,ex,ey,ez,dbint,
     &  de,dde,bmovecut,ebeam

      integer istatus
      integer :: icharge=-1

      character(32) cbmodel

*KEEP,phyconparam.
      include 'phyconparam.cmn'
*KEND.
      double precision util_atan2

      istatus=0

      if (abs(angle/2.0d0-edge(1)).gt.1.0d-12 .or.
     &    abs(angle/2.0d0-edge(2)).gt.1.0d-12) then
        stop "*** Error in csbend: Edge-angle .ne. half the deflection angle ***"
      endif

      if (angle*dlength.eq.0.0d0) then
        istatus=-2
        return
      endif

      de=0.0d0

      ang=angle*grarad1
      edg1=edge(1)*grarad1
      edg2=edge(2)*grarad1
      ex=cos(edg1)
      ez=sin(edg1)
      gamma=ebeam/emassg1 !GeV
      beta=dsqrt((1.0d0-1.0d0/gamma)*(1.0d0+1.0d0/gamma))
      elmom=emassg1*dsqrt((gamma-1.0d0)*(gamma+1.0d0)) !GeV
      brho=elmom*1.0d9/clight1
      rho=dlength/ang
      strength=brho/rho
      posi=0.0d0
      gap=2.0d0*hgap
      dtim0=ds/clight1/beta
      dtim=dtim0

      if (cbmodel.eq."hard-edge") then

        posi(1,3)=-sin(ang/2.0d0)*rho
        posi(3,3)=-(1.0d0-cos(ang/2.0d0))*rho

        posi(1,1)=-posi(1,3)
        posi(3,1)=posi(3,3)

        posi(4,2)=ang*rho/2.0d0
        posi(4,3)=ang*rho

        seclen=2.0d0*abs(posi(1,3))

        dtim=rho*ang/2.0d0/clight1/beta
        dde=powcon1*strength**2*gamma*ebeam*dtim !GeV
        de=de+dde

      else if (cbmodel.eq."linear"
     &    .or.cbmodel.eq."cubic-spline"
     &    .or.cbmodel.eq."quintic-spline"
     &    ) then

        x2=0.0d0
        y2=0.0d0
        z2=0.0d0
        s2=0.0d0

        vx2=cos(edg1)*beta*clight1
        vy2=0.0d0
        vz2=sin(edg1)*beta*clight1

        if (cbmodel.eq."linear") then
          fringe=6.0d0*fint*gap
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
          stop "*** Model not yet defined in csbend ***"
        endif

        y2b=0.0d0

        do while (x2.lt.fringe+3.0d0*ds)

          x1=x2
          y1=y2
          z1=z2
          s1=s2

          vx1=vx2
          vy1=vy2
          vz1=vz2

          x2b=x1+vx1*dtim/2.0d0
          y2b=y1+vy1*dtim/2.0d0
          z2b=z1+vz1*dtim/2.0d0

          bz2=0.0d0

          if (x2b.lt.fringe) then
            if (cbmodel.eq."linear") then
              bx2=y2b/fringe
              by2=x2b/fringe
            else if (cbmodel.eq."cubic-spline") then
              x2=x2b*x2b
              x3=x2*x2b
              bx2=(2.0d0*fa*x2b+3.0d0*fb*x2)*y2b
              by2=fa*x2+fb*x3+y2b**2*(-fa-3.0d0*fb*x2b)
            else if (cbmodel.eq."quintic-spline") then

              x2=x2b*x2b
              x3=x2*x2b
              x4=x2*x2
              x5=x3*x2
              y2=y2b*y2b
              y3=y2*y2b

              bx2=y2b*(3.0d0*Fa*x2+4.0d0*fb*x3+5.0d0*fc*x4)
     &          +y3*(-fa-4.0d0*fb*x2b-10.0d0*fc*x2) !This term is not Maxwell conform
c             by2=(fa*x3+fb*x4+fc*x5)+y2*x2b*(3.0d0*fa+6.0d0*fb*x2b+10.0d0*fc*x2) ! The sign seems to be wrong in the manual
              by2=(fa*x3+fb*x4+fc*x5)-y2*x2b*(3.0d0*fa+6.0d0*fb*x2b+10.0d0*fc*x2)

            endif

            bx2=bx2*strength
            by2=by2*strength

          else
            bx2=0.0d0
            by2=strength
          endif

          dde=powcon1*by2**2*gamma*ebeam*dtim !GeV
          de=de+dde

          call bmovetayl(x1,y1,z1,vx1,vy1,vz1,bx2,by2,bz2,dtim,
     &      x2,y2,z2,vx2,vy2,vz2,vxp,vyp,vzp,gamma,icharge,bmovecut,
     &      0,0,dgamma)

          s2=s2+ds

        enddo !fringe

        ang1=util_atan2(vz1,vx1)
        ang2=util_atan2(vz2,vx2)

        x2=x2+rho*sin(ang2)
        z2=z2+rho*(1.0d0-cos(ang2))
        s2=s2+rho*ang2

        dtim=rho*ang2/clight1/beta
        dde=powcon1*strength**2*gamma*ebeam*dtim !GeV
        de=de+dde

        posi(1,1)=-x2
        posi(3,1)=-z2
        posi(1,3)=x2
        posi(3,3)=-z2
        posi(4,1)=0.0d0
        posi(4,2)=s2
        posi(4,3)=2.0d0*s2

        seclen=2.0d0*x2

      else
        stop "*** Model not yet defined in csbend ***"
      endif !model

      return
      end
