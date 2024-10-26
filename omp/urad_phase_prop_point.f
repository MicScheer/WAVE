*CMZ :          07/08/2024  08.52.56  by  Michael Scheer
*CMZ :  4.01/05 15/04/2024  21.58.40  by  Michael Scheer
*CMZ :  4.01/04 28/12/2023  15.30.57  by  Michael Scheer
*CMZ :  4.01/02 12/05/2023  17.13.05  by  Michael Scheer
*CMZ :  4.01/00 21/02/2023  16.51.29  by  Michael Scheer
*-- Author : Michael Scheer
      subroutine urad_phase_prop_point(sourcepoint,field,
     &  nzprop,nyprop,xprop,yprop,zprop,pinw,pinh,epho,fprop)

      implicit none

      integer nzprop,nyprop

      complex*16 :: czero=(0.0d0,0.0d0),cone=(1.0d0,0.0d0)

      double complex expom,apolh,apol45,apoll,apolr
      double complex :: field(3),cjfield(3),fprop(3,nzprop,nyprop)

      double precision sourcepoint(3),xprop,yprop(nyprop),zprop(nzprop),pinw,pinh,epho

      double precision dx,dx2,dy,dz,y,z,omc,domc,dzy2,eps(6),
     &  dr,drred,x,xobs,yobs,zobs,darlambda1,ans,stok1,stok2,stok3,stok4,zmin,ymin,fsum

      integer :: iy,iz,n,jy,jz,iobs,ieps,iobfr,jobs,jobfr,nobsv

*KEEP,phyconparam.
      include 'phyconparam.cmn'
*KEND.

      fprop=(0.0d0,0.0d0)

      if (nyprop.gt.1) then
        dy=yprop(2)-yprop(1)
      else
        dy=pinh/1000.0d0
      endif

      if (nzprop.gt.1) then
        dz=zprop(2)-zprop(1)
      else
        dz=pinw/1000.0d0
      endif

      omc=epho/(hbarev1*clight1)
      x=xprop

      xobs=sourcepoint(1)
      yobs=sourcepoint(2)
      zobs=sourcepoint(3)

      daRLAMBDA1=dz*dy*epho/WTOE1*1.0D9   !1/lambda[m]=1/(wtoe1/freq*1.e-9)
      cjfield=dconjg(field)

      do iy=1,nyprop

        y=yprop(iy)

        do iz=1,nzprop

          z=zprop(iz)

          dx=xobs-x
          dx2=dx*dx
          DY=YOBS-y
          DZ=ZOBS-z
          DZY2=DZ*DZ+DY*DY

C     TO MAKE SURE THAT TAYLOR-EXPANSION IS VALID

          IF (DZY2.GT.0.01D0*dx2) THEN
            WRITE(6,*)'*** ERROR IN URAD_phase_prop_mc ***'
            WRITE(6,*)'CHECK INPUT FILE AND INCREASE PinX'
            WRITE(6,*)'*** PROGRAM ABORTED ***'
            STOP
          ENDIF

          EPS(1)=DZY2/dx2
          DO IEPS=2,6
            EPS(IEPS)=EPS(IEPS-1)*EPS(1)
          ENDDO !IEPS

          ans=-0.0205078125D0*eps(6)+0.02734375D0*eps(5)
     &      -0.0390625D0*eps(4)+
     &      0.0625D0*eps(3)-0.125D0*eps(2)+0.5D0*eps(1)

          DR=DABS(dx*(ANS+1.0D0))
          DRRED=-DABS(dx*ANS)

          IF (DR.NE.0.0d0) THEN
            EXPOM=CDEXP(DCMPLX(0.0d0,DRRED*OMC))/DR
          ELSE
            EXPOM=1.0D0
          ENDIF

          if (dx.gt.0.0d0) then
            fprop(1:3,iz,iy)=fprop(1:3,iz,iy)+field*expom*darlambda1
          else
            fprop(1:3,iz,iy)=dconjg(fprop(1:3,iz,iy))+cjfield*expom*darlambda1
          endif

c          write(77,*)sourcepoint,z,y,dreal(fprop(1:3,iz,iy)),dimag(fprop(1:3,iz,iy))

        ENDDO !nzprop
      enddo !nyprop

      return
      end
