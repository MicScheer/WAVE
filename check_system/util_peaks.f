*CMZ : 00.00/11 21/02/2011  10.51.26  by  Michael Scheer
*CMZ : 00.00/07 28/01/2008  18.09.30  by  Michael Scheer
*CMZ : 00.00/06 25/01/2008  15.07.26  by  Michael Scheer
*-- Author :    Michael Scheer   25/01/2008
      subroutine util_peaks(ndim,x,y,xpeaks,ypeaks,npeaks,pkmin,nsmooth,smooth,sigma)

c estimates peaks of spectrum. Should be positive definite for estimation
c of peak widths, since peak is assumed to be gaussian, for which approximation
c by parabola gives only right result, if sig=sqrt(-a0/(2*a2)).

      implicit none

      double precision x(ndim),y(ndim),xpeaks(ndim),ypeaks(ndim),
     &  pkmin,smooth(ndim),sigma(ndim),thresh,sm,s0,sp,x0,det,dymax

      double precision fmxtot,a1,a2,dxm,dxp,dxm2,dxp2,dxmax,sum
      integer ndim,i,npeaks,nsmooth,n1,n2,nsmooth1


      fmxtot=-1.0e30
      nsmooth=min(nsmooth,ndim)
      nsmooth=nsmooth/2
      nsmooth1=2*nsmooth+1

      sum=0.0d0
      do i=1,nsmooth1
        sum=sum+y(i)
        smooth(i)=y(i)
      enddo
      do i=ndim-nsmooth,ndim
        smooth(i)=y(i)
      enddo

      n1=0
      n2=nsmooth1
      i=nsmooth+1
      do while (n2.lt.ndim)

        i=i+1
        n1=n1+1
        n2=n2+1

        sum=sum-y(n1)+y(n2)
        smooth(i)=sum/nsmooth1

        if (smooth(i).gt.fmxtot) fmxtot=smooth(i)

      enddo

c look for all maxima above fmxtot*vpkmin
      print*,' '
      print*,'Peaks (center, height, width, distance to previous peak):'
      print*,' '

      npeaks=0
      thresh=fmxtot*pkmin

      do i=2,ndim-1

c calculate s=a0+a1*(x-x0)+a2*(x-x0)**2
c change system: (x0,s0)->(0,0), i.e.
c calculate s=a1*dx+a2*dx**2
c  ds/dx=a1+2*a2*dx_max =! 0, dx_max=-a1/2/a2

        x0=x(i)
        dxm=x(i-1)-x0
        dxp=x(i+1)-x0

        s0=smooth(i)
        sm=smooth(i-1)-s0
        sp=smooth(i+1)-s0

        if (
     &      s0.ge.thresh
     &      .and.
     &      sm.lt.0.0d0.and.sp.le.0.0d0
     &      ) then

          npeaks=npeaks+1


c           sm=a1*dxm+a2*dxm**2
c           sp=a1*dxp+a2*dxp**2

c          (dxm dxm2) (a1) = (sm)
c          (dxp dxp2) (a2) = (sp)


          dxm2=dxm*dxm
          dxp2=dxp*dxp

          det=dxm*dxp2-dxp*dxm2

          if (det.ne.0.0d0) then
            a1=(sm*dxp2-sp*dxm2)/det
            a2=(sp*dxm-sm*dxp)/det
          else
            a2=0.0d0
          endif

          if (a2.lt.0.0d0) then
            dxmax=-a1/(2.0d0*a2)
            dymax=(a1+a2*dxmax)*dxmax
            xpeaks(npeaks)=x0+dxmax
            ypeaks(npeaks)=s0+dymax
            sigma(npeaks)=sqrt(-abs(s0)/(2.0d0*a2))
          else
            xpeaks(npeaks)=x(i)
            ypeaks(npeaks)=smooth(i)
            sigma(npeaks)=-9999.0d0
          endif

        endif

      enddo

      return
      end
