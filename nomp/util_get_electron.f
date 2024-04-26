*CMZ :  4.01/05 15/04/2024  15.07.21  by  Michael Scheer
*CMZ :  4.01/02 11/05/2023  12.10.27  by  Michael Scheer
*-- Author :    Michael Scheer   08/05/2023
      subroutine util_get_electron(xbeta,betah,alphah,betav,alphav,emith,emitv,
     &  disph,dispph,dispv,disppv,
     &  espread,bunchlen,xelec,yelec,zelec,ypelec,zpelec,deelec,modebunch)

      implicit none

      double precision xbeta,betah,alphah,betav,alphav,emith,emitv,
     &  ebeam,xelec,yelec,zelec,ypelec,zpelec,deelec,
     &  s0h,beta0h,gamma0h,espread,bunchlen,
     &  s1h,beta1h,betap1h,alpha1h,gamma1h,phase1h,
     &  s2h,beta2h,betap2h,alpha2h,gamma2h,phase2h,
     &  s0v,beta0v,gamma0v,
     &  s1v,beta1v,betap1v,alpha1v,gamma1v,pvase1v,
     &  s2v,beta2v,betap2v,alpha2v,gamma2v,pvase2v,
     &  sigz,sigzp,sigy,sigyp,
     &  disph,dispph,dispv,disppv

      integer,parameter :: ng=6
      real rr(2),g(ng)

      integer :: modebunch

      if (emith.ne.0.0d0) then

        s1h=xbeta
        s2h=xelec
        beta1h=betah
        betap1h=-2.0d0*alphah

        call util_beta_function_drift(
     &    s0h,beta0h,gamma0h,
     &    s1h,beta1h,betap1h,alpha1h,gamma1h,phase1h,
     &    s2h,beta2h,betap2h,alpha2h,gamma2h,phase2h)

        sigz=sqrt(emith*beta0h)
        sigzp=sqrt(emith/beta0h)

      else

        sigz=0.0d0
        sigzp=0.0d0

      endif

      if (emitv.ne.0.0d0) then

        s1v=xbeta
        s2v=xelec
        beta1v=betav
        betap1v=-2.0d0*alphav

        call util_beta_function_drift(
     &    s0v,beta0v,gamma0v,
     &    s1v,beta1v,betap1v,alpha1v,gamma1v,pvase1v,
     &    s2v,beta2v,betap2v,alpha2v,gamma2v,pvase2v)

        sigy=sqrt(emitv*beta0v)
        sigyp=sqrt(emitv/beta0v)

      else

        sigy=0.0d0
        sigyp=0.0d0

      endif

      call util_random_gauss(ng,g,rr)

      deelec=g(1)*espread
      xelec=xelec+g(2)*bunchlen

      zelec=g(3)*sigz+deelec*disph
      zpelec=g(4)*sigzp+deelec*dispph

      yelec=g(5)*sigy+deelec*dispv
      ypelec=g(6)*sigyp+deelec*disppv

      zelec=zelec+(xelec-s0h)*zpelec
      yelec=yelec+(xelec-s0v)*ypelec

      return
      end
