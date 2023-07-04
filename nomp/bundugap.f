*CMZ :  4.01/02 02/05/2023  14.50.26  by  Michael Scheer
*CMZ :  4.00/11 27/03/2021  12.33.38  by  Michael Scheer
*CMZ :  4.00/07 06/05/2020  14.51.17  by  Michael Scheer
*-- Author :    Michael Scheer   05/05/2020
      subroutine bundugap(xin,yin,zin,bxout,byout,bzout,axout,ayout,azout,kini)

! Calculates Beff according to NIM A 60711 of Johannes Bahrdt and Efim Gluskin
! and uses it for the field model of subroutine bhalbasy2

      implicit none

*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,undugap.
      include 'undugap.cmn'
*KEEP,phycon.
      include 'phycon.cmn'
*KEND.

      integer kini,nundugappol
      integer :: ical=0,lung

      double precision x,y,z,bx,by,bz,ax,ay,az,r,beff,
     &  xin,yin,zin,bxout,byout,bzout,axout,ayout,azout,
     &  zkz,yky,xkx,x2,dsnzkz,dsnxkx,dshyky,dcszkz,dcsxkx,dchyky,
     &  xlundugap2,xkundugap2,bxh,byh,bzh,axh,ayh,azh,ahwmod,xl,yl,zl,
     &  ylundugap2,ykundugap2,zl2,
     &  zlundugap2,zkundugap2,poll,overhang,polwid,polhig

      double precision wlen1,totlen,totlen2,park,eharm1

      double precision atilde,btilde,ctilde,dtilde
      parameter (atilde=5.939, btilde=-11.883, ctilde=16.354,dtilde=-8.550) ! Eq. 5

      save

C--- K-VALUES

      if (kini.gt.0) ical=0

      if (ical.eq.0) then

        if (zlundugap.le.0.0d0) then
          write(lungfo,*)""
          write(lungfo,*)"*** error in bundugap: zero or negative period-lenght ***"
          write(lungfo,*)""
          write(lungfo,*)"*** Program WAVE aborted ***"
          print*,""
          print*,"*** Error in bundugap: Zero or negative period-lenght ***"
          print*,""
          print*,"*** Program WAVE aborted ***"
          stop
        endif

        r = undufullgap/zlundugap

        Beff = undugapa * exp(undugapb*r+undugapc*r*r)
        Poll = undufullgap * atilde * exp(r*(btilde+r*(ctilde+dtilde*r)))
        PolWid = 0.040 * r * 3.0
        PolHig = 0.030 * r * 3.0
        Overhang = 0.005 * r * 3.0

        xl = xlundugap
        yl = ylundugap
        zl = zlundugap

C--- K-VALUES

        XKundugap=0.0D0
        YKundugap=0.0D0
        ZKundugap=0.0D0

        IF (Zl.NE.0.0D0) ZKundugap=2.0d0*PI1/Zl
        IF (Yl.NE.0.0D0) YKundugap=2.0d0*PI1/Yl
        IF (Xl.NE.0.0D0) XKundugap=2.0d0*PI1/Xl

C--- ADJUST K-VALUES

        YKundugap=DSQRT(ZKundugap**2+XKundugap**2)
        Yl=2.0d0*PI1/YKundugap

C--- BENDING RADIUS AND DEVICE LENGTH

        PARK=ECHARGE1*DABS(beff)*Zl/(2.*PI1*EMASSKG1*CLIGHT1)
        WLEN1=(1+PARK**2/2.)/2./DMYGAMMA**2*Zl*1.0d9

        TOTLEN=Zl*((undugappol-1.0D0)/2.0D0+1.0D0)
        TOTLEN2=TOTLEN/2.0D0

        IF (WLEN1.NE.0.00) EHARM1=WTOE1/WLEN1

        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)
     &    '     Parameters of undulator model based on Beff(Gap) - Fit:'
        WRITE(LUNGFO,*)
        WRITE(LUNGFO,*)
     &    '     Fit parameters:  ',SNGL(undugapa),sngl(undugapb),sngl(undugapc)
        WRITE(LUNGFO,*)
     &    '     Full Gap and pol-length [mm]:  ',sngl(undufullgap*1000.), sngl(Poll*1000.)
        WRITE(LUNGFO,*)
     &    '     Pol-width and -height [mm]:  ',sngl(PolWid*1000.),sngl(PolHig*1000.)
        WRITE(LUNGFO,*)
     &    '     Overhang of magnets [mm]:  ',sngl(Overhang*1000.)
        WRITE(LUNGFO,*)
     &    '     Width and height of  magnets [mm]:  ',
     &    sngl(polwid+2.0d0*Overhang*1000.),
     &    sngl(polhig+Overhang*1000.)
        WRITE(LUNGFO,*)
     &    '     l0, l0x, l0y [mm]: ',
     &    SNGL(Zlundugap*1000.),SNGL(Xlundugap*1000.),SNGL(Ylundugap*1000.)
        WRITE(LUNGFO,*)
     &    '     Beff [T] and deflection parameter K: ',sngl(beff),sngl(park)
        WRITE(LUNGFO,*)
     &    '     wavelength [nm] and energy of first harmonic: ',sngl(wlen1),sngl(eharm1)

        write(lungfo,*)

        open(newunit=lung,file='bundugap.par')
        write(lung,*)"   UMAGLY_H=",sngl((polhig+overhang)*1000.0d0)
        write(lung,*)"   UMAGLZ_H=",sngl((polwid+2.0d0*overhang)*1000.0d0)
        write(lung,*)"   UPOLLY_H=",sngl((polhig*1000.0d0))
        write(lung,*)"   UPOLLZ_H=",sngl((polwid*1000.0d0))
        close(lung)

        ical=1
      endif !ical

      Nundugappol=undugappol
      AHWMOD=-ISIGN(1,-(MOD(Nundugappol,4)-2))/2.0d0

      X=XIN

      IF (DABS(XIN).GT.TOTLEN2) THEN
        BXOUT=0.0d0
        BYOUT=0.0d0
        BZOUT=0.0d0
        AXOUT=0.0d0
        AYOUT=0.0d0
        AZOUT=0.0d0
        RETURN
      ENDIF

      IF (DABS(X).LE.TOTLEN2-Zl/2.0d0) THEN

        XKX=XKundugap*(-ZIN)
        YKY=YKundugap*YIN
        ZKZ=ZKundugap*X

        DSNXKX=DSIN(XKX)
        DCSXKX=DCOS(XKX)
        DSHYKY=DSINH(YKY)
        DCHYKY=DSQRT(1.0d0+DSHYKY*DSHYKY)
        DSNZKZ=DSIN(ZKZ)
        DCSZKZ=DCOS(ZKZ)

        BXH=-XKundugap/YKundugap*beff*DSNXKX*DSHYKY*DCSZKZ
        BYH=                 beff*DCSXKX*DCHYKY*DCSZKZ
        BZH=-ZKundugap/YKundugap*beff*DCSXKX*DSHYKY*DSNZKZ

        AXH=beff/ZKundugap*                    DCSXKX*DCHYKY*DSNZKZ
        AYH=beff/ZKundugap*XKundugap/YKundugap*DSNXKX*DSHYKY*DSNZKZ
        AZH=0.0d0

        BZOUT=-BXH
        BYOUT=BYH
        BXOUT=BZH

        AZOUT=-AXH
        AYOUT=AYH
        AXOUT=AZH

        RETURN

      ELSE

        XKundugap2=XKundugap

        ZKundugap2=ZKundugap
        Zl2=2.0d0*PI1/ZKundugap2
        YKundugap2=DSQRT(ZKundugap2**2+XKundugap2**2)
        Yl=2.0d0*PI1/YKundugap2

        X2=X+TOTLEN2+Zl/2.0d0

        XKX=XKundugap2*(-ZIN)
        YKY=YKundugap2*YIN
        ZKZ=ZKundugap2*(X2)

        DSNXKX=DSIN(XKX)
        DCSXKX=DCOS(XKX)
        DSHYKY=DSINH(YKY)
        DCHYKY=DSQRT(1.0d0+DSHYKY*DSHYKY)
        DSNZKZ=DSIN(ZKZ)
        DCSZKZ=DCOS(ZKZ)

        BXH=-XKundugap2/YKundugap2*beff*DSNXKX*DSHYKY*DCSZKZ
        BYH=                 beff*DCSXKX*DCHYKY*DCSZKZ
        BZH=-ZKundugap2/YKundugap2*beff*DCSXKX*DSHYKY*DSNZKZ

        AXH=beff/ZKundugap2*                    DCSXKX*DCHYKY*DSNZKZ
        AYH=beff/ZKundugap2*XKundugap2/YKundugap2*DSNXKX*DSHYKY*DSNZKZ
        AZH=0.0d0

        ZKundugap2=ZKundugap*2.0d0
        Zl2=2.0d0*PI1/ZKundugap2
        YKundugap2=DSQRT(ZKundugap2**2+XKundugap2**2)
        Yl=2.0d0*PI1/YKundugap2

        XKX=XKundugap2*(-ZIN)
        YKY=YKundugap2*YIN
        ZKZ=ZKundugap2*(X2)

        DSNXKX=DSIN(XKX)
        DCSXKX=DCOS(XKX)
        DSHYKY=DSINH(YKY)
        DCHYKY=DSQRT(1.0d0+DSHYKY*DSHYKY)
        DSNZKZ=DSIN(ZKZ)
        DCSZKZ=DCOS(ZKZ)

        BXH=BXH-XKundugap2/YKundugap2*beff*DSNXKX*DSHYKY*DCSZKZ
        BYH=BYH+                      beff*DCSXKX*DCHYKY*DCSZKZ
        BZH=BZH-ZKundugap2/YKundugap2*beff*DCSXKX*DSHYKY*DSNZKZ

        AXH=AXH+beff/ZKundugap2*                    DCSXKX*DCHYKY*DSNZKZ
        AYH=AYH+beff/ZKundugap2*XKundugap2/YKundugap2*DSNXKX*DSHYKY*DSNZKZ
        AZH=0.0d0

        BZOUT=BXH*AHWMOD
        BYOUT=-BYH*AHWMOD
        BXOUT=-BZH*AHWMOD

        AZOUT=AXH*AHWMOD
        AYOUT=-AYH*AHWMOD
        AXOUT=-AZH*AHWMOD
      endif

      return
      end subroutine bundugap
