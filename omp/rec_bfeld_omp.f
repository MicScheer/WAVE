*CMZ :  4.00/07 04/06/2020  17.43.02  by  Michael Scheer
*CMZ :  3.03/02 04/03/2016  18.10.04  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.24  by  Michael Scheer
*CMZ :  2.66/07 12/05/2010  13.34.28  by  Michael Scheer
*CMZ :  2.57/05 23/10/2009  09.19.41  by  Michael Scheer
*CMZ :  2.42/04 30/01/2004  09.51.57  by  Michael Scheer
*CMZ :  2.16/04 17/07/2000  15.36.32  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.35  by  Michael Scheer
*CMZ :  2.14/02 20/04/2000  14.28.53  by  Michael Scheer
*CMZ :  1.03/06 09/06/98  14.43.04  by  Michael Scheer
*CMZ :  1.00/00 02/06/97  10.56.40  by  Michael Scheer
*CMZ : 00.01/08 22/06/95  11.22.17  by  Michael Scheer
*CMZ : 00.01/07 09/03/95  16.17.03  by  Michael Scheer
*CMZ : 00.00/00 03/03/95  15.50.11  by  Johannes Bahrdt
*-- Author : Michael Scheer
c****************************************************************
      subroutine REC_bfeld_omp(xb,yb,zb,bxout,byout,bzout)
*KEEP,gplhint.
*KEND.

      use omp_lib
      use ompmod

C *** INSIDE THIS ROUTINE UNITS OF XB,YB,ZB ... ARE [mm]

      IMPLICIT NONE

      DOUBLE PRECISION XB,YB,ZB,BX,BY,BZ,XBOLD,
     &  bxth(1024),byth(1024),bzth(1024),bxout,byout,bzout,
     &  xlenth(1024),ylenth(1024),zlenth(1024),
     &  dxth(1024),dyth(1024),dzth(1024),
     &  winrecth(1024),bcth(1024),
     &  thetath(1024),phith(1024)
      DOUBLE PRECISION XMUE0,PI,SMALL,GR
      DOUBLE PRECISION XXLEN,XXXLEN
      DOUBLE PRECISION YYLEN,YYYLEN
      DOUBLE PRECISION ZZLEN,ZZZLEN
      DOUBLE PRECISION X,Y,Z,X1,Y1,Z1,X2,Y2,Z2
      DOUBLE PRECISION X1P,Y1P,Z1P,X2P,Y2P,Z2P
      DOUBLE PRECISION X11,Y11,Z11,X12,Y12,Z12
      DOUBLE PRECISION C1,S1,S2,S3,S4,S5,S6,S7,S8,RR1,RR2
      DOUBLE PRECISION BXX,BYY,BZZ,BXXX,BYYY,BZZZ,CURRENT
      DOUBLE PRECISION COSPHI,SINPHI,SINTHE,COSTHE

c--------new 25.10.2002 Johannes Bahrdt
      DOUBLE PRECISION zz1(4),zz2(4),zz1exp(4),zz2exp(4)

      INTEGER IMAG,IMAGOLD,IMAG1,I,IZ1,IZ2,imagend,nthreads,ith

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,klotz.
      include 'klotz.cmn'
*KEEP,phycon.
      include 'phycon.cmn'
*KEND.

      DATA IMAG1/1/
      DATA IMAGOLD/1/
      DATA XBOLD/1.D30/
      DATA xmue0/1.2566d-6/
      DATA SMALL/1.0d-8/

      nthreads=min(mthreads,1024)
      nthreads=max(1,nthreads)

      pi=PI1
      gr=pi/180.d0

      bx=0.
      by=0.
      bz=0.

c************ Loop over magnets ******************************

      IF (XB.LT.XBOLD) THEN
        IMAGOLD=1
        IMAG1=1
      ENDIF

      do imag=IMAGOLD,IMAGTOT
        IF (DX(IMAG).LT.XB-WINREC*500.D0) cycle
        IF (BC(IMAG).EQ.0.0) cycle
        IF (DX(IMAG).GT.XB+WINREC*500.D0) exit
        IMAG1=IMAG
        exit
      enddo

      imagold=imag1
      imagend=imagtot
      do imag=IMAGOLD,IMAGTOT
        IF (BC(IMAG).EQ.0.0) cycle
        IF (DX(IMAG).GT.XB+WINREC*500.D0) then
          imagend=imag
          exit
        endif
      enddo

!$OMP PARALLEL NUM_THREADS(nthreads) DEFAULT(PRIVATE)
!$OMP& SHARED(bxth,byth,bzth,xlen,ylen,zlen,dx,dy,dz,winrec,bc,theta,phi,xb,yb,zb,imagold,imagend,XMUE0,PI,SMALL,GR)
      ith=OMP_GET_THREAD_NUM()+1
      bxth(ith)=0.0d0
      byth(ith)=0.0d0
      bzth(ith)=0.0d0
!$OMP DO

      do imag=IMAGOLD,imagend

        IF (BC(IMAG).EQ.0.0) cycle

        IF (DABS(THETA(IMAG)-  0.D0*GR).LT.SMALL) THEN
          COSTHE=1.0D0
          SINTHE=0.0D0
        ELSE IF (DABS(THETA(IMAG)-  90.D0*GR).LT.SMALL) THEN
          COSTHE=0.0D0
          SINTHE=1.0D0
        ELSE IF (DABS(THETA(IMAG)+  90.D0*GR).LT.SMALL) THEN
          COSTHE=0.0D0
          SINTHE=-1.0D0
        ELSE IF (DABS(THETA(IMAG)-  180.D0*GR).LT.SMALL) THEN
          COSTHE=-1.0D0
          SINTHE=0.0D0
        ELSE IF (DABS(THETA(IMAG)+  180.D0*GR).LT.SMALL) THEN
          COSTHE=-1.0D0
          SINTHE=0.0D0
        ELSE IF (DABS(THETA(IMAG)-  270.D0*GR).LT.SMALL) THEN
          COSTHE=0.0D0
          SINTHE=-1.0D0
        ELSE IF (DABS(THETA(IMAG)+  270.D0*GR).LT.SMALL) THEN
          COSTHE=0.0D0
          SINTHE=-1.0D0
        ELSE
          WRITE(6,*) THETA(IMAG)
          WRITE(6,*) '*** ERROR IN REC_BFELD: THETA WRONG (+/- 0.,90.,180.,270.) ***'
          STOP
        ENDIF

        IF (DABS(PHI(IMAG)-  0.D0*GR).LT.SMALL) THEN
          COSPHI=1.0D0
          SINPHI=0.0D0
        ELSE IF (DABS(PHI(IMAG)-  90.D0*GR).LT.SMALL) THEN
          COSPHI=0.0D0
          SINPHI=1.0D0
        ELSE IF (DABS(PHI(IMAG)+  90.D0*GR).LT.SMALL) THEN
          COSPHI=0.0D0
          SINPHI=-1.0D0
        ELSE IF (DABS(PHI(IMAG)-  180.D0*GR).LT.SMALL) THEN
          COSPHI=-1.0D0
          SINPHI=0.0D0
        ELSE IF (DABS(PHI(IMAG)+  180.D0*GR).LT.SMALL) THEN
          COSPHI=-1.0D0
          SINPHI=0.0D0
        ELSE IF (DABS(PHI(IMAG)-  270.D0*GR).LT.SMALL) THEN
          COSPHI=0.0D0
          SINPHI=-1.0D0
        ELSE IF (DABS(PHI(IMAG)+  270.D0*GR).LT.SMALL) THEN
          COSPHI=0.0D0
          SINPHI=-1.0D0
          WRITE(6,*) PHI(IMAG)
          WRITE(6,*) '*** ERROR IN REC_BFELD: PHI WRONG (+/- 0.,90.,180.,270.) ***'
          STOP
        ENDIF

c     Drehung um phi(imag)

        XXLEN=XLEN(IMAG)*cosphi+ZLEN(IMAG)*sinphi
        YYLEN=YLEN(IMAG)
        ZZLEN=ZLEN(IMAG)*cosphi-XLEN(IMAG)*sinphi

c     Drehung um theta(imag)

        XXXLEN=XXLEN*costhe-YYLEN*sinthe
        YYYLEN=YYLEN*costhe+XXLEN*sinthe
        ZZZLEN=ZZLEN

        XXXLEN=DABS(XXXLEN)
        YYYLEN=DABS(YYYLEN)
        ZZZLEN=DABS(ZZZLEN)

        X1=-XXXLEN/2.
        Y1=-YYYLEN/2.
        Z1=-ZZZLEN/2.

        X2=+XXXLEN/2.
        Y2=+YYYLEN/2.
        Z2=+ZZZLEN/2.

c     Verschiebung des Klotzes in den Koordinatenursprung
c     Neuberechnung des Aufpunktes

        x11=xb-dx(IMAG)
        y11=yb-dy(IMAG)
        z11=zb-dz(IMAG)

c     Drehung des Koordinatensystems
c     Neuberechnung des Aufpunktes
c     Drehung um phi(imag)

        x12=x11*cosphi+z11*sinphi
        y12=y11
        z12=z11*cosphi-x11*sinphi

c     Drehung um theta(imag)

        x=x12*costhe-y12*sinthe
        y=y12*costhe+x12*sinthe
        z=z12

        X1P=X-X1
        Y1P=Y-Y1
        Z1P=Z-Z1
C
        X2P=X-X2
        Y2P=Y-Y2
        Z2P=Z-Z2

        if(dabs(x1p).lt.small)x1p=small
        if(dabs(y1p).lt.small)y1p=small
        if(dabs(z1p).lt.small)z1p=small
        if(dabs(x2p).lt.small)x2p=small
        if(dabs(y2p).lt.small)y2p=small
        if(dabs(z2p).lt.small)z2p=small

C
        current=BC(IMAG)/XMUE0
        C1=(XMUE0/(4.D0*PI))*current
C
        S1=DSQRT(X2P*X2P+Y2P*Y2P+Z2P*Z2P)
        S2=DSQRT(X2P*X2P+Y1P*Y1P+Z2P*Z2P)
        S3=DSQRT(X2P*X2P+Y2P*Y2P+Z1P*Z1P)
        S4=DSQRT(X2P*X2P+Y1P*Y1P+Z1P*Z1P)
        S5=DSQRT(X1P*X1P+Y2P*Y2P+Z2P*Z2P)
        S6=DSQRT(X1P*X1P+Y1P*Y1P+Z2P*Z2P)
        S7=DSQRT(X1P*X1P+Y2P*Y2P+Z1P*Z1P)
        S8=DSQRT(X1P*X1P+Y1P*Y1P+Z1P*Z1P)
C
C  BX
C
        rr1=(Z1P+S3)*(Z2P+S2)*(Z2P+S5)*(Z1P+S8)
        rr2=(Z2P+S1)*(Z1P+S4)*(Z1P+S7)*(Z2P+S6)

        if((dabs(rr1).lt.small).and.(dabs(rr2).lt.small))then
          bxx=0.d0
        else
          BXX=C1*DLOG(rr1/rr2)
        endif

C
C    BY
C
        BYY=-C1*( DATAN(Z2P*Y2P/(X2P*S1))
     &    - DATAN(Z2P*Y1P/(X2P*S2))
     &    - DATAN(Z1P*Y2P/(X2P*S3))
     &    + DATAN(Z1P*Y1P/(X2P*S4)))
     &    +C1*( DATAN(Z2P*Y2P/(X1P*S5))
     &    - DATAN(Z2P*Y1P/(X1P*S6))
     &    - DATAN(Z1P*Y2P/(X1P*S7))
     &    + DATAN(Z1P*Y1P/(X1P*S8)))
     &    +C1*( DATAN(X2P*Y2P/(Z1P*S3))
     &    - DATAN(X1P*Y2P/(Z1P*S7))
     &    - DATAN(X2P*Y1P/(Z1P*S4))
     &    + DATAN(X1P*Y1P/(Z1P*S8)))
     &    -C1*( DATAN(X2P*Y2P/(Z2P*S1))
     &    - DATAN(X1P*Y2P/(Z2P*S5))
     &    - DATAN(X2P*Y1P/(Z2P*S2))
     &    + DATAN(X1P*Y1P/(Z2P*S6)))
C
C   BZ
C

c--------new 25.10.2002 Johannes Bahrdt
c     rr1=(X2P+S3)*(X1P+S8)*(X2P+S2)*(X1P+S5)
c     rr2=(X2P+S4)*(X1P+S7)*(X2P+S1)*(X1P+S6)
        zz1(1)=X2P+S3
        zz1(2)=X1P+S8
        zz1(3)=X2P+S2
        zz1(4)=X1P+S5

        zz2(1)=X2P+S4
        zz2(2)=X1P+S7
        zz2(3)=X2P+S1
        zz2(4)=X1P+S6

        zz1exp(1)=1.d0/x2p
        zz1exp(2)=1.d0/x1p
        zz1exp(3)=1.d0/x2p
        zz1exp(4)=1.d0/x1p

        zz2exp(1)=1.d0/x2p
        zz2exp(2)=1.d0/x1p
        zz2exp(3)=1.d0/x2p
        zz2exp(4)=1.d0/x1p

        iz1=0
        iz2=0
        do i=1,4
          if(dabs(zz1(i)).lt.small)then
            iz1=iz1+1
            zz1(i)=zz1exp(i)
          endif
          if(dabs(zz2(i)).lt.small)then
            iz2=iz2+1
            zz2(i)=zz2exp(i)
          endif
        enddo

        if(iz1.eq.iz2)then
          rr1=zz1(1)*zz1(2)*zz1(3)*zz1(4)
          rr2=zz2(1)*zz2(2)*zz2(3)*zz2(4)
        endif
        if(iz1.gt.iz2)then
          rr1=0.
          rr2=1.
        endif
        if(iz1.lt.iz2)then
          rr1=0.
          rr2=0.
          WRITE(6,*)
          WRITE(6,*) '*** WARNING IN REC_BFELD'
          WRITE(6,*) 'BZ SET TO ZERO'
          WRITE(6,*)
        endif

c-------------- end changes 25.10.2002

        rr1=(X2P+S3)*(X1P+S8)*(X2P+S2)*(X1P+S5)
        rr2=(X2P+S4)*(X1P+S7)*(X2P+S1)*(X1P+S6)

        if((dabs(rr1).lt.small).and.(dabs(rr2).lt.small))then
          bzz=0.d0
        else
          BZZ=C1*DLOG(rr1/rr2)
        endif

c     Rueckdrehung des Koordinatensystems
c     Drehung um theta(imag)
        Bxxx=Bxx*costhe+Byy*sinthe
        Byyy=Byy*costhe-Bxx*sinthe
        Bzzz=Bzz

c     Rueckdrehung um phi(imag)
        Bxth(ith)=Bxth(ith)+Bxxx*cosphi-Bzzz*sinphi
        Byth(ith)=Byth(ith)+Byyy
        Bzth(ith)=Bzth(ith)+Bzzz*cosphi+Bxxx*sinphi

c      print*,xb,imag,bxxx,byyy,bzzz

      enddo

!$OMP END DO
!$OMP END PARALLEL

      do ith=1,nthreads
        bxout=bxout+bxth(ith)
        byout=byout+byth(ith)
        bzout=bzout+bzth(ith)
      enddo

      XBOLD=XB

      return
      end
