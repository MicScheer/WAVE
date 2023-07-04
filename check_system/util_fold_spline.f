*CMZ : 00.00/15 12/10/2013  12.22.24  by  Michael Scheer
*CMZ : 00.00/07 07/06/2011  14.38.25  by  Michael Scheer
*CMZ : 00.00/02 26/03/97  10.21.42  by  Michael Scheer
*CMZ : 00.00/01 26/01/95  11.40.18  by  Michael Scheer
*CMZ : 00.00/00 10/01/95  15.27.28  by  Michael Scheer
*-- Author :
      SUBROUTINE UTIL_FOLD_SPLINE(NF,XF,F,NG,XG,G,FG,ISTAT)

C--- SUBROUTINE TO EVALUATE THE FOLDED FUNCTION FG(X0)=INT{F(XF)*G(XF-X0),DXF}
c--  NG ODD!

      IMPLICIT NONE

      INTEGER NF,NG
      double precision XF(NF),F(NF),XG(NG),G(NG),FG,X0

      INTEGER NSTEP,ISTEP

      INTEGER JSTEP,NFGSTEP,IG,ISTAT

      double precision, dimension(:), allocatable::
     &  xw,fw,gw,xgf,xgff,ws1,ws2,ws3,ws4,y2f,y2g

      double precision XGCEN,XX
      double precision ANS,FOLD,XP,XM,Y1P,Y2P,Y1M,Y2M,C1P,C2P,C1M,C2M

C******************************************************
      STOP 'ERROR SR UTIL_FOLD_SPLINE: NOCH NICHT FERTIG'
C******************************************************

      ISTAT=-1

      IF ((NG/2)*2.EQ.NG) THEN
        ISTAT=1
        NG=NG-1
      ENDIF

      nstep=NF+NG

      allocate(xw(nstep))
      allocate(fw(nstep))
      allocate(gw(nstep))
      allocate(xgf(nstep))
      allocate(xgff(nstep))
      allocate(ws1(nstep))
      allocate(ws2(nstep))
      allocate(ws3(nstep))
      allocate(ws4(nstep))
      allocate(y2f(nstep))
      allocate(y2g(nstep))

C- SPLINES OF FUNCTIONS F AND G

      CALL UTIL_SPLINE_COEF(XF,F,NF,9999.0d0,9999.0d0,Y2F,WS1,WS2,WS3,WS4)
      CALL UTIL_SPLINE_COEF(XG,G,NG,9999.0d0,9999.0d0,Y2G,WS1,WS2,WS3,WS4)

C- DETERMINE RANGE OF INTEGRATION

      XGCEN=(XG(1)+XG(NG))/2.
      DO IG=1,NG
          XGF(IG)=X0+(XG(IG)-XGCEN)
      ENDDO

      NFGSTEP=NF+NG

      DO ISTEP=1,NF
          XGFF(ISTEP)=XF(ISTEP)
      ENDDO
      DO ISTEP=1,NG
          XGFF(ISTEP+NF)=XGF(ISTEP)
      ENDDO

C- SORT GRID POINTS

      CALL UTIL_SORT(NFGSTEP,XGFF)

C- CHECK ON DOUBLE COUNTING

100   IF (NFGSTEP.LT.1)
     &  STOP '*** ERROR SR UTIL_FOLD_SPLINE: STUPID INPUT ***'
      DO ISTEP=1,NFGSTEP-1
          IF (XGFF(ISTEP).EQ.XGFF(ISTEP+1)) THEN
          DO JSTEP=ISTEP+1,NFGSTEP-1
         XGFF(JSTEP)=XGFF(JSTEP+1)
          ENDDO
          NFGSTEP=NFGSTEP-1
          GOTO 100
          ENDIF
      ENDDO

C- FILL SPLINE BUFFER

      DO ISTEP=1,NFGSTEP
            XX=XGFF(ISTEP)
          IF (
     &          XX.LT.DMIN1(XF(1),XF(NF))
     &      .OR.XX.GT.DMAX1(XF(1),XF(NF))) THEN
              FW(ISTEP)=0.0D0
          ELSE
              CALL UTIL_SPLINE_INTER(XF,F,Y2F,NF,XX,FW(ISTEP),-1)
          ENDIF
          XX=(XGF(ISTEP)-X0)+XGCEN
          IF (
     &          XX.LT.DMIN1(XG(1),XG(NG))
     &      .OR.XX.GT.DMAX1(XG(1),XG(NG))) THEN
              GW(ISTEP)=0.0D0
          ELSE
              CALL UTIL_SPLINE_INTER(XG,G,Y2G,NG,XX,GW(ISTEP),-1)
          ENDIF
      ENDDO

      CALL UTIL_SPLINE_COEF(XGFF,FW,NFGSTEP,9999.0d0,9999.0d0,Y2F,WS1,WS2,WS3,WS4)
      CALL UTIL_SPLINE_COEF(XGFF,GW,NFGSTEP,9999.0d0,9999.0d0,Y2G,WS1,WS2,WS3,WS4)

C- ACTUAL INTEGRATION

      FOLD=0.0D0
      DO ISTEP=1,NFGSTEP-1
          XP=XGFF(ISTEP+1)
          XM=XGFF(ISTEP)
          Y1P=FW(ISTEP+1)
          Y2P=GW(ISTEP+1)
          Y1M=FW(ISTEP)
          Y2M=GW(ISTEP)
          C1P=Y2F(ISTEP+1)
          C2P=Y2G(ISTEP+1)
          C1M=Y2F(ISTEP)
          C2M=Y2G(ISTEP)
      ANS=(-XP**3*C1M-XP**3*C1P-XP**3*C2M-XP**3*C2P+3.*XP
     . **2*XM*C1M+3.*XP**2*XM*C1P+3.*XP**2*XM*C2M+3.*XP**2
     . *XM*C2P-3.*XP*XM**2*C1M-3.*XP*XM**2*C1P-3.*XP*XM**2
     . *C2M-3.*XP*XM**2*C2P+12.*XP*Y1M+12.*XP*Y1P+12.*XP*
     . Y2M+12.*XP*Y2P+XM**3*C1M+XM**3*C1P+XM**3*C2M+XM**3*
     . C2P-12.*XM*Y1M-12.*XM*Y1P-12.*XM*Y2M-12.*XM*Y2P)/
     . 24.
          FOLD=FOLD+ANS
      ENDDO

      FG=FOLD

      deallocate(xw)
      deallocate(fw)
      deallocate(gw)
      deallocate(xgf)
      deallocate(xgff)
      deallocate(ws1)
      deallocate(ws2)
      deallocate(ws3)
      deallocate(ws4)
      deallocate(y2f)
      deallocate(y2g)

      RETURN
      END
