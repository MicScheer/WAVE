*CMZ :  2.54/05 17/05/2005  12.13.14  by  Michael Scheer
*-- Author :
      PROGRAM BMAP_ASCII_TO_BIN

      IMPLICIT NONE

        DOUBLE PRECISION BMXMIN,BMXMAX,BMYMIN,BMYMAX,BMZMIN,BMZMAX

      REAL*4 X,Y,Z,BX,BY,BZ,XO,YO,ZO,DX,DY,DZ

        INTEGER ICODEBM,I,J,MODE
      INTEGER IBMESSX,IBMESSY,IBMESSZ
        INTEGER NBMESSX,NBMESSY,NBMESSZ

        BYTE IC

        CHARACTER(128) FILEB0,FILEASCII
        CHARACTER(65) COMMENT
        CHARACTER C

        EQUIVALENCE (IC,C)

        PRINT *,'Enter ASCII field map file:'
        ACCEPT '(A)',FILEASCII
        PRINT *,'Enter binary output filename:'
        ACCEPT '(A)',FILEB0

        OPEN(UNIT=21,FILE=FILEB0,STATUS='NEW'
     &    ,FORM='UNFORMATTED',ERR=9998)
        OPEN(UNIT=20,FILE=FILEASCII,STATUS='OLD'
     &    ,FORM='FORMATTED',ERR=9999)

        PRINT *,'Enter mode for input format'
        PRINT *,'(0 for standard format of BMESS, or non-zero for'
        PRINT *,'column format x,y,z,Bx,By,Bz:'
        ACCEPT *,MODE

        IF (MODE.EQ.0) THEN

          READ(20,'(A)')COMMENT
          DO I=1,65
            C=COMMENT(I:I)
            IF (IC.NE.9.AND.IC.NE.32) THEN !TAB AND BLANK
              DO J=I,65
                C=COMMENT(J:J)
                IF (IC.EQ.9.OR.IC.EQ.32) GOTO 9
              ENDDO
            ENDIF
          ENDDO
9         READ(COMMENT(I:(J-1)),*) ICODEBM
          COMMENT=COMMENT(J:65)
          WRITE(21)ICODEBM,COMMENT

          READ(20,*)NBMESSX,BMXMIN,BMXMAX
          READ(20,*)NBMESSY,BMYMIN,BMYMAX
          READ(20,*)NBMESSZ,BMZMIN,BMZMAX

          WRITE(21)NBMESSX,BMXMIN,BMXMAX
          WRITE(21)NBMESSY,BMYMIN,BMYMAX
          WRITE(21)NBMESSZ,BMZMIN,BMZMAX

          DO IBMESSX=1,NBMESSX
            DO IBMESSY=1,NBMESSY
              DO IBMESSZ=1,NBMESSZ

                READ(20,*)BX,BY,BZ
                WRITE(21)BX,BY,BZ

              ENDDO
            ENDDO
          ENDDO

        ELSE !IF (MODE.EQ.0) THEN

          BMXMIN=1.D30
          BMXMAX=-1.D30
          BMYMIN=1.D30
          BMYMAX=-1.D30
          BMZMIN=1.D30
          BMZMAX=-1.D30

          PRINT*,'Enter run number to be written to file:'
          READ(5,*)ICODEBM

          PRINT*,'Enter comment to be written to file:'
          READ(5,'(A)')COMMENT

          WRITE(21)ICODEBM,COMMENT

          READ(20,*,END=90)XO,YO,ZO,BX,BY,BZ
          REWIND(20)

          DX=-9999.
          DY=-9999.
          DZ=-9999.

10        CONTINUE

          READ(20,*,END=90)X,Y,Z,BX,BY,BZ

          IF (X.LT.BMXMIN) BMXMIN=X
          IF (X.GT.BMXMAX) BMXMAX=X
          IF (Y.LT.BMYMIN) BMYMIN=Y
          IF (Y.GT.BMYMAX) BMYMAX=Y
          IF (Z.LT.BMZMIN) BMZMIN=Z
          IF (Z.GT.BMZMAX) BMZMAX=Z

          IF (DX.EQ.-9999.AND.X.NE.XO) DX=X-XO
          IF (DY.EQ.-9999.AND.Y.NE.YO) DY=Y-YO
          IF (DZ.EQ.-9999.AND.Z.NE.XO) DZ=Z-ZO

          X=XO
          Y=YO
          Z=ZO

          GOTO 10
90        CONTINUE

          NBMESSX=NINT((BMXMAX-BMXMIN)/DX)+1
          NBMESSY=NINT((BMYMAX-BMYMIN)/DY)+1
          NBMESSZ=NINT((BMZMAX-BMZMIN)/DZ)+1

          PRINT*
          PRINT*,'Grid (NX,XMIN,XMAX,DX,NY,...):'
          PRINT*,NBMESSX,BMXMIN,BMXMAX,DX
          PRINT*,NBMESSY,BMYMIN,BMYMAX,DY
          PRINT*,NBMESSZ,BMZMIN,BMZMAX,DZ
          PRINT*

          WRITE(21)NBMESSX,BMXMIN,BMXMAX
          WRITE(21)NBMESSY,BMYMIN,BMYMAX
          WRITE(21)NBMESSZ,BMZMIN,BMZMAX

          REWIND(20)

          DO IBMESSX=1,NBMESSX
            DO IBMESSY=1,NBMESSY
              DO IBMESSZ=1,NBMESSZ

                READ(20,*)X,Y,Z,BX,BY,BZ
                WRITE(21)BX,BY,BZ

              ENDDO
            ENDDO
          ENDDO

        ENDIF !(MODE.EQ.0) THEN

      CLOSE(20)
        CLOSE(21)

        STOP

9998    PRINT *,'***Error while opening output file'
        STOP

9999    PRINT *,'***Error while opening input file'
        STOP

      END
