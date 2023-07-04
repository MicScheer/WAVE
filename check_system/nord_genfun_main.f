*CMZ :  3.03/04 27/11/2017  17.13.28  by  Michael Scheer
*-- Author :    Michael Scheer   27/11/2017

      program nord_genfun_main

c +PATCH,//WAVE/MAIN
c +DECK,nord_genfun_main.

      implicit none

      integer nkoef,nordng,i,j,k,l

      print*,"Enter NORDNG:"
      read(5,*)nordng

      NKOEF=0

      DO I=1,NORDNG
        DO J=1,NORDNG
          DO K=1,NORDNG
            DO L=1,NORDNG

              IF (I-1 + J-1 + K-1 + L-1 .LT. NORDNG
     &            .AND.
     &            I-1 + J-1 + K-1 + L-1 .GT. 0
     &            ) THEN
                NKOEF=NKOEF+1
              ENDIF

            ENDDO
          ENDDO
        ENDDO
      ENDDO

      print*,"Number of coefficients:",nkoef

      end
