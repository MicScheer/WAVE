*CMZ :  3.07/00 15/03/2019  15.04.58  by  Michael Scheer
*CMZ :  3.00/00 11/03/2013  15.12.11  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.24  by  Michael Scheer
*CMZ :  2.51/02 29/04/2010  11.46.31  by  Michael Scheer
*CMZ :  2.17/00 03/11/2000  14.00.48  by  Michael Scheer
*CMZ :  2.16/08 23/10/2000  16.27.20  by  Michael Scheer
*CMZ :  2.13/03 12/01/2000  16.31.33  by  Michael Scheer
*CMZ : 00.01/02 21/11/94  11.20.25  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.56.43  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.14.09  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE wfoldsto_omp
*KEEP,gplhint.
*KEND.

*KEEP,spectf90u.
      include 'spectf90u.cmn'
*KEEP,sourcef90u.
      include 'sourcef90u.cmn'
*KEEP,observf90u.
      include 'observf90u.cmn'
*KEEP,wfoldf90u.
      include 'wfoldf90u.cmn'
*KEND.

      use omp_lib
      use wobsvmod
      use circpinmod

C--- CALCULATES FOLDING OF STOKES INTENSITIY WITH ELECTRON PHASE SPACE
C    DISTRIBUTIONS (GAUSSIAN DISTRIBUTION)

      IMPLICIT NONE

*KEEP,cmpara.
      include 'cmpara.cmn'
*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,observf90.
      include 'observf90.cmn'
*KEEP,freqs.
      include 'freqs.cmn'
*KEEP,sourcef90.
      include 'sourcef90.cmn'
*KEEP,spect.
      include 'spect.cmn'
*KEND.

      INTEGER ISTOK,kfreq,IZ,IY,i
      integer ith,kfold_th,kpincirc_th,krphi_th,kspecanaf_th,kfst

C--- CALCULATE FOURIER-COEFFICIENTS OF GAUSSIAN

C     CALL WGFOUR !ALREADY DONE, USE VALUES OF ISIGSTO

      kfold_th=ifold
      kpincirc_th=ipincirc
      krphi_th=irphi
      kspecanaf_th=ispecanaf

      DO ISTOK=1,4

!$OMP PARALLEL NUM_THREADS(mthreads) DEFAULT(PRIVATE)
!$OMP& SHARED(mthreads,obsv,obsvz,obsvy,
!$OMP& iw_circ,if1dim,ipin,ispecanaf,ifold,ipincirc,irphi,stokes,stokesf,pinr,pincen)
!$OMP& FIRSTPRIVATE(nfreq,nobsv,nsource,nobsvz,nobsvy,kfold_th,
!$OMP& kpincirc_th,krphi_th,kspecanaf_th,mobsv,mobsvy,mobsvz,istok)

      ith=OMP_GET_THREAD_NUM()+1

      allocate(
     &  x_th(max(nobsvy,nobsvz)),
     &  wobsv1_th(max(nobsvy,nobsvz)),wobsv2_th(max(nobsvy,nobsvz)),
     &  wobsv3_th(max(nobsvy,nobsvz)),wobsv4_th(max(nobsvy,nobsvz)),
     &  wobsv5_th(max(nobsvy,nobsvz)),wobsv6_th(max(nobsvy,nobsvz)),
     &  wobsv7_th(max(nobsvy,nobsvz)))

      if (ipincirc.ne.0) then
        allocate(fphir_th(nobsv))
      endif

!$OMP DO

      DO kfreq=1,NFREQ

c      do kfst=0,nfreq*4-1

c        istok=mod(kfst,4)+1
c        kfreq=kfst/4+1

c       print*,"ith, kfreq,istok:",ith,kfreq,istok

C--- CALCULATE 2D POLYNOMIALS OF INTENSITY DISTRIBUTION

        IF (ISPECANAF.EQ.0.AND.IFOLD.EQ.-2) CALL WPOLY2ST(ISTOK,kfreq)

C--- PERFORM FOLDING

        CALL WFOLISTO_omp(ISTOK,kfreq) !MUST ALSO BE CALLED FOR ISPECANAF
c        print*,"stokesf(1,18):",istok,stokesf(1,18)

C--- CALCULATE INTEGRATED INTENSITY

        IF(ISPECANAF.NE.0) CALL SPECANAF    !MUST BE CALLED EACH TIME

      enddo !kfst

!$OMP END DO

        deallocate(x_th,
     &    wobsv1_th,wobsv2_th,wobsv3_th,wobsv4_th,wobsv5_th,
     &    wobsv6_th,wobsv7_th)

        if (ipincirc.ne.0) then
          deallocate(fphir_th)
        endif

!$OMP END PARALLEL

      enddo !istok

      call blenstoffreq_omp

      RETURN
      END
