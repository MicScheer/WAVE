*CMZ :  3.05/24 04/12/2018  14.44.48  by  Michael Scheer
*CMZ :  3.05/09 02/08/2018  12.01.11  by  Michael Scheer
*CMZ :  3.05/08 26/07/2018  11.14.56  by  Michael Scheer
*CMZ :  3.05/07 19/07/2018  16.08.52  by  Michael Scheer
*CMZ :  3.05/06 19/07/2018  15.12.45  by  Michael Scheer
*CMZ :  3.05/05 13/07/2018  14.18.58  by  Michael Scheer
*-- Author :    Michael Scheer   09/07/2018
      program mrad_main

      use mrad_mradmod

      implicit none

      integer :: istatus=-9999, isilent=0, iwritefiles=2

      double precision
     &  elecsin(ndim1elecs_p_m,1), obsv(3,1), phener(1), ds

      double precision ::
     &  xf=-1.0d0,
     &  yf=0.0d0,
     &  zf=0.0d0,
     &  efx=1.0d0,
     &  efy=0.0d0,
     &  efz=0.0d0

      integer ::
     &  mstep=2,
     &  mthstep=1,
     &  melec=1,
     &  mobsv=1,
     &  mphener=1,
     &  modeobsrndm=1,
     &  itrackback=0,
     &  ieneloss=0,
     &  ivelofield=1,
     &  nloop=1

      call mrad_master(isilent,iwritefiles,istatus)

      if (istatus.ne.0) then
        print*,"*** Call to mrad_master returned error status:",istatus
      endif


      print*
      print*,"              --- mrad_main finished ---"
      print*

      end

      include 'mrad_master.f'
