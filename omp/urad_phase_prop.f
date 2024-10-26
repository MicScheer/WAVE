*CMZ :          18/10/2024  14.11.15  by  Michael Scheer
*CMZ :  4.01/05 15/04/2024  11.54.00  by  Michael Scheer
*CMZ :  4.01/04 28/12/2023  15.30.57  by  Michael Scheer
*CMZ :  4.01/02 12/05/2023  17.13.05  by  Michael Scheer
*CMZ :  4.01/00 21/02/2023  16.51.29  by  Michael Scheer
*-- Author : Michael Scheer
      subroutine urad_phase_prop(mthreads)

      use omp_lib
      use uradphasemod

      implicit none

      complex*16 :: cph00,ci=(0.0d0,1.0d0)

      integer :: mthreads,ktime=1,kfreq,icbrill,iobfr,iobsv

      if (ktime.eq.1) call util_zeit_kommentar_delta(6,'Entered urad_phase_prop',1)

      if (ifieldprop_u.gt.0) then

        if (modepin_u.ne.1) then
          call urad_phase_prop_classic(mthreads)
        else
          call urad_phase_prop_mc(mthreads)
        endif

      else !(ifieldprop_u.gt.0)

        call urad_phase_prop_geo

      endif !(ifieldprop_u.gt.0)

      icbrill=nobsvprop_u/2+1

      if (globphaseprop_u.eq.9999.0d0) then
        do kfreq=1,nepho_u
          iobfr=icbrill+nobsvprop_u*(kfreq-1)
          cph00=aradprop_u(3,iobfr)/abs(aradprop_u(3,iobfr))
          do iobsv=1,nobsvprop_u
            iobfr=iobsv+nobsvprop_u*(kfreq-1)
            aradprop_u(:,iobfr)=aradprop_u(:,iobfr)/cph00
          enddo
        enddo
      else
        aradprop_u=aradprop_u*exp(ci*globphaseprop_u)
      endif

      if (ktime.eq.1) call util_zeit_kommentar_delta(6,'Leaving urad_phase_prop',0)

      return
      end
