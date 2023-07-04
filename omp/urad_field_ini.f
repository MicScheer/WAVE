*CMZ :  4.01/02 09/05/2023  13.11.31  by  Michael Scheer
*CMZ :  4.01/00 10/02/2023  13.52.47  by  Michael Scheer
*-- Author :    Michael Scheer   10/02/2023
      subroutine urad_field_ini(perl,shift,b0v,b0h,modewave)

      implicit none

*KEEP,AMPLI.

      INTEGER IAMPDIMP
      PARAMETER (IAMPDIMP=10000)

      INTEGER IMAMPLI,IAMPSKIP,IAMPTERM,IAMPREAD,modeph
     &  ,IAMPOBSV,IAMPCOMP,IAMPREP,IAMPSUP,IAMPSEED,iamppin,iamppincirc
     &  ,iampcoh,iampincoh,nampbunchharm,noemitph,noespreadph,mampthreads

      DOUBLE PRECISION AMPFREQ,ampr2corr,AMPPHI(IAMPDIMP),AMPSHIFT(IAMPDIMP)
     &  ,AMPSCALE(IAMPDIMP),AMPRAN,FACAMPLI,phrERROR,
     &  phrshift,phrb0h,phrb0v,phrperl,phdx,phdy,phdz,phrxbeta,
     &  phrbetah,phralphah,phrbetav,phralphav,phrespread,phremith,phremitv,
     &  ampcohsig,ampbunchlen,ampbunchcharge,ampbunchp0,ampbunchr56,
     &  phrdisph,phrdispph,phrdispv,phrdisppv,phrbunlen

      COMMON/AMPLIC/
     &  ampfreq,ampr2corr,AMPPHI,AMPSHIFT,AMPSCALE,AMPRAN,FACAMPLI,phrERROR
     &  ,ampcohsig,ampbunchlen,ampbunchcharge,ampbunchp0,ampbunchr56,phrxbeta,
     &  phrshift,phrb0h,phrb0v,phrperl,phdx,phdy,phdz,
     &  phrbetah,phralphah,phrbetav,phralphav,phrespread,phremith,phremitv,
     &  phrdisph,phrdispph,phrdispv,phrdisppv,phrbunlen,
     &  IMAMPLI,IAMPSKIP,IAMPTERM,IAMPREAD,modeph
     &  ,IAMPOBSV,IAMPCOMP,IAMPREP,IAMPSUP,IAMPSEED,iamppin,iamppincirc
     &  ,iampcoh,iampincoh,nampbunchharm,noemitph,noespreadph,mampthreads

      NAMELIST/AMPLIN/
     &  IMAMPLI,IAMPSKIP,IAMPTERM,IAMPREAD
     & ,IAMPCOMP,IAMPREP
     & ,AMPSHIFT,AMPSCALE,ampfreq,ampr2corr,AMPPHI,IAMPSUP,AMPRAN,IAMPSEED
     &  ,ampcohsig,ampbunchlen,ampbunchcharge,ampbunchp0,ampbunchr56
     &  ,iampcoh,iampincoh,nampbunchharm

      NAMELIST/PHASEREPN/phrERROR,noemitph,noespreadph,modeph,phrxbeta,
     &  phrshift,phrperl,phrb0h,phrb0v,phdx,phdy,phdz,
     &  phrdisph,phrdispph,phrdispv,phrdisppv,phrbunlen,
     &  phrbetah,phralphah,phrbetav,phralphav,phrespread,phremith,phremitv
*KEND.

      double precision shift,perl,b0h,b0v
      integer modewave

      phrshift=shift
      phrperl=perl
      phrb0h=b0h
      phrb0v=b0v

      return
      end
