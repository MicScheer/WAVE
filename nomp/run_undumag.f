*CMZ :  4.01/02 07/05/2023  12.03.20  by  Michael Scheer
*CMZ :  4.00/17 15/11/2022  10.06.37  by  Michael Scheer
*CMZ :  4.00/16 23/07/2022  09.11.30  by  Michael Scheer
*CMZ :  4.00/15 05/07/2022  13.36.34  by  Michael Scheer
*CMZ :  4.00/14 10/02/2022  22.33.53  by  Michael Scheer
*CMZ :  4.00/11 29/04/2021  19.23.47  by  Michael Scheer
*CMZ :  4.00/10 25/09/2020  11.24.30  by  Michael Scheer
*CMZ :  4.00/07 04/06/2020  21.18.28  by  Michael Scheer
*CMZ :  4.00/04 17/05/2019  14.17.20  by  Michael Scheer
*CMZ :  4.00/03 16/04/2019  09.27.10  by  Michael Scheer
*CMZ :  4.00/01 12/04/2019  13.07.58  by  Michael Scheer
*-- Author :    Michael Scheer   11/04/2019
      subroutine run_undumag(kbundumag,lungfo)

      !use waveenv

      implicit none

*KEEP,waveenv.
      include 'waveenv.cmn'
*KEEP,klotz.
      include 'klotz.cmn'
*KEEP,undumagc.
      include 'undumagc.cmn'
*KEND.

      double precision xend,x,perlen

      integer kbundumag,luntmp,lunclc,lunmat,istat,lunnam,lund,ki,ke,kio,keo,
     &  kip,kep,kim,kem,indi,inde,nwords,ipos(2,10),
     &  iutil_fexist,lungfo,i,k,k1end,k2end,k1,k2,ixsym,nendpol,nendmag

      character(250) cval,cend
      character(16) c16
      character(32) c32
      character(4) cind,cindo,cpol,cmag
      character(2048) cline,cline1,chnamtmp,chclctmp,chundutmp

      if (FracDivFe_h.eq.0.0d0) FracDivFe_h=1.0d0
      if (fracdivfez_h.eq.0.0d0) fracdivfez_h=1.0d0

      if (kbundumag.gt.2) then
        chundutmp=trim(chwavehome) // chpathsep // "undumag" // chpathsep // trim(chundumag_h)
      else
        chundutmp=trim(chwavehome) // chpathsep // "undumag" // chpathsep // trim(chundumag)
      endif

      if (iutil_fexist(trim(chundutmp)).eq.0) goto 92

      if (kbundumag.eq.1 .or. kbundumag.eq.3 .or. kbundumag.eq.0) then

        return

      else if (kbundumag.eq.2.or.kbundumag.eq.4) then

        if (iutil_fexist(trim(chundutmp)).eq.0) goto 92

        chnamtmp=trim(chwavehome) // chpathsep // "undumag" // chpathsep // trim(chundunam)

        open(newunit=lund,file=trim(chnamtmp),err=90)
        open(newunit=lunnam,file="undumag.nam")

        if (kbundumag.eq.4) then
          nxmapu=nxmapu_h
          xmapminu=xmapminu_h
          xmapmaxu=xmapmaxu_h
          nymapu=nymapu_h
          ymapminu=ymapminu_h
          ymapmaxu=ymapmaxu_h
          nzmapu=nzmapu_h
          zmapminu=zmapminu_h
          zmapmaxu=zmapmaxu_h
        endif

        do while (.true.)

          read(lund,'(a)',end=9)cline
          call util_string_split_sep(cline,10,nwords,ipos,'=',istat)
          if (nwords.gt.0) then
            call util_string_split_pos_1(cline,k,'!',istat)
            c32=adjustl(cline(ipos(1,1):ipos(2,1)))
            call util_lower_case(c32)
          else
            c32=''
          endif

          if (c32.eq.'nuthreads') then

            write(cval,*) muthreads
            call util_string_trim(cval,k1,k2)
            call util_string_split_pos_1(cline,k,'!',istat)
            if (istat.eq.0) then
              write(lunnam,'(a)')" nuthreads=" // cval(k1:k2) // " " // cline(k:len_trim(cline))
            else
              write(lunnam,'(a)')" nuthreads=" // cval(k1:k2)
            endif

          else if (c32.eq.'ebeam') then

            write(cval,*) uebeam
            call util_string_trim(cval,k1,k2)
            call util_string_split_pos_1(cline,k,'!',istat)
            if (istat.eq.0) then
              write(lunnam,'(a)')" ebeam=" // cval(k1:k2) // " " // cline(k:len_trim(cline))
            else
              write(lunnam,'(a)')" ebeam=" // cval(k1:k2)
            endif

          else if (c32.eq.'dxmap') then

            write(cval,*) dxmapu
            call util_string_trim(cval,k1,k2)
            call util_string_split_pos_1(cline,k,'!',istat)
            if (istat.eq.0) then
              write(lunnam,'(a)')" dxmap=" // cval(k1:k2) // " " // cline(k:len_trim(cline))
            else
              write(lunnam,'(a)')" dxmap=" // cval(k1:k2)
            endif

          else if (c32.eq.'kmapmode') then
            write(lunnam,'(a)')" kmapmode=1"

          else if (c32.eq.'kmapnohead') then
            write(lunnam,'(a)')" kmapnohead=1"

          else if (c32.eq.'knointmap') then
            write(lunnam,'(a)')" knointmap=1"

          else if (c32.eq.'nxmap') then

            if (kbundumag.eq.2.and.nxmapu.lt.3) nxmapu=3

            write(cval,*) nxmapu
            call util_string_trim(cval,k1,k2)
            call util_string_split_pos_1(cline,k,'!',istat)
            if (istat.eq.0) then
              write(lunnam,'(a)')" nxmap=" // cval(k1:k2) // " " // cline(k:len_trim(cline))
            else
              write(lunnam,'(a)')" nxmap=" // cval(k1:k2)
            endif

          else if (c32.eq.'xmapmin') then

            write(cval,*) xmapminu
            call util_string_trim(cval,k1,k2)
            call util_string_split_pos_1(cline,k,'!',istat)
            if (istat.eq.0) then
              write(lunnam,'(a)')" xmapmin=" // cval(k1:k2) // " " // cline(k:len_trim(cline))
            else
              write(lunnam,'(a)')" xmapmin=" // cval(k1:k2)
            endif

          else if (c32.eq.'xmapmax') then

            write(cval,*) xmapmaxu
            call util_string_trim(cval,k1,k2)
            call util_string_split_pos_1(cline,k,'!',istat)
            if (istat.eq.0) then
              write(lunnam,'(a)')" xmapmax=" // cval(k1:k2) // " " // cline(k:len_trim(cline))
            else
              write(lunnam,'(a)')" xmapmax=" // cval(k1:k2)
            endif

          else if (c32.eq.'nymap') then

            if (nymapu.lt.3) nymapu=3

            write(cval,*) nymapu
            call util_string_trim(cval,k1,k2)
            call util_string_split_pos_1(cline,k,'!',istat)
            if (istat.eq.0) then
              write(lunnam,'(a)')" nymap=" // cval(k1:k2) // " " // cline(k:len_trim(cline))
            else
              write(lunnam,'(a)')" nymap=" // cval(k1:k2)
            endif

          else if (c32.eq.'ymapmin') then

            write(cval,*) ymapminu
            call util_string_trim(cval,k1,k2)
            call util_string_split_pos_1(cline,k,'!',istat)
            if (istat.eq.0) then
              write(lunnam,'(a)')" ymapmin=" // cval(k1:k2) // " " // cline(k:len_trim(cline))
            else
              write(lunnam,'(a)')" ymapmin=" // cval(k1:k2)
            endif

          else if (c32.eq.'ymapmax') then

            write(cval,*) ymapmaxu
            call util_string_split_pos_1(cline,k,'!',istat)
            call util_string_trim(cval,k1,k2)
            if (istat.eq.0) then
              write(lunnam,'(a)')" ymapmax=" // cval(k1:k2) // " " // cline(k:len_trim(cline))
            else
              write(lunnam,'(a)')" ymapmax=" // cval(k1:k2)
            endif

          else if (c32.eq.'nzmap') then

            if (nzmapu.lt.3) nzmapu=3

            write(cval,*) nzmapu
            call util_string_split_pos_1(cline,k,'!',istat)
            call util_string_trim(cval,k1,k2)
            if (istat.eq.0) then
              write(lunnam,'(a)')" nzmap=" // cval(k1:k2) // " " // cline(k:len_trim(cline))
            else
              write(lunnam,'(a)')" nzmap=" // cval(k1:k2)
            endif

          else if (c32.eq.'zmapmin') then

            write(cval,*) zmapminu
            call util_string_split_pos_1(cline,k,'!',istat)
            call util_string_trim(cval,k1,k2)
            if (istat.eq.0) then
              write(lunnam,'(a)')" zmapmin=" // cval(k1:k2) // " " // cline(k:len_trim(cline))
            else
              write(lunnam,'(a)')" zmapmin=" // cval(k1:k2)
            endif

          else if (c32.eq.'zmapmax') then

            write(cval,*) zmapmaxu
            call util_string_split_pos_1(cline,k,'!',istat)
            call util_string_trim(cval,k1,k2)
            if (istat.eq.0) then
              write(lunnam,'(a)')" zmapmax=" // cval(k1:k2) // " " // cline(k:len_trim(cline))
            else
              write(lunnam,'(a)')" zmapmax=" // cval(k1:k2)
            endif

          else if (c32.eq.'randox') then
            write(cval,*) urandox
            call util_string_split_pos_1(cline,k,'!',istat)
            call util_string_trim(cval,k1,k2)
            if (istat.eq.0) then
              write(lunnam,'(a)')" randox=" // cval(k1:k2) // " " // cline(k:len_trim(cline))
            else
              write(lunnam,'(a)')" randox=" // cval(k1:k2)
            endif

          else if (c32.eq.'randoy') then
            write(cval,*) urandoy
            call util_string_split_pos_1(cline,k,'!',istat)
            call util_string_trim(cval,k1,k2)
            if (istat.eq.0) then
              write(lunnam,'(a)')" randoy=" // cval(k1:k2) // " " // cline(k:len_trim(cline))
            else
              write(lunnam,'(a)')" randoy=" // cval(k1:k2)
            endif

          else if (c32.eq.'randoz') then
            write(cval,*) urandox
            call util_string_split_pos_1(cline,k,'!',istat)
            call util_string_trim(cval,k1,k2)
            if (istat.eq.0) then
              write(lunnam,'(a)')" randoz=" // cval(k1:k2) // " " // cline(k:len_trim(cline))
            else
              write(lunnam,'(a)')" randoz=" // cval(k1:k2)
            endif

          else if (c32.eq.'window') then
            write(cval,*) uwwindow
            call util_string_split_pos_1(cline,k,'!',istat)
            call util_string_trim(cval,k1,k2)
            if (istat.eq.0) then
              write(lunnam,'(a)')" window=" // cval(k1:k2) // " " // cline(k:len_trim(cline))
            else
              write(lunnam,'(a)')" window=" // cval(k1:k2)
            endif
          else if (c32.eq.'corrtiny') then
            if (kbundumag.eq.2) then
              write(cval,*) ucorrtiny
            else if (kbundumag.eq.4) then
              write(cval,*) ucorrtiny_h
            endif
            call util_string_split_pos_1(cline,k,'!',istat)
            call util_string_trim(cval,k1,k2)
            if (istat.eq.0) then
              write(lunnam,'(a)')" corrtiny=" // cval(k1:k2) // " " // cline(k:len_trim(cline))
            else
              write(lunnam,'(a)')" corrtiny=" // cval(k1:k2)
            endif

          else if (c32.eq.'ixsym') then

            call util_string_split_pos_1(cline,k,'!',istat)

            if (kbundumag.eq.4.and.umagspac_h.ne.upolspac_h.and.ixsym.ne.0) then
              print*,''
              print*,"*** Error in run_undumag: In undumag_nam.tmp is ixsym not zero, but"
              print*,"*** UMAGSPAC_H and UPOLSPAC_H are different!"
              print*,''
              stop "*** Program WAVE aborted ***"
            endif

            if (kbundumag.eq.4.and.ixsym_h.eq.1) then
              write(lunnam,'(a)')" ixsym=1 " // cline(k:len_trim(cline))
            else
              write(lunnam,'(a)')" ixsym=0 " // cline(k:len_trim(cline))
            endif

          else if (c32.eq.'iysym') then

            call util_string_split_pos_1(cline,k,'!',istat)

            if (kbundumag.eq.4.and.iysym_h.eq.1) then
              write(lunnam,'(a)')" iysym=1 " // cline(k:len_trim(cline))
            else
              write(lunnam,'(a)')" iysym=0 " // cline(k:len_trim(cline))
            endif

          else if (c32.eq.'izsym') then

            call util_string_split_pos_1(cline,k,'!',istat)

            if (kbundumag.eq.4.and.izsym_h.eq.1) then
              write(lunnam,'(a)')" izsym=1 " // cline(k:len_trim(cline))
            else
              write(lunnam,'(a)')" izsym=0 " // cline(k:len_trim(cline))
            endif

          else
            write(lunnam,'(a)')trim(cline)
          endif
        enddo

9       flush(lunnam)
        close(lunnam)
        close(lund)

        if (kbundumag.eq.2) then

          if (uairgap.lt.0.0d0) then
            print*,"*** Error in run_undumag: Negative airgap (UAIRGAP) in namelist UNDUMAGN ***"
            stop "*** Program WAVE aborted ***"
          endif

          open(newunit=lunmat,file="undumag_mu.dat")
          write(lunmat,*)umupar,uksiper," ! mu_Par and ksi_Per"
          flush(lunmat)
          close(lunmat)

          open(newunit=lunclc,file="undumag.clc")


          write(lunclc,'(a)') "* Written by WAVE"
          write(lunclc,'(a)') "* AppleII"
          write(lunclc,'(a)') " "
          write(lunclc,'(a)') "& User_Comment"
          write(cval,*)kwrun
          call util_string_trim(cval,k1,k2)
          write(lunclc,'(a)') trim(chwcom) // " (Run " // cval(k1:k2) // ")"
          write(lunclc,'(a)') " "
          write(lunclc,'(a)') ' '

          write(lunclc,'(a)') '*---------- Variables'
          write(lunclc,'(a)') ' '

          write(c16,'(I16)')nunduper
          write(lunclc,'(a)') "$nPeriods=" // adjustl(c16)
          write(c16,'(g16.6)')umaglx
          write(lunclc,'(a)') "$LxMag=" // adjustl(c16)
          write(c16,'(g16.6)')umagly
          write(lunclc,'(a)') "$LyMag=" // adjustl(c16)
          write(c16,'(g16.6)')umaglz
          write(lunclc,'(a)') "$LzMag=" // adjustl(c16)

          write(c16,'(g16.6)') uairgap
          write(lunclc,'(a)') '$AirGap=' // adjustl(c16)
          write(lunclc,'(a)') '$PerLen=4.*($LxMag+$AirGap)'

          write(c16,'(g16.6)')ucoating
          write(lunclc,'(a)') '$Mcoating=' // adjustl(c16)
          write(c16,'(g16.6)')undugap
          write(lunclc,'(a)') '$FullGap=' // adjustl(c16)
          write(lunclc,'(a)') '$matrec=1'
          write(c16,'(g16.7)')umagbc
          write(lunclc,'(a)') '$Br='  // adjustl(c16)
          write(c16,'(g16.7)')umupar
          write(lunclc,'(a)') '$Mu='  // adjustl(c16)
          write(c16,'(g16.7)')uksiper
          write(lunclc,'(a)') '$KsiPerp=' // adjustl(c16)
          write(c16,'(g16.7)')US2SHIFT
          write(lunclc,'(a)') '$S2Shift=' // adjustl(c16)
          write(c16,'(g16.7)')US3SHIFT
          write(lunclc,'(a)') '$S3Shift=' // adjustl(c16)
          write(c16,'(g16.7)')undusplit
          write(lunclc,'(a)') '$zSlit=' // adjustl(c16)
          write(lunclc,'(a)') '$xMagCen=0.0'
          write(c16,'(I16)') nudivx
          write(lunclc,'(a)') '$nMagDivX=' // adjustl(c16)
          write(c16,'(I16)') nuhdivx
          write(lunclc,'(a)') '$nHalfMagDivX=' // adjustl(c16)
          write(c16,'(I16)') nudivy
          write(lunclc,'(a)') '$nMagDivY=' // adjustl(c16)
          write(c16,'(I16)') nudivz
          write(lunclc,'(a)') '$nMagDivZ=' // adjustl(c16)

          write(lunclc,'(a)') '$E1Br = $Br / 4.'
          write(lunclc,'(a)') '$E2Br = - $Br * 3. / 4.'
          write(lunclc,'(a)') '$PerLen = 4. * ( $LxMag + $AirGap )'
          write(lunclc,'(a)') '$HalfGap = $FullGap / 2.'
          write(lunclc,'(a)') '$HalfPerLen = $PerLen / 2.'
          write(lunclc,'(a)') '$LxHalfMag = $LxMag / 2.'
          write(lunclc,'(a)') '$yMagCen = - $HalfGap - $LyMag / 2.'
          write(lunclc,'(a)') '$zMagCen = - $LzMag / 2. - $zSlit / 2.'
          write(lunclc,'(a)') '$hS3Shift = $S3Shift / 2.'
          write(lunclc,'(a)') '$hS2Shift = $S2Shift / 2.'

          write(lunclc,'(a)') '$x1LRMagCen = $xMagCen + $LxHalfMag / 2. - $hS3Shift - $hS2Shift'
          write(lunclc,'(a)') '$x2LRMagCen = $x1LRMagCen + $LxHalfMag / 2. + $AirGap + $LxMag / 2.'
          write(lunclc,'(a)') '$x3LRMagCen = $x2LRMagCen + $LxHalfMag + $AirGap + $LxHalfMag / 2.'
          write(lunclc,'(a)') '$x4LRMagCen = $x3LRMagCen + $LxHalfMag'
          write(lunclc,'(a)') '$x5LRMagCen = $x4LRMagCen + $LxHalfMag / 2. + $AirGap + $LxMag / 2.'
          write(lunclc,'(a)') '$x6LRMagCen = $x5LRMagCen + $LxMag / 2. + $AirGap + $LxHalfMag / 2.'
          write(lunclc,'(a)') '$x1LLMagCen = $x1LRMagCen + $S3Shift + $S2Shift'
          write(lunclc,'(a)') '$x2LLMagCen = $x2LRMagCen + $S3Shift + $S2Shift'
          write(lunclc,'(a)') '$x3LLMagCen = $x3LRMagCen + $S3Shift + $S2Shift'
          write(lunclc,'(a)') '$x4LLMagCen = $x4LRMagCen + $S3Shift + $S2Shift'
          write(lunclc,'(a)') '$x5LLMagCen = $x5LRMagCen + $S3Shift + $S2Shift'
          write(lunclc,'(a)') '$x6LLMagCen = $x6LRMagCen + $S3Shift + $S2Shift'
          write(lunclc,'(a)') '$x1ULMagCen = $x1LRMagCen'
          write(lunclc,'(a)') '$x2ULMagCen = $x2LRMagCen'
          write(lunclc,'(a)') '$x3ULMagCen = $x3LRMagCen'
          write(lunclc,'(a)') '$x4ULMagCen = $x4LRMagCen'
          write(lunclc,'(a)') '$x5ULMagCen = $x5LRMagCen'
          write(lunclc,'(a)') '$x6ULMagCen = $x6LRMagCen'
          write(lunclc,'(a)') '$x1URMagCen = $x1LLMagCen - $S2Shift * 2.'
          write(lunclc,'(a)') '$x2URMagCen = $x2LLMagCen - $S2Shift * 2.'
          write(lunclc,'(a)') '$x3URMagCen = $x3LLMagCen - $S2Shift * 2.'
          write(lunclc,'(a)') '$x4URMagCen = $x4LLMagCen - $S2Shift * 2.'
          write(lunclc,'(a)') '$x5URMagCen = $x5LLMagCen - $S2Shift * 2.'
          write(lunclc,'(a)') '$x6URMagCen = $x6LLMagCen - $S2Shift * 2.'
          write(lunclc,'(a)') '$yUMagCen = - $yMagCen'
          write(lunclc,'(a)') '$zLLMagCen = - $zMagCen'
          write(lunclc,'(a)') '$zULMagCen = - $zMagCen'
          write(lunclc,'(a)') '$yModCen = - 2. * $yMagCen'
          write(lunclc,'(a)') '$zModCen = - 2. * $zMagCen'

          write(lunclc,'(a)') ' '
          write(lunclc,'(a)') '$white=0'
          write(lunclc,'(a)') '$black=1'
          write(lunclc,'(a)') '$red=2'
          write(lunclc,'(a)') '$green=3'
          write(lunclc,'(a)') '$blue=4'
          write(lunclc,'(a)') '$yellow=5'
          write(lunclc,'(a)') '$magenta=6'
          write(lunclc,'(a)') '$cyan=7'
          write(lunclc,'(a)') ' '

           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '* Upstream endpoles'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '$x1E1LLMagCen = $x1LLMagCen - $PerLen'
           write(lunclc,'(a)') '$x2E1LLMagCen = $x2LLMagCen - $PerLen'
           write(lunclc,'(a)') '$x3E1LLMagCen = $x3LLMagCen - $PerLen'
           write(lunclc,'(a)') '$x1E1ULMagCen = $x1ULMagCen - $PerLen'
           write(lunclc,'(a)') '$x2E1ULMagCen = $x2ULMagCen - $PerLen'
           write(lunclc,'(a)') '$x3E1ULMagCen = $x3ULMagCen - $PerLen'
           write(lunclc,'(a)') '$x1E1LRMagCen = $x1LRMagCen - $PerLen'
           write(lunclc,'(a)') '$x2E1LRMagCen = $x2LRMagCen - $PerLen'
           write(lunclc,'(a)') '$x3E1LRMagCen = $x3LRMagCen - $PerLen'
           write(lunclc,'(a)') '$x1E1URMagCen = $x1URMagCen - $PerLen'
           write(lunclc,'(a)') '$x2E1URMagCen = $x2URMagCen - $PerLen'
           write(lunclc,'(a)') '$x3E1URMagCen = $x3URMagCen - $PerLen'
           write(lunclc,'(a)') '$x1E2LLMagCen = $x1LLMagCen - $HalfPerLen'
           write(lunclc,'(a)') '$x2E2LLMagCen = $x2LLMagCen - $HalfPerLen'
           write(lunclc,'(a)') '$x3E2LLMagCen = $x3LLMagCen - $HalfPerLen'
           write(lunclc,'(a)') '$x1E2ULMagCen = $x1ULMagCen - $HalfPerLen'
           write(lunclc,'(a)') '$x2E2ULMagCen = $x2ULMagCen - $HalfPerLen'
           write(lunclc,'(a)') '$x3E2ULMagCen = $x3ULMagCen - $HalfPerLen'
           write(lunclc,'(a)') '$x1E2LRMagCen = $x1LRMagCen - $HalfPerLen'
           write(lunclc,'(a)') '$x2E2LRMagCen = $x2LRMagCen - $HalfPerLen'
           write(lunclc,'(a)') '$x3E2LRMagCen = $x3LRMagCen - $HalfPerLen'
           write(lunclc,'(a)') '$x1E2URMagCen = $x1URMagCen - $HalfPerLen'
           write(lunclc,'(a)') '$x2E2URMagCen = $x2URMagCen - $HalfPerLen'
           write(lunclc,'(a)') '$x3E2URMagCen = $x3URMagCen - $HalfPerLen'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '* Downstream endpoles'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '$dxED0 = ( $nPeriods + 1.0 ) * $PerLen'
           write(lunclc,'(a)') '$dxED1 = ( $nPeriods + 2.0 ) * $PerLen'
           write(lunclc,'(a)') '$dxED2 = ( $nPeriods + 1.5 ) * $PerLen'
           write(lunclc,'(a)') '$x1D0LLMagCen = $x1E1LLMagCen + $dxED0'
           write(lunclc,'(a)') '$x2D0LLMagCen = $x2E1LLMagCen + $dxED0'
           write(lunclc,'(a)') '$x3D0LLMagCen = $x3E1LLMagCen + $dxED0'
           write(lunclc,'(a)') '$x1D0LRMagCen = $x1E1LRMagCen + $dxED0'
           write(lunclc,'(a)') '$x2D0LRMagCen = $x2E1LRMagCen + $dxED0'
           write(lunclc,'(a)') '$x3D0LRMagCen = $x3E1LRMagCen + $dxED0'
           write(lunclc,'(a)') '$x1D0ULMagCen = $x1E1ULMagCen + $dxED0'
           write(lunclc,'(a)') '$x2D0ULMagCen = $x2E1ULMagCen + $dxED0'
           write(lunclc,'(a)') '$x3D0ULMagCen = $x3E1ULMagCen + $dxED0'
           write(lunclc,'(a)') '$x1D0URMagCen = $x1E1URMagCen + $dxED0'
           write(lunclc,'(a)') '$x2D0URMagCen = $x2E1URMagCen + $dxED0'
           write(lunclc,'(a)') '$x3D0URMagCen = $x3E1URMagCen + $dxED0'
           write(lunclc,'(a)') '$x1D1LLMagCen = $x1E1LLMagCen + $dxED1'
           write(lunclc,'(a)') '$x2D1LLMagCen = $x2E1LLMagCen + $dxED1'
           write(lunclc,'(a)') '$x3D1LLMagCen = $x3E1LLMagCen + $dxED1'
           write(lunclc,'(a)') '$x1D1LRMagCen = $x1E1LRMagCen + $dxED1'
           write(lunclc,'(a)') '$x2D1LRMagCen = $x2E1LRMagCen + $dxED1'
           write(lunclc,'(a)') '$x3D1LRMagCen = $x3E1LRMagCen + $dxED1'
           write(lunclc,'(a)') '$x1D1ULMagCen = $x1E1ULMagCen + $dxED1'
           write(lunclc,'(a)') '$x2D1ULMagCen = $x2E1ULMagCen + $dxED1'
           write(lunclc,'(a)') '$x3D1ULMagCen = $x3E1ULMagCen + $dxED1'
           write(lunclc,'(a)') '$x1D1URMagCen = $x1E1URMagCen + $dxED1'
           write(lunclc,'(a)') '$x2D1URMagCen = $x2E1URMagCen + $dxED1'
           write(lunclc,'(a)') '$x3D1URMagCen = $x3E1URMagCen + $dxED1'
           write(lunclc,'(a)') '$x1D2LLMagCen = $x1E1LLMagCen + $dxED2'
           write(lunclc,'(a)') '$x2D2LLMagCen = $x2E1LLMagCen + $dxED2'
           write(lunclc,'(a)') '$x3D2LLMagCen = $x3E1LLMagCen + $dxED2'
           write(lunclc,'(a)') '$x1D2LRMagCen = $x1E1LRMagCen + $dxED2'
           write(lunclc,'(a)') '$x2D2LRMagCen = $x2E1LRMagCen + $dxED2'
           write(lunclc,'(a)') '$x3D2LRMagCen = $x3E1LRMagCen + $dxED2'
           write(lunclc,'(a)') '$x1D2ULMagCen = $x1E1ULMagCen + $dxED2'
           write(lunclc,'(a)') '$x2D2ULMagCen = $x2E1ULMagCen + $dxED2'
           write(lunclc,'(a)') '$x3D2ULMagCen = $x3E1ULMagCen + $dxED2'
           write(lunclc,'(a)') '$x1D2URMagCen = $x1E1URMagCen + $dxED2'
           write(lunclc,'(a)') '$x2D2URMagCen = $x2E1URMagCen + $dxED2'
           write(lunclc,'(a)') '$x3D2URMagCen = $x3E1URMagCen + $dxED2'
           write(lunclc,'(a)') '$colormag = $red'
           write(lunclc,'(a)') '$e1colormag = $magenta'
           write(lunclc,'(a)') '$e2colormag = $green'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '*---------- Magnets'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '* Lower right girder'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '& Magnet'
           write(lunclc,'(a)') 'Block mag1 HMag1 $colormag                 !key, name, mother, color'
           write(lunclc,'(a)') '$x1LRMagCen $yMagCen $zMagCen              !position'
           write(lunclc,'(a)') '$Br 1.0 0.0 0.0 $matrec                    !length bc and components of mag. vector, material index' //
     &      ''
           write(lunclc,'(a)') '$LxHalfMag $LyMag $LzMag                   !dimension'
           write(lunclc,'(a)') '$nHalfMagDivX $nMagDivY $nMagDivZ 1. 1.    !segmention'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '& Magnet'
           write(lunclc,'(a)') 'Block mag2 Mag2 $colormag                  !key, name, mother, color'
           write(lunclc,'(a)') '$x2LRMagCen $yMagCen $zMagCen              !position'
           write(lunclc,'(a)') '$Br 0.0 1.0 0.0 $matrec                    !length bc and components of mag. vector, material index' //
     &      ''
           write(lunclc,'(a)') '$LxMag $LyMag $LzMag                       !dimension'
           write(lunclc,'(a)') '$nMagDivX $nMagDivY $nMagDivZ 1. 1.        !segmention'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '& Magnet'
           write(lunclc,'(a)') 'Block mag3 HMag3 $colormag                 !key, name, mother, color'
           write(lunclc,'(a)') '$x3LRMagCen $yMagCen $zMagCen              !position'
           write(lunclc,'(a)') '$Br -1.0 0.0 0.0 $matrec                   !length bc and components of mag. vector, material index' //
     &      ''
           write(lunclc,'(a)') '$LxHalfMag $LyMag $LzMag                   !dimension'
           write(lunclc,'(a)') '$nHalfMagDivX $nMagDivY $nMagDivZ 1. 1.    !segmention'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '& Magnet'
           write(lunclc,'(a)') 'Block mag4 HMag4 $colormag                 !key, name, mother, color'
           write(lunclc,'(a)') '$x4LRMagCen $yMagCen $zMagCen              !position'
           write(lunclc,'(a)') '$Br -1.0 0.0 0.0 $matrec                   !length bc and components of mag. vector, material index' //
     &      ''
           write(lunclc,'(a)') '$LxHalfMag $LyMag $LzMag                   !dimension'
           write(lunclc,'(a)') '$nHalfMagDivX $nMagDivY $nMagDivZ 1. 1.    !segmention'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '& Magnet'
           write(lunclc,'(a)') 'Block mag5 Mag5 $colormag                  !key, name, mother, color'
           write(lunclc,'(a)') '$x5LRMagCen $yMagCen $zMagCen              !position'
           write(lunclc,'(a)') '$Br 0.0 -1.0 0.0 $matrec                   !length bc and components of mag. vector, material index' //
     &      ''
           write(lunclc,'(a)') '$LxMag $LyMag $LzMag                       !dimension'
           write(lunclc,'(a)') '$nMagDivX $nMagDivY $nMagDivZ 1. 1.        !segmention'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '& Magnet'
           write(lunclc,'(a)') 'Block mag6 HMag6 $colormag                 !key, name, mother, color'
           write(lunclc,'(a)') '$x6LRMagCen $yMagCen $zMagCen              !position'
           write(lunclc,'(a)') '$Br 1.0 0.0 0.0 $matrec                    !length bc and components of mag. vector, material index' //
     &      ''
           write(lunclc,'(a)') '$LxHalfMag $LyMag $LzMag                   !dimension'
           write(lunclc,'(a)') '$nHalfMagDivX $nMagDivY $nMagDivZ 1. 1.    !segmention'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '* Upper right girder'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '& Magnet'
           write(lunclc,'(a)') 'Block mag7 HMag7 $colormag                 !key, name, mother, color'
           write(lunclc,'(a)') '$x1URMagCen $yUMagCen $zMagCen             !position'
           write(lunclc,'(a)') '$Br -1.0 0.0 0.0 $matrec                   !length bc and components of mag. vector, material index' //
     &      ''
           write(lunclc,'(a)') '$LxHalfMag $LyMag $LzMag                   !dimension'
           write(lunclc,'(a)') '$nHalfMagDivX $nMagDivY $nMagDivZ 1. 1.    !segmention'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '& Magnet'
           write(lunclc,'(a)') 'Block mag8 Mag8 $colormag                  !key, name, mother, color'
           write(lunclc,'(a)') '$x2URMagCen $yUMagCen $zMagCen             !position'
           write(lunclc,'(a)') '$Br 0.0 1.0 0.0 $matrec                    !length bc and components of mag. vector, material index' //
     &      ''
           write(lunclc,'(a)') '$LxMag $LyMag $LzMag                       !dimension'
           write(lunclc,'(a)') '$nMagDivX $nMagDivY $nMagDivZ 1. 1.        !segmention'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '& Magnet'
           write(lunclc,'(a)') 'Block mag9 HMag9 $colormag                 !key, name, mother, color'
           write(lunclc,'(a)') '$x3URMagCen $yUMagCen $zMagCen             !position'
           write(lunclc,'(a)') '$Br 1.0 0.0 0.0 $matrec                    !length bc and components of mag. vector, material index' //
     &      ''
           write(lunclc,'(a)') '$LxHalfMag $LyMag $LzMag                   !dimension'
           write(lunclc,'(a)') '$nHalfMagDivX $nMagDivY $nMagDivZ 1. 1.    !segmention'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '& Magnet'
           write(lunclc,'(a)') 'Block mag10 HMag10 $colormag               !key, name, mother, color'
           write(lunclc,'(a)') '$x4URMagCen $yUMagCen $zMagCen             !position'
           write(lunclc,'(a)') '$Br 1.0 0.0 0.0 $matrec                    !length bc and components of mag. vector, material index' //
     &      ''
           write(lunclc,'(a)') '$LxHalfMag $LyMag $LzMag                   !dimension'
           write(lunclc,'(a)') '$nHalfMagDivX $nMagDivY $nMagDivZ 1. 1.    !segmention'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '& Magnet'
           write(lunclc,'(a)') 'Block mag11 Mag11 $colormag                !key, name, mother, color'
           write(lunclc,'(a)') '$x5URMagCen $yUMagCen $zMagCen             !position'
           write(lunclc,'(a)') '$Br 0.0 -1.0 0.0 $matrec                   !length bc and components of mag. vector, material index' //
     &      ''
           write(lunclc,'(a)') '$LxMag $LyMag $LzMag                       !dimension'
           write(lunclc,'(a)') '$nMagDivX $nMagDivY $nMagDivZ 1. 1.        !segmention'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '& Magnet'
           write(lunclc,'(a)') 'Block mag12 HMag12 $colormag               !key, name, mother, color'
           write(lunclc,'(a)') '$x6URMagCen $yUMagCen $zMagCen             !position'
           write(lunclc,'(a)') '$Br -1.0 0.0 0.0 $matrec                   !length bc and components of mag. vector, material index' //
     &      ''
           write(lunclc,'(a)') '$LxHalfMag $LyMag $LzMag                   !dimension'
           write(lunclc,'(a)') '$nHalfMagDivX $nMagDivY $nMagDivZ 1. 1.    !segmention'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '* Lower left girder'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '& Magnet'
           write(lunclc,'(a)') 'Block mag13 HMag13 $colormag               !key, name, mother, color'
           write(lunclc,'(a)') '$x1LLMagCen $yMagCen $zLLMagCen            !position'
           write(lunclc,'(a)') '$Br 1.0 0.0 0.0 $matrec                    !length bc and components of mag. vector, material index' //
     &      ''
           write(lunclc,'(a)') '$LxHalfMag $LyMag $LzMag                   !dimension'
           write(lunclc,'(a)') '$nHalfMagDivX $nMagDivY $nMagDivZ 1. 1.    !segmention'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '& Magnet'
           write(lunclc,'(a)') 'Block mag14 Mag14 $colormag                !key, name, mother, color'
           write(lunclc,'(a)') '$x2LLMagCen $yMagCen $zLLMagCen            !position'
           write(lunclc,'(a)') '$Br 0.0 1.0 0.0 $matrec                    !length bc and components of mag. vector, material index' //
     &      ''
           write(lunclc,'(a)') '$LxMag $LyMag $LzMag                       !dimension'
           write(lunclc,'(a)') '$nMagDivX $nMagDivY $nMagDivZ 1. 1.        !segmention'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '& Magnet'
           write(lunclc,'(a)') 'Block mag15 HMag15 $colormag               !key, name, mother, color'
           write(lunclc,'(a)') '$x3LLMagCen $yMagCen $zLLMagCen            !position'
           write(lunclc,'(a)') '$Br -1.0 0.0 0.0 $matrec                   !length bc and components of mag. vector, material index' //
     &      ''
           write(lunclc,'(a)') '$LxHalfMag $LyMag $LzMag                   !dimension'
           write(lunclc,'(a)') '$nHalfMagDivX $nMagDivY $nMagDivZ 1. 1.    !segmention'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '& Magnet'
           write(lunclc,'(a)') 'Block mag16 HMag16 $colormag               !key, name, mother, color'
           write(lunclc,'(a)') '$x4LLMagCen $yMagCen $zLLMagCen            !position'
           write(lunclc,'(a)') '$Br -1.0 0.0 0.0 $matrec                   !length bc and components of mag. vector, material index' //
     &      ''
           write(lunclc,'(a)') '$LxHalfMag $LyMag $LzMag                   !dimension'
           write(lunclc,'(a)') '$nHalfMagDivX $nMagDivY $nMagDivZ 1. 1.    !segmention'
           write(lunclc,'(a)') '& Magnet'
           write(lunclc,'(a)') 'Block mag17 Mag17 $colormag                !key, name, mother, color'
           write(lunclc,'(a)') '$x5LLMagCen $yMagCen $zLLMagCen            !position'
           write(lunclc,'(a)') '$Br 0.0 -1.0 0.0 $matrec                   !length bc and components of mag. vector, material index' //
     &      ''
           write(lunclc,'(a)') '$LxMag $LyMag $LzMag                       !dimension'
           write(lunclc,'(a)') '$nMagDivX $nMagDivY $nMagDivZ 1. 1.        !segmention'
           write(lunclc,'(a)') '& Magnet'
           write(lunclc,'(a)') 'Block mag18 HMag18 $colormag               !key, name, mother, color'
           write(lunclc,'(a)') '$x6LLMagCen $yMagCen $zLLMagCen            !position'
           write(lunclc,'(a)') '$Br 1.0 0.0 0.0 $matrec                    !length bc and components of mag. vector, material index' //
     &      ''
           write(lunclc,'(a)') '$LxHalfMag $LyMag $LzMag                   !dimension'
           write(lunclc,'(a)') '$nHalfMagDivX $nMagDivY $nMagDivZ 1. 1.    !segmention'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '* Upper right girder'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '& Magnet'
           write(lunclc,'(a)') 'Block mag19 HMag19 $colormag               !key, name, mother, color'
           write(lunclc,'(a)') '$x1LRMagCen $yUMagCen $zULMagCen           !position'
           write(lunclc,'(a)') '$Br -1.0 0.0 0.0 $matrec                   !length bc and components of mag. vector, material index' //
     &      ''
           write(lunclc,'(a)') '$LxHalfMag $LyMag $LzMag                   !dimension'
           write(lunclc,'(a)') '$nHalfMagDivX $nMagDivY $nMagDivZ 1. 1.    !segmention'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '& Magnet'
           write(lunclc,'(a)') 'Block mag20 Mag20 $colormag                !key, name, mother, color'
           write(lunclc,'(a)') '$x2LRMagCen $yUMagCen $zULMagCen           !position'
           write(lunclc,'(a)') '$Br 0.0 1.0 0.0 $matrec                    !length bc and components of mag. vector, material index' //
     &      ''
           write(lunclc,'(a)') '$LxMag $LyMag $LzMag                       !dimension'
           write(lunclc,'(a)') '$nMagDivX $nMagDivY $nMagDivZ 1. 1.        !segmention'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '& Magnet'
           write(lunclc,'(a)') 'Block mag21 HMag21 $colormag               !key, name, mother, color'
           write(lunclc,'(a)') '$x3LRMagCen $yUMagCen $zULMagCen           !position'
           write(lunclc,'(a)') '$Br 1.0 0.0 0.0 $matrec                    !length bc and components of mag. vector, material index' //
     &      ''
           write(lunclc,'(a)') '$LxHalfMag $LyMag $LzMag                   !dimension'
           write(lunclc,'(a)') '$nHalfMagDivX $nMagDivY $nMagDivZ 1. 1.    !segmention'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '& Magnet'
           write(lunclc,'(a)') 'Block mag22 HMag22 $colormag               !key, name, mother, color'
           write(lunclc,'(a)') '$x4LRMagCen $yUMagCen $zULMagCen           !position'
           write(lunclc,'(a)') '$Br 1.0 0.0 0.0 $matrec                    !length bc and components of mag. vector, material index' //
     &      ''
           write(lunclc,'(a)') '$LxHalfMag $LyMag $LzMag                   !dimension'
           write(lunclc,'(a)') '$nHalfMagDivX $nMagDivY $nMagDivZ 1. 1.    !segmention'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '& Magnet'
           write(lunclc,'(a)') 'Block mag23 Mag23 $colormag                !key, name, mother, color'
           write(lunclc,'(a)') '$x5LRMagCen $yUMagCen $zULMagCen           !position'
           write(lunclc,'(a)') '$Br 0.0 -1.0 0.0 $matrec                   !length bc and components of mag. vector, material index' //
     &      ''
           write(lunclc,'(a)') '$LxMag $LyMag $LzMag                       !dimension'
           write(lunclc,'(a)') '$nMagDivX $nMagDivY $nMagDivZ 1. 1.        !segmention'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '& Magnet'
           write(lunclc,'(a)') 'Block mag24 HMag24 $colormag               !key, name, mother, color'
           write(lunclc,'(a)') '$x6LRMagCen $yUMagCen $zULMagCen           !position'
           write(lunclc,'(a)') '$Br -1.0 0.0 0.0 $matrec                   !length bc and components of mag. vector, material index' //
     &      ''
           write(lunclc,'(a)') '$LxHalfMag $LyMag $LzMag                   !dimension'
           write(lunclc,'(a)') '$nHalfMagDivX $nMagDivY $nMagDivZ 1. 1.    !segmention'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '* Lower right girder'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '& Special_Magnet'
           write(lunclc,'(a)') 'Block mag25 Mag25 $e1colormag              !key, name, mother, color'
           write(lunclc,'(a)') '$x1E1LRMagCen $yMagCen $zMagCen            !position'
           write(lunclc,'(a)') '$E1Br 1.0 0.0 0.0 $matrec                  !length bc and components of mag. vector, material index' //
     &      ''
           write(lunclc,'(a)') '$LxHalfMag $LyMag $LzMag                   !dimension'
           write(lunclc,'(a)') '$nHalfMagDivX $nMagDivY $nMagDivZ 1. 1.    !segmention'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '& Special_Magnet'
           write(lunclc,'(a)') 'Block mag26 Mag26 $e1colormag              !key, name, mother, color'
           write(lunclc,'(a)') '$x2E1LRMagCen $yMagCen $zMagCen            !position'
           write(lunclc,'(a)') '$E1Br 0.0 1.0 0.0 $matrec                  !length bc and components of mag. vector, material index' //
     &      ''
           write(lunclc,'(a)') '$LxMag $LyMag $LzMag                       !dimension'
           write(lunclc,'(a)') '$nMagDivX $nMagDivY $nMagDivZ 1. 1.        !segmention'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '& Special_Magnet'
           write(lunclc,'(a)') 'Block mag27 HMag27 $e1colormag             !key, name, mother, color'
           write(lunclc,'(a)') '$x3E1LRMagCen $yMagCen $zMagCen            !position'
           write(lunclc,'(a)') '$E1Br -1.0 0.0 0.0 $matrec                 !length bc and components of mag. vector, material index' //
     &      ''
           write(lunclc,'(a)') '$LxHalfMag $LyMag $LzMag                   !dimension'
           write(lunclc,'(a)') '$nHalfMagDivX $nMagDivY $nMagDivZ 1. 1.    !segmention'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '* Upper right girder'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '& Special_Magnet'
           write(lunclc,'(a)') 'Block mag28 HMag28 $e1colormag             !key, name, mother, color'
           write(lunclc,'(a)') '$x1E1URMagCen $yUMagCen $zMagCen           !position'
           write(lunclc,'(a)') '$E1Br -1.0 0.0 0.0 $matrec                 !length bc and components of mag. vector, material index' //
     &      ''
           write(lunclc,'(a)') '$LxHalfMag $LyMag $LzMag                   !dimension'
           write(lunclc,'(a)') '$nHalfMagDivX $nMagDivY $nMagDivZ 1. 1.    !segmention'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '& Special_Magnet'
           write(lunclc,'(a)') 'Block mag29 Mag29 $e1colormag              !key, name, mother, color'
           write(lunclc,'(a)') '$x2E1URMagCen $yUMagCen $zMagCen           !position'
           write(lunclc,'(a)') '$E1Br 0.0 1.0 0.0 $matrec                  !length bc and components of mag. vector, material index' //
     &      ''
           write(lunclc,'(a)') '$LxMag $LyMag $LzMag                       !dimension'
           write(lunclc,'(a)') '$nMagDivX $nMagDivY $nMagDivZ 1. 1.        !segmention'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '& Special_Magnet'
           write(lunclc,'(a)') 'Block mag30 HMag30 $e1colormag             !key, name, mother, color'
           write(lunclc,'(a)') '$x3E1URMagCen $yUMagCen $zMagCen           !position'
           write(lunclc,'(a)') '$E1Br 1.0 0.0 0.0 $matrec                  !length bc and components of mag. vector, material index' //
     &      ''
           write(lunclc,'(a)') '$LxHalfMag $LyMag $LzMag                   !dimension'
           write(lunclc,'(a)') '$nHalfMagDivX $nMagDivY $nMagDivZ 1. 1.    !segmention'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '* Lower left girder'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '& Special_Magnet'
           write(lunclc,'(a)') 'Block mag31 HMag31 $e1colormag              !key, name, mother, color'
           write(lunclc,'(a)') '$x1E1LLMagCen $yMagCen $zLLMagCen           !position'
           write(lunclc,'(a)') '$E1Br 1.0 0.0 0.0 $matrec                   !length bc and components of mag. vector, material index' //
     &      ''
           write(lunclc,'(a)') '$LxHalfMag $LyMag $LzMag                    !dimension'
           write(lunclc,'(a)') '$nHalfMagDivX $nMagDivY $nMagDivZ 1. 1.     !segmention'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '& Special_Magnet'
           write(lunclc,'(a)') 'Block mag32 Mag32 $e1colormag               !key, name, mother, color'
           write(lunclc,'(a)') '$x2E1LLMagCen $yMagCen $zLLMagCen           !position'
           write(lunclc,'(a)') '$E1Br 0.0 1.0 0.0 $matrec                   !length bc and components of mag. vector, material index' //
     &      ''
           write(lunclc,'(a)') '$LxMag $LyMag $LzMag                        !dimension'
           write(lunclc,'(a)') '$nMagDivX $nMagDivY $nMagDivZ 1. 1.         !segmention'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '& Special_Magnet'
           write(lunclc,'(a)') 'Block mag33 HMag33 $e1colormag              !key, name, mother, color'
           write(lunclc,'(a)') '$x3E1LLMagCen $yMagCen $zLLMagCen           !position'
           write(lunclc,'(a)') '$E1Br -1.0 0.0 0.0 $matrec                  !length bc and components of mag. vector, material index' //
     &      ''
           write(lunclc,'(a)') '$LxHalfMag $LyMag $LzMag                    !dimension'
           write(lunclc,'(a)') '$nHalfMagDivX $nMagDivY $nMagDivZ 1. 1.     !segmention'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '* Lower right girder'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '& Special_Magnet'
           write(lunclc,'(a)') 'Block mag34 HMag34 $e1colormag              !key, name, mother, color'
           write(lunclc,'(a)') '$x1E1LRMagCen $yUMagCen $zULMagCen          !position'
           write(lunclc,'(a)') '$E1Br -1.0 0.0 0.0 $matrec                  !length bc and components of mag. vector, material index' //
     &      ''
           write(lunclc,'(a)') '$LxHalfMag $LyMag $LzMag                    !dimension'
           write(lunclc,'(a)') '$nHalfMagDivX $nMagDivY $nMagDivZ 1. 1.     !segmention'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '& Special_Magnet'
           write(lunclc,'(a)') 'Block mag35 Mag35 $e1colormag               !key, name, mother, color'
           write(lunclc,'(a)') '$x2E1LRMagCen $yUMagCen $zULMagCen          !position'
           write(lunclc,'(a)') '$E1Br 0.0 1.0 0.0 $matrec                   !length bc and components of mag. vector, material index' //
     &      ''
           write(lunclc,'(a)') '$LxMag $LyMag $LzMag                        !dimension'
           write(lunclc,'(a)') '$nMagDivX $nMagDivY $nMagDivZ 1. 1.         !segmention'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '& Special_Magnet'
           write(lunclc,'(a)') 'Block mag36 HMag36 $e1colormag              !key, name, mother, color'
           write(lunclc,'(a)') '$x3E1LRMagCen $yUMagCen $zULMagCen          !position'
           write(lunclc,'(a)') '$E1Br 1.0 0.0 0.0 $matrec                   !length bc and components of mag. vector, material index' //
     &      ''
           write(lunclc,'(a)') '$LxHalfMag $LyMag $LzMag                    !dimension'
           write(lunclc,'(a)') '$nHalfMagDivX $nMagDivY $nMagDivZ 1. 1.     !segmention'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '& Special_Magnet'
           write(lunclc,'(a)') 'Block mag37 HMag37 $e2colormag              !key, name, mother, color'
           write(lunclc,'(a)') '$x1E2LRMagCen $yMagCen $zMagCen             !position'
           write(lunclc,'(a)') '$E2Br 1.0 0.0 0.0 $matrec                   !length bc and components of mag. vector, material index' //
     &      ''
           write(lunclc,'(a)') '$LxHalfMag $LyMag $LzMag                    !dimension'
           write(lunclc,'(a)') '$nHalfMagDivX $nMagDivY $nMagDivZ 1. 1.     !segmention'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '& Special_Magnet'
           write(lunclc,'(a)') 'Block mag38 Mag38 $e2colormag               !key, name, mother, color'
           write(lunclc,'(a)') '$x2E2LRMagCen $yMagCen $zMagCen             !position'
           write(lunclc,'(a)') '$E2Br 0.0 1.0 0.0 $matrec                   !length bc and components of mag. vector, material index' //
     &      ''
           write(lunclc,'(a)') '$LxMag $LyMag $LzMag                        !dimension'
           write(lunclc,'(a)') '$nMagDivX $nMagDivY $nMagDivZ 1. 1.         !segmention'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '& Special_Magnet'
           write(lunclc,'(a)') 'Block mag39 HMag39 $e2colormag              !key, name, mother, color'
           write(lunclc,'(a)') '$x3E2LRMagCen $yMagCen $zMagCen             !position'
           write(lunclc,'(a)') '$E2Br -1.0 0.0 0.0 $matrec                  !length bc and components of mag. vector, material index' //
     &      ''
           write(lunclc,'(a)') '$LxHalfMag $LyMag $LzMag                    !dimension'
           write(lunclc,'(a)') '$nHalfMagDivX $nMagDivY $nMagDivZ 1. 1.     !segmention'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '* Upper right girder'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '& Special_Magnet'
           write(lunclc,'(a)') 'Block mag40 HMag40 $e2colormag              !key, name, mother, color'
           write(lunclc,'(a)') '$x1E2URMagCen $yUMagCen $zMagCen            !position'
           write(lunclc,'(a)') '$E2Br -1.0 0.0 0.0 $matrec                  !length bc and components of mag. vector, material index' //
     &      ''
           write(lunclc,'(a)') '$LxHalfMag $LyMag $LzMag                    !dimension'
           write(lunclc,'(a)') '$nHalfMagDivX $nMagDivY $nMagDivZ 1. 1.     !segmention'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '& Special_Magnet'
           write(lunclc,'(a)') 'Block mag41 Mag41 $e2colormag               !key, name, mother, color'
           write(lunclc,'(a)') '$x2E2URMagCen $yUMagCen $zMagCen            !position'
           write(lunclc,'(a)') '$E2Br 0.0 1.0 0.0 $matrec                   !length bc and components of mag. vector, material index' //
     &      ''
           write(lunclc,'(a)') '$LxMag $LyMag $LzMag                        !dimension'
           write(lunclc,'(a)') '$nMagDivX $nMagDivY $nMagDivZ 1. 1.         !segmention'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '& Special_Magnet'
           write(lunclc,'(a)') 'Block mag42 HMag42 $e2colormag              !key, name, mother, color'
           write(lunclc,'(a)') '$x3E2URMagCen $yUMagCen $zMagCen            !position'
           write(lunclc,'(a)') '$E2Br 1.0 0.0 0.0 $matrec                   !length bc and components of mag. vector, material index' //
     &      ''
           write(lunclc,'(a)') '$LxHalfMag $LyMag $LzMag                    !dimension'
           write(lunclc,'(a)') '$nHalfMagDivX $nMagDivY $nMagDivZ 1. 1.     !segmention'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '* Lower left girder'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '& Special_Magnet'
           write(lunclc,'(a)') 'Block mag43 HMag43 $e2colormag              !key, name, mother, color'
           write(lunclc,'(a)') '$x1E2LLMagCen $yMagCen $zLLMagCen           !position'
           write(lunclc,'(a)') '$E2Br 1.0 0.0 0.0 $matrec                   !length bc and components of mag. vector, material index' //
     &      ''
           write(lunclc,'(a)') '$LxHalfMag $LyMag $LzMag                    !dimension'
           write(lunclc,'(a)') '$nHalfMagDivX $nMagDivY $nMagDivZ 1. 1.     !segmention'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '& Special_Magnet'
           write(lunclc,'(a)') 'Block mag44 Mag44 $e2colormag               !key, name, mother, color'
           write(lunclc,'(a)') '$x2E2LLMagCen $yMagCen $zLLMagCen           !position'
           write(lunclc,'(a)') '$E2Br 0.0 1.0 0.0 $matrec                   !length bc and components of mag. vector, material index' //
     &      ''
           write(lunclc,'(a)') '$LxMag $LyMag $LzMag                        !dimension'
           write(lunclc,'(a)') '$nMagDivX $nMagDivY $nMagDivZ 1. 1.         !segmention'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '& Special_Magnet'
           write(lunclc,'(a)') 'Block mag45 HMag45 $e2colormag              !key, name, mother, color'
           write(lunclc,'(a)') '$x3E2LLMagCen $yMagCen $zLLMagCen           !position'
           write(lunclc,'(a)') '$E2Br -1.0 0.0 0.0 $matrec                  !length bc and components of mag. vector, material index' //
     &      ''
           write(lunclc,'(a)') '$LxHalfMag $LyMag $LzMag                    !dimension'
           write(lunclc,'(a)') '$nHalfMagDivX $nMagDivY $nMagDivZ 1. 1.     !segmention'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '* Lower right girder'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '& Special_Magnet'
           write(lunclc,'(a)') 'Block mag46 HMag46 $e2colormag               !key, name, mother, color'
           write(lunclc,'(a)') '$x1E2LRMagCen $yUMagCen $zULMagCen           !position'
           write(lunclc,'(a)') '$E2Br -1.0 0.0 0.0 $matrec                   !length bc and components of mag. vector, material index' //
     &      ''
           write(lunclc,'(a)') '$LxHalfMag $LyMag $LzMag                     !dimension'
           write(lunclc,'(a)') '$nHalfMagDivX $nMagDivY $nMagDivZ 1. 1.      !segmention'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '& Special_Magnet'
           write(lunclc,'(a)') 'Block mag47 Mag47 $e2colormag                !key, name, mother, color'
           write(lunclc,'(a)') '$x2E2LRMagCen $yUMagCen $zULMagCen           !position'
           write(lunclc,'(a)') '$E2Br 0.0 1.0 0.0 $matrec                    !length bc and components of mag. vector, material index' //
     &      ''
           write(lunclc,'(a)') '$LxMag $LyMag $LzMag                         !dimension'
           write(lunclc,'(a)') '$nMagDivX $nMagDivY $nMagDivZ 1. 1.          !segmention'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '& Special_Magnet'
           write(lunclc,'(a)') 'Block mag48 HMag48 $e2colormag               !key, name, mother, color'
           write(lunclc,'(a)') '$x3E2LRMagCen $yUMagCen $zULMagCen           !position'
           write(lunclc,'(a)') '$E2Br 1.0 0.0 0.0 $matrec                    !length bc and components of mag. vector, material index' //
     &      ''
           write(lunclc,'(a)') '$LxHalfMag $LyMag $LzMag                     !dimension'
           write(lunclc,'(a)') '$nHalfMagDivX $nMagDivY $nMagDivZ 1. 1.      !segmention'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '& Special_Magnet'
           write(lunclc,'(a)') 'Block mag49 HMag49 $e2colormag               !key, name, mother, color'
           write(lunclc,'(a)') '$x1D0LRMagCen $yMagCen $zMagCen              !position'
           write(lunclc,'(a)') '$Br 1.0 0.0 0.0 $matrec                      !length bc and components of mag. vector, material index' //
     &      ''
           write(lunclc,'(a)') '$LxHalfMag $LyMag $LzMag                     !dimension'
           write(lunclc,'(a)') '$nHalfMagDivX $nMagDivY $nMagDivZ 1. 1.      !segmention'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '& Special_Magnet'
           write(lunclc,'(a)') 'Block mag50 Mag50 $e2colormag                !key, name, mother, color'
           write(lunclc,'(a)') '$x2D0LRMagCen $yMagCen $zMagCen              !position'
           write(lunclc,'(a)') '$Br 0.0 1.0 0.0 $matrec                      !length bc and components of mag. vector, material index' //
     &      ''
           write(lunclc,'(a)') '$LxMag $LyMag $LzMag                         !dimension'
           write(lunclc,'(a)') '$nMagDivX $nMagDivY $nMagDivZ 1. 1.          !segmention'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '& Special_Magnet'
           write(lunclc,'(a)') 'Block mag51 HMag51 $e2colormag               !key, name, mother, color'
           write(lunclc,'(a)') '$x3D0LRMagCen $yMagCen $zMagCen              !position'
           write(lunclc,'(a)') '$Br -1.0 0.0 0.0 $matrec                     !length bc and components of mag. vector, material index' //
     &      ''
           write(lunclc,'(a)') '$LxHalfMag $LyMag $LzMag                     !dimension'
           write(lunclc,'(a)') '$nHalfMagDivX $nMagDivY $nMagDivZ 1. 1.      !segmention'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '& Special_Magnet'
           write(lunclc,'(a)') 'Block mag52 HMag52 $e1colormag               !key, name, mother, color'
           write(lunclc,'(a)') '$x1D1LRMagCen $yMagCen $zMagCen              !position'
           write(lunclc,'(a)') '$E1Br 1.0 0.0 0.0 $matrec                    !length bc and components of mag. vector, material index' //
     &      ''
           write(lunclc,'(a)') '$LxHalfMag $LyMag $LzMag                     !dimension'
           write(lunclc,'(a)') '$nHalfMagDivX $nMagDivY $nMagDivZ 1. 1.      !segmention'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '& Special_Magnet'
           write(lunclc,'(a)') 'Block mag53 Mag53 $e1colormag                !key, name, mother, color'
           write(lunclc,'(a)') '$x2D1LRMagCen $yMagCen $zMagCen              !position'
           write(lunclc,'(a)') '$E1Br 0.0 1.0 0.0 $matrec                    !length bc and components of mag. vector, material index' //
     &      ''
           write(lunclc,'(a)') '$LxMag $LyMag $LzMag                         !dimension'
           write(lunclc,'(a)') '$nMagDivX $nMagDivY $nMagDivZ 1. 1.          !segmention'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '& Special_Magnet'
           write(lunclc,'(a)') 'Block mag54 HMag54 $e1colormag               !key, name, mother, color'
           write(lunclc,'(a)') '$x3D1LRMagCen $yMagCen $zMagCen              !position'
           write(lunclc,'(a)') '$E1Br -1.0 0.0 0.0 $matrec                   !length bc and components of mag. vector, material index' //
     &      ''
           write(lunclc,'(a)') '$LxHalfMag $LyMag $LzMag                     !dimension'
           write(lunclc,'(a)') '$nHalfMagDivX $nMagDivY $nMagDivZ 1. 1.      !segmention'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '* Upper right girder'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '& Special_Magnet'
           write(lunclc,'(a)') 'Block mag55 HMag55 $e2colormag               !key, name, mother, color'
           write(lunclc,'(a)') '$x1D0URMagCen $yUMagCen $zMagCen             !position'
           write(lunclc,'(a)') '$Br -1.0 0.0 0.0 $matrec                     !length bc and components of mag. vector, material index' //
     &      ''
           write(lunclc,'(a)') '$LxHalfMag $LyMag $LzMag                     !dimension'
           write(lunclc,'(a)') '$nHalfMagDivX $nMagDivY $nMagDivZ 1. 1.      !segmention'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '& Special_Magnet'
           write(lunclc,'(a)') 'Block mag56 Mag56 $e2colormag                !key, name, mother, color'
           write(lunclc,'(a)') '$x2D0URMagCen $yUMagCen $zMagCen             !position'
           write(lunclc,'(a)') '$Br 0.0 1.0 0.0 $matrec                      !length bc and components of mag. vector, material index' //
     &      ''
           write(lunclc,'(a)') '$LxMag $LyMag $LzMag                         !dimension'
           write(lunclc,'(a)') '$nMagDivX $nMagDivY $nMagDivZ 1. 1.          !segmention'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '& Special_Magnet'
           write(lunclc,'(a)') 'Block mag57 HMag57 $e2colormag               !key, name, mother, color'
           write(lunclc,'(a)') '$x3D0URMagCen $yUMagCen $zMagCen             !position'
           write(lunclc,'(a)') '$Br 1.0 0.0 0.0 $matrec                      !length bc and components of mag. vector, material index' //
     &      ''
           write(lunclc,'(a)') '$LxHalfMag $LyMag $LzMag                     !dimension'
           write(lunclc,'(a)') '$nHalfMagDivX $nMagDivY $nMagDivZ 1. 1.      !segmention'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '& Special_Magnet'
           write(lunclc,'(a)') 'Block mag58 HMag58 $e1colormag               !key, name, mother, color'
           write(lunclc,'(a)') '$x1D1URMagCen $yUMagCen $zMagCen             !position'
           write(lunclc,'(a)') '$E1Br -1.0 0.0 0.0 $matrec                   !length bc and components of mag. vector, material index' //
     &      ''
           write(lunclc,'(a)') '$LxHalfMag $LyMag $LzMag                     !dimension'
           write(lunclc,'(a)') '$nHalfMagDivX $nMagDivY $nMagDivZ 1. 1.      !segmention'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '& Special_Magnet'
           write(lunclc,'(a)') 'Block mag59 Mag59 $e1colormag                !key, name, mother, color'
           write(lunclc,'(a)') '$x2D1URMagCen $yUMagCen $zMagCen             !position'
           write(lunclc,'(a)') '$E1Br 0.0 1.0 0.0 $matrec                    !length bc and components of mag. vector, material index' //
     &      ''
           write(lunclc,'(a)') '$LxMag $LyMag $LzMag                         !dimension'
           write(lunclc,'(a)') '$nMagDivX $nMagDivY $nMagDivZ 1. 1.          !segmention'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '& Special_Magnet'
           write(lunclc,'(a)') 'Block mag60 HMag60 $e1colormag               !key, name, mother, color'
           write(lunclc,'(a)') '$x3D1URMagCen $yUMagCen $zMagCen             !position'
           write(lunclc,'(a)') '$E1Br 1.0 0.0 0.0 $matrec                    !length bc and components of mag. vector, material index' //
     &      ''
           write(lunclc,'(a)') '$LxHalfMag $LyMag $LzMag                     !dimension'
           write(lunclc,'(a)') '$nHalfMagDivX $nMagDivY $nMagDivZ 1. 1.      !segmention'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '* Lower left girder'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '& Special_Magnet'
           write(lunclc,'(a)') 'Block mag61 HMag61 $e2colormag               !key, name, mother, color'
           write(lunclc,'(a)') '$x1D0LLMagCen $yMagCen $zLLMagCen            !position'
           write(lunclc,'(a)') '$Br 1.0 0.0 0.0 $matrec                      !length bc and components of mag. vector, material index' //
     &      ''
           write(lunclc,'(a)') '$LxHalfMag $LyMag $LzMag                     !dimension'
           write(lunclc,'(a)') '$nHalfMagDivX $nMagDivY $nMagDivZ 1. 1.      !segmention'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '& Special_Magnet'
           write(lunclc,'(a)') 'Block mag62 Mag62 $e2colormag                !key, name, mother, color'
           write(lunclc,'(a)') '$x2D0LLMagCen $yMagCen $zLLMagCen            !position'
           write(lunclc,'(a)') '$Br 0.0 1.0 0.0 $matrec                      !length bc and components of mag. vector, material index' //
     &      ''
           write(lunclc,'(a)') '$LxMag $LyMag $LzMag                         !dimension'
           write(lunclc,'(a)') '$nMagDivX $nMagDivY $nMagDivZ 1. 1.          !segmention'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '& Special_Magnet'
           write(lunclc,'(a)') 'Block mag63 HMag63 $e2colormag               !key, name, mother, color'
           write(lunclc,'(a)') '$x3D0LLMagCen $yMagCen $zLLMagCen            !position'
           write(lunclc,'(a)') '$Br -1.0 0.0 0.0 $matrec                     !length bc and components of mag. vector, material index' //
     &      ''
           write(lunclc,'(a)') '$LxHalfMag $LyMag $LzMag                     !dimension'
           write(lunclc,'(a)') '$nHalfMagDivX $nMagDivY $nMagDivZ 1. 1.      !segmention'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '& Special_Magnet'
           write(lunclc,'(a)') 'Block mag64 HMag64 $e1colormag               !key, name, mother, color'
           write(lunclc,'(a)') '$x1D1LLMagCen $yMagCen $zLLMagCen            !position'
           write(lunclc,'(a)') '$E1Br 1.0 0.0 0.0 $matrec                    !length bc and components of mag. vector, material index' //
     &      ''
           write(lunclc,'(a)') '$LxHalfMag $LyMag $LzMag                     !dimension'
           write(lunclc,'(a)') '$nHalfMagDivX $nMagDivY $nMagDivZ 1. 1.      !segmention'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '& Special_Magnet'
           write(lunclc,'(a)') 'Block mag65 Mag65 $e1colormag                !key, name, mother, color'
           write(lunclc,'(a)') '$x2D1LLMagCen $yMagCen $zLLMagCen            !position'
           write(lunclc,'(a)') '$E1Br 0.0 1.0 0.0 $matrec                    !length bc and components of mag. vector, material index' //
     &      ''
           write(lunclc,'(a)') '$LxMag $LyMag $LzMag                         !dimension'
           write(lunclc,'(a)') '$nMagDivX $nMagDivY $nMagDivZ 1. 1.          !segmention'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '& Special_Magnet'
           write(lunclc,'(a)') 'Block mag66 HMag66 $e1colormag               !key, name, mother, color'
           write(lunclc,'(a)') '$x3D1LLMagCen $yMagCen $zLLMagCen            !position'
           write(lunclc,'(a)') '$E1Br -1.0 0.0 0.0 $matrec                   !length bc and components of mag. vector, material index' //
     &      ''
           write(lunclc,'(a)') '$LxHalfMag $LyMag $LzMag                     !dimension'
           write(lunclc,'(a)') '$nHalfMagDivX $nMagDivY $nMagDivZ 1. 1.      !segmention'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '* Upper left girder'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '& Special_Magnet'
           write(lunclc,'(a)') 'Block mag67 HMag67 $e2colormag               !key, name, mother, color'
           write(lunclc,'(a)') '$x1D0ULMagCen $yUMagCen $zULMagCen           !position'
           write(lunclc,'(a)') '$Br -1.0 0.0 0.0 $matrec                     !length bc and components of mag. vector, material index' //
     &      ''
           write(lunclc,'(a)') '$LxHalfMag $LyMag $LzMag                     !dimension'
           write(lunclc,'(a)') '$nHalfMagDivX $nMagDivY $nMagDivZ 1. 1.      !segmention'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '& Special_Magnet'
           write(lunclc,'(a)') 'Block mag68 Mag68 $e2colormag                !key, name, mother, color'
           write(lunclc,'(a)') '$x2D0ULMagCen $yUMagCen $zULMagCen           !position'
           write(lunclc,'(a)') '$Br 0.0 1.0 0.0 $matrec                      !length bc and components of mag. vector, material index' //
     &      ''
           write(lunclc,'(a)') '$LxMag $LyMag $LzMag                         !dimension'
           write(lunclc,'(a)') '$nMagDivX $nMagDivY $nMagDivZ 1. 1.          !segmention'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '& Special_Magnet'
           write(lunclc,'(a)') 'Block mag69 HMag69 $e2colormag               !key, name, mother, color'
           write(lunclc,'(a)') '$x3D0ULMagCen $yUMagCen $zULMagCen           !position'
           write(lunclc,'(a)') '$Br 1.0 0.0 0.0 $matrec                      !length bc and components of mag. vector, material index' //
     &      ''
           write(lunclc,'(a)') '$LxHalfMag $LyMag $LzMag                     !dimension'
           write(lunclc,'(a)') '$nHalfMagDivX $nMagDivY $nMagDivZ 1. 1.      !segmention'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '* Lower right girder'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '& Special_Magnet'
           write(lunclc,'(a)') 'Block mag70 HMag70 $e1colormag               !key, name, mother, color'
           write(lunclc,'(a)') '$x1D1LRMagCen $yUMagCen $zULMagCen           !position'
           write(lunclc,'(a)') '$E1Br -1.0 0.0 0.0 $matrec                   !length bc and components of mag. vector, material index' //
     &      ''
           write(lunclc,'(a)') '$LxHalfMag $LyMag $LzMag                     !dimension'
           write(lunclc,'(a)') '$nHalfMagDivX $nMagDivY $nMagDivZ 1. 1.      !segmention'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '& Special_Magnet'
           write(lunclc,'(a)') 'Block mag71 Mag71 $e1colormag                !key, name, mother, color'
           write(lunclc,'(a)') '$x2D1LRMagCen $yUMagCen $zULMagCen           !position'
           write(lunclc,'(a)') '$E1Br 0.0 1.0 0.0 $matrec                    !length bc and components of mag. vector, material index' //
     &      ''
           write(lunclc,'(a)') '$LxMag $LyMag $LzMag                         !dimension'
           write(lunclc,'(a)') '$nMagDivX $nMagDivY $nMagDivZ 1. 1.          !segmention'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '& Special_Magnet'
           write(lunclc,'(a)') 'Block mag72 HMag72 $e1colormag               !key, name, mother, color'
           write(lunclc,'(a)') '$x3D1LRMagCen $yUMagCen $zULMagCen           !position'
           write(lunclc,'(a)') '$E1Br 1.0 0.0 0.0 $matrec                    !length bc and components of mag. vector, material index' //
     &      ''
           write(lunclc,'(a)') '$LxHalfMag $LyMag $LzMag                     !dimension'
           write(lunclc,'(a)') '$nHalfMagDivX $nMagDivY $nMagDivZ 1. 1.      !segmention'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '& Special_Magnet'
           write(lunclc,'(a)') 'Block mag73 HMag73 $e2colormag               !key, name, mother, color'
           write(lunclc,'(a)') '$x1D2LRMagCen $yMagCen $zMagCen              !position'
           write(lunclc,'(a)') '$E2Br 1.0 0.0 0.0 $matrec                    !length bc and components of mag. vector, material index' //
     &      ''
           write(lunclc,'(a)') '$LxHalfMag $LyMag $LzMag                     !dimension'
           write(lunclc,'(a)') '$nHalfMagDivX $nMagDivY $nMagDivZ 1. 1.      !segmention'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '& Special_Magnet'
           write(lunclc,'(a)') 'Block mag74 Mag74 $e2colormag                !key, name, mother, color'
           write(lunclc,'(a)') '$x2D2LRMagCen $yMagCen $zMagCen              !position'
           write(lunclc,'(a)') '$E2Br 0.0 1.0 0.0 $matrec                    !length bc and components of mag. vector, material index' //
     &      ''
           write(lunclc,'(a)') '$LxMag $LyMag $LzMag                         !dimension'
           write(lunclc,'(a)') '$nMagDivX $nMagDivY $nMagDivZ 1. 1.          !segmention'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '& Special_Magnet'
           write(lunclc,'(a)') 'Block mag75 HMag75 $e2colormag               !key, name, mother, color'
           write(lunclc,'(a)') '$x3D2LRMagCen $yMagCen $zMagCen              !position'
           write(lunclc,'(a)') '$E2Br -1.0 0.0 0.0 $matrec                   !length bc and components of mag. vector, material index' //
     &      ''
           write(lunclc,'(a)') '$LxHalfMag $LyMag $LzMag                     !dimension'
           write(lunclc,'(a)') '$nHalfMagDivX $nMagDivY $nMagDivZ 1. 1.      !segmention'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '* Upper right girder'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '& Special_Magnet'
           write(lunclc,'(a)') 'Block mag76 HMag76 $e2colormag               !key, name, mother, color'
           write(lunclc,'(a)') '$x1D2URMagCen $yUMagCen $zMagCen             !position'
           write(lunclc,'(a)') '$E2Br -1.0 0.0 0.0 $matrec                   !length bc and components of mag. vector, material index' //
     &      ''
           write(lunclc,'(a)') '$LxHalfMag $LyMag $LzMag                     !dimension'
           write(lunclc,'(a)') '$nHalfMagDivX $nMagDivY $nMagDivZ 1. 1.      !segmention'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '& Special_Magnet'
           write(lunclc,'(a)') 'Block mag77 Mag77 $e2colormag                !key, name, mother, color'
           write(lunclc,'(a)') '$x2D2URMagCen $yUMagCen $zMagCen             !position'
           write(lunclc,'(a)') '$E2Br 0.0 1.0 0.0 $matrec                    !length bc and components of mag. vector, material index' //
     &      ''
           write(lunclc,'(a)') '$LxMag $LyMag $LzMag                         !dimension'
           write(lunclc,'(a)') '$nMagDivX $nMagDivY $nMagDivZ 1. 1.          !segmention'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '& Special_Magnet'
           write(lunclc,'(a)') 'Block mag78 HMag78 $e2colormag               !key, name, mother, color'
           write(lunclc,'(a)') '$x3D2URMagCen $yUMagCen $zMagCen             !position'
           write(lunclc,'(a)') '$E2Br 1.0 0.0 0.0 $matrec                    !length bc and components of mag. vector, material index' //
     &      ''
           write(lunclc,'(a)') '$LxHalfMag $LyMag $LzMag                     !dimension'
           write(lunclc,'(a)') '$nHalfMagDivX $nMagDivY $nMagDivZ 1. 1.      !segmention'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '* Lower left girder'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '& Special_Magnet'
           write(lunclc,'(a)') 'Block mag79 HMag79 $e2colormag               !key, name, mother, color'
           write(lunclc,'(a)') '$x1D2LLMagCen $yMagCen $zLLMagCen            !position'
           write(lunclc,'(a)') '$E2Br 1.0 0.0 0.0 $matrec                    !length bc and components of mag. vector, material index' //
     &      ''
           write(lunclc,'(a)') '$LxHalfMag $LyMag $LzMag                     !dimension'
           write(lunclc,'(a)') '$nHalfMagDivX $nMagDivY $nMagDivZ 1. 1.      !segmention'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '& Special_Magnet'
           write(lunclc,'(a)') 'Block mag80 Mag80 $e2colormag                !key, name, mother, color'
           write(lunclc,'(a)') '$x2D2LLMagCen $yMagCen $zLLMagCen            !position'
           write(lunclc,'(a)') '$E2Br 0.0 1.0 0.0 $matrec                    !length bc and components of mag. vector, material index' //
     &      ''
           write(lunclc,'(a)') '$LxMag $LyMag $LzMag                         !dimension'
           write(lunclc,'(a)') '$nMagDivX $nMagDivY $nMagDivZ 1. 1.          !segmention'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '& Special_Magnet'
           write(lunclc,'(a)') 'Block mag81 HMag81 $e2colormag               !key, name, mother, color'
           write(lunclc,'(a)') '$x3D2LLMagCen $yMagCen $zLLMagCen            !position'
           write(lunclc,'(a)') '$E2Br -1.0 0.0 0.0 $matrec                   !length bc and components of mag. vector, material index' //
     &      ''
           write(lunclc,'(a)') '$LxHalfMag $LyMag $LzMag                     !dimension'
           write(lunclc,'(a)') '$nHalfMagDivX $nMagDivY $nMagDivZ 1. 1.      !segmention'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '* Lower right girder'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '& Special_Magnet'
           write(lunclc,'(a)') 'Block mag82 HMag82 $e2colormag               !key, name, mother, color'
           write(lunclc,'(a)') '$x1D2LRMagCen $yUMagCen $zULMagCen           !position'
           write(lunclc,'(a)') '$E2Br -1.0 0.0 0.0 $matrec                   !length bc and components of mag. vector, material index' //
     &      ''
           write(lunclc,'(a)') '$LxHalfMag $LyMag $LzMag                     !dimension'
           write(lunclc,'(a)') '$nHalfMagDivX $nMagDivY $nMagDivZ 1. 1.      !segmention'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '& Special_Magnet'
           write(lunclc,'(a)') 'Block mag83 Mag83 $e2colormag                !key, name, mother, color'
           write(lunclc,'(a)') '$x2D2LRMagCen $yUMagCen $zULMagCen           !position'
           write(lunclc,'(a)') '$E2Br 0.0 1.0 0.0 $matrec                    !length bc and components of mag. vector, material index' //
     &      ''
           write(lunclc,'(a)') '$LxMag $LyMag $LzMag                         !dimension'
           write(lunclc,'(a)') '$nMagDivX $nMagDivY $nMagDivZ 1. 1.          !segmention'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '& Special_Magnet'
           write(lunclc,'(a)') 'Block mag84 HMag84 $e2colormag               !key, name, mother, color'
           write(lunclc,'(a)') '$x3D2LRMagCen $yUMagCen $zULMagCen           !position'
           write(lunclc,'(a)') '$E2Br 1.0 0.0 0.0 $matrec                    !length bc and components of mag. vector, material index' //
     &      ''
           write(lunclc,'(a)') '$LxHalfMag $LyMag $LzMag                     !dimension'
           write(lunclc,'(a)') '$nHalfMagDivX $nMagDivY $nMagDivZ 1. 1.      !segmention'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '& Module'
           write(lunclc,'(a)') '0. 0. 0. 1 1        !offset of module, number and number of associated module'
           write(lunclc,'(a)') '$nPeriods           !number of arrays within module'
           write(lunclc,'(a)') '$PerLen 1. 0. 0. 0. !spacing and direction of arrangement, rotation angle'
           write(lunclc,'(a)') '1. 1. 1.            !scaling and symmetry of magnetization vector'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') '& Materials'
           write(lunclc,'(a)') '1                       ! number of material files'
           write(lunclc,'(a)') '1 1 1 undumag_mu.dat    ! number, type, mode, and filename'
           write(lunclc,'(a)') ''
           write(lunclc,'(a)') ''


          flush(lunclc)
          close(lunclc)
          close(luntmp)

        else if (kbundumag.eq.4) then

          open(newunit=lunmat,file="undumag_mu.dat")
          write(lunmat,*)umupar,uksiper," ! mu_Par and ksi_Per"
          flush(lunmat)
          close(lunmat)

          if (uperlen_h.gt.0.0d0) then
            umaglx_h=uperlen_h/2.0d0-umagspac_h-upollx_h-upolspac_h
          else
            uperlen_h=2.0d0*(umaglx_h+umagspac_h+upollx_h+upolspac_h)
          endif

          open(newunit=lunclc,file="undumag.clc")

          write(lunclc,'(a)') "*! Lines written by WAVE"
          write(lunclc,'(a)') "*! Hybrid"
          write(lunclc,'(a)') " "
          write(lunclc,'(a)') "& User_Comment"
          write(cval,*)kwrun
          call util_string_trim(cval,k1,k2)
          write(lunclc,'(a)') trim(chwcom) // " (Run " // cval(k1:k2) // ")"
          write(lunclc,'(a)') " "

          write(lunclc,'(a)') " "

          write(cval,*)nperiod_h
          call util_string_trim(cval,k1,k2)
          write(lunclc,'(a)') "$nPeriods=" // cval(k1:k2)

          write(cval,'(g12.6)')undugap_h
          call util_string_trim(cval,k1,k2)
          write(lunclc,'(a)') "$FullGap=" // cval(k1:k2)

          write(lunclc,'(a)') " "
          write(cval,'(g12.6)')abs(umagbc_h)
          call util_string_trim(cval,k1,k2)
          write(lunclc,'(a)') "$Br=" // cval(k1:k2)
          write(cval,'(g12.6)')abs(umupar_h)
          call util_string_trim(cval,k1,k2)
          write(lunclc,'(a)') "$Mu=" // cval(k1:k2)
          write(cval,'(g12.6)')abs(uksiper_h)
          call util_string_trim(cval,k1,k2)
          write(lunclc,'(a)') "$KsiPerp=" // cval(k1:k2)
          write(lunclc,'(a)') "$MagMat=1"
          write(lunclc,'(a)') " "

          write(cval,'(g12.6)')umaglx_h
          call util_string_trim(cval,k1,k2)
          write(lunclc,'(a)') "$LxMag=" // cval(k1:k2)
          write(cval,'(g12.6)')umagly_h
          call util_string_trim(cval,k1,k2)
          write(lunclc,'(a)') "$LyMag=" // cval(k1:k2)
          write(cval,'(g12.6)')umaglz_h
          call util_string_trim(cval,k1,k2)
          write(lunclc,'(a)') "$LzMagFull=" // cval(k1:k2)
          write(cval,'(g12.6)')UMAGCH_H
          call util_string_trim(cval,k1,k2)
          write(lunclc,'(a)') "$ChamfM=" // cval(k1:k2)
          write(lunclc,'(a)') " "

          write(cval,'(g12.6)')UCOATING_H
          call util_string_trim(cval,k1,k2)
          write(lunclc,'(a)') "$MCoating=" // cval(k1:k2)
          write(lunclc,'(a)') " "

          write(cval,'(g12.6)')UMAGSPAC_H
          call util_string_trim(cval,k1,k2)
          write(lunclc,'(a)') "$AirGap=" // cval(k1:k2)
          write(lunclc,'(a)') " "

          write(cval,*)numdivx_h
          call util_string_trim(cval,k1,k2)
          write(lunclc,'(a)') "$nMagDivX=" // cval(k1:k2)
          write(cval,*)numdivy_h
          call util_string_trim(cval,k1,k2)
          write(lunclc,'(a)') "$nMagDivY=" // cval(k1:k2)
          write(cval,*)max(numdivz_h/2,1)
          call util_string_trim(cval,k1,k2)
          write(lunclc,'(a)') "$nMagDivZHalf=" // cval(k1:k2)

          write(lunclc,'(a)') " "

          write(lunclc,'(a)') "$IronMat=2"
          write(lunclc,'(a)') " "

          write(cval,'(g12.6)')upollx_h
          call util_string_trim(cval,k1,k2)
          write(lunclc,'(a)') "$LxPol=" // cval(k1:k2)
          write(cval,'(g12.6)')upolly_h
          call util_string_trim(cval,k1,k2)
          write(lunclc,'(a)') "$LyPol=" // cval(k1:k2)
          write(cval,'(g12.6)')upollz_h
          call util_string_trim(cval,k1,k2)
          write(lunclc,'(a)') "$LzPolFull=" // cval(k1:k2)
          write(cval,'(g12.6)')UPOLCH_H
          call util_string_trim(cval,k1,k2)
          write(lunclc,'(a)') "$ChamfP=" // cval(k1:k2)

          write(lunclc,'(a)') " "
          write(cval,'(g12.6)')UPOLSPAC_H
          call util_string_trim(cval,k1,k2)
          write(lunclc,'(a)') "$KeeperGap=" // cval(k1:k2) // "                 ! Be careful with ixsym, if keeper- and airgap are different"
          write(lunclc,'(a)') " "

          write(cval,*)nupdivx_h
          call util_string_trim(cval,k1,k2)
          write(lunclc,'(a)') "$nPolDivX=" // cval(k1:k2)
          if (nupdivx_h/2*2.eq.nupdivx_h) then
            write(cval,*)nupdivx_h/2
          else
            write(cval,*)nupdivx_h/2+1
          endif
          write(lunclc,'(a)') "$nPolDivXHalf=" // cval(k1:k2)
          write(cval,*)nupdivy_h
          call util_string_trim(cval,k1,k2)
          write(lunclc,'(a)') "$nPolDivY=" // cval(k1:k2)
          write(cval,*)max(nupdivz_h/2,1)
          call util_string_trim(cval,k1,k2)
          write(lunclc,'(a)') "$nPolDivZHalf=" // cval(k1:k2)
          write(lunclc,'(a)') " "
          write(cval,'(g12.6)')max(FracDivFe_h,1.0d0)
          call util_string_trim(cval,k1,k2)
          write(lunclc,'(a)') "$FracDivFeY=" // cval(k1:k2)
          write(cval,'(g12.6)')max(FracDivFeZ_h,1.0d0)
          call util_string_trim(cval,k1,k2)
          write(lunclc,'(a)') "$FracDivFeZ=" // cval(k1:k2)
          write(lunclc,'(a)') " "

          nendpol=0
          nendmag=0
          xend=0.0d0
          perlen=2.0d0*(upollx_h+umaglx_h+umagspac_h+upolspac_h)

          do i=1,nspec_h

            if (abs(msmag_h(i)).eq.1) then

              nendmag=nendmag+1
              write(cend,*)nendmag

              call util_string_trim(cend,k1end,k2end)

              cline="$LxEndMag" // cend(k1end:k2end) // "="
              call util_string_trim(cline,k1,k2)
              cline=cline(k1:k2)
              write(cval,'(g12.6)')usmaglx_h(i)
              call util_string_trim(cval,k1,k2)
              write(lunclc,'(a)')trim(cline) // cval(k1:k2)

              cline="$LyEndMag" // cend(k1end:k2end) // "="
              call util_string_trim(cline,k1,k2)
              cline=cline(k1:k2)
              write(cval,'(g12.6)')usmagly_h(i)
              call util_string_trim(cval,k1,k2)
              write(lunclc,'(a)')trim(cline) // cval(k1:k2)

              cline="$LzEndMagFull" // cend(k1end:k2end) // "="
              call util_string_trim(cline,k1,k2)
              cline=cline(k1:k2)
              write(cval,'(g12.6)')usmaglz_h(i)
              call util_string_trim(cval,k1,k2)
              write(lunclc,'(a)')trim(cline) // cval(k1:k2)

              cline="$ChamfEndMag" // cend(k1end:k2end) // "="
              call util_string_trim(cline,k1,k2)
              cline=cline(k1:k2)
              write(cval,'(g12.6)')usmagch_h(i)
              call util_string_trim(cval,k1,k2)
              write(lunclc,'(a)')trim(cline) // cval(k1:k2)

              cline="$SpacerEndMag" // cend(k1end:k2end) // "="
              call util_string_trim(cline,k1,k2)
              cline=cline(k1:k2)
              write(cval,'(g12.6)')usmagspac_h(i)
              call util_string_trim(cval,k1,k2)
              write(lunclc,'(a)')trim(cline) // cval(k1:k2)

              cline="$YoffsetEndMag" // cend(k1end:k2end) // "="
              call util_string_trim(cline,k1,k2)
              cline=cline(k1:k2)
              write(cval,'(g12.6)')usmagdy_h(i)
              call util_string_trim(cval,k1,k2)
              write(lunclc,'(a)')trim(cline) // cval(k1:k2)

            else

              nendpol=nendpol+1
              write(cval,*)nendpol
              call util_string_trim(cend,k1end,k2end)

              cline="$LxEndPol" // cend(k1end:k2end) // "="
              call util_string_trim(cline,k1,k2)
              cline=cline(k1:k2)
              write(cval,'(g12.6)')usmaglx_h(i)
              call util_string_trim(cval,k1,k2)
              write(lunclc,'(a)')trim(cline) // cval(k1:k2)

              cline="$LyEndPol" // cend(k1end:k2end) // "="
              call util_string_trim(cline,k1,k2)
              cline=cline(k1:k2)
              write(cval,'(g12.6)')usmagly_h(i)
              call util_string_trim(cval,k1,k2)
              write(lunclc,'(a)')trim(cline) // cval(k1:k2)

              cline="$LzEndPolFull" // cend(k1end:k2end) // "="
              call util_string_trim(cline,k1,k2)
              cline=cline(k1:k2)
              write(cval,'(g12.6)')usmaglz_h(i)
              call util_string_trim(cval,k1,k2)
              write(lunclc,'(a)')trim(cline) // cval(k1:k2)

              cline="$ChamfEndPol" // cend(k1end:k2end) // "="
              call util_string_trim(cline,k1,k2)
              cline=cline(k1:k2)
              write(cval,'(g12.6)')usmagch_h(i)
              call util_string_trim(cval,k1,k2)
              write(lunclc,'(a)')trim(cline) // cval(k1:k2)

              cline="$SpacerEndPol" // cend(k1end:k2end) // "="
              call util_string_trim(cline,k1,k2)
              cline=cline(k1:k2)
              write(cval,'(g12.6)')usmagspac_h(i)
              call util_string_trim(cval,k1,k2)
              write(lunclc,'(a)')trim(cline) // cval(k1:k2)

              cline="$YoffsetEndPol" // cend(k1end:k2end) // "="
              call util_string_trim(cline,k1,k2)
              cline=cline(k1:k2)
              write(cval,'(g12.6)')usmagdy_h(i)
              call util_string_trim(cval,k1,k2)
              write(lunclc,'(a)')trim(cline) // cval(k1:k2)

            endif

            write(lunclc,'(a)') " "

          enddo !nspec_h

          write(lunclc,'(a)') " "

          write(lunclc,'(a)') "$Mcol=2"
          write(lunclc,'(a)') "$Pcol=4"
          write(lunclc,'(a)') " "
          write(lunclc,'(a)') "$PerLen = 2. * ( $KeeperGap + $LxPol + $AirGap + $LxMag ) "
          write(lunclc,'(a)') " "
          write(lunclc,'(a)') "$LzMag = $LzMagFull / 2."
          write(lunclc,'(a)') "$LzPol = $LzPolFull / 2."
          write(lunclc,'(a)') " "
          write(lunclc,'(a)') "$xShift  = ( $nPeriods - 1 ) * $PerLen"
          write(lunclc,'(a)') " "
          write(lunclc,'(a)') "$xMag1 = - $xShift - $LxPol / 2. - $AirGap - $LxMag / 2."
          write(lunclc,'(a)') "$xMag2 = $xMag1 - $PerLen / 2."
          write(lunclc,'(a)') " "
          write(lunclc,'(a)') "$yMag = - $FullGap / 2. - $LyMag / 2."
          write(lunclc,'(a)') "$zMag  = - $LzMag / 2."
          write(lunclc,'(a)') " "
          write(lunclc,'(a)') "$xPol1 = - $xShift -  $PerLen / 2."
          write(lunclc,'(a)') "$xPol2 = $xPol1 - $PerLen / 2."
          write(lunclc,'(a)') " "
          write(lunclc,'(a)') "$yPol = - $FullGap / 2. - $LyPol / 2."
          write(lunclc,'(a)') "$zPol  = - $LzPol / 2."
          write(lunclc,'(a)') " "
          write(lunclc,'(a)') " "
          write(lunclc,'(a)') "! Central pole, see under special magnets"
          write(lunclc,'(a)') "$LxHalfPol = $LxPol / 2."
          write(lunclc,'(a)') "$xHPol1 = - $LxPol / 4."
          write(lunclc,'(a)') " "

          nendpol=0
          nendmag=0
          xend=-dble(nperiod_h)*perlen-upollx_h/2.0d0-umagspac_h
          xend=xend-usmaglx_h(1)-usmagspac_h(1)

          write(lunclc,'(a)')
     &      "$xSpec0  = - $nPeriods * $PerLen - $LxPol / 2. - $AirGap"
          cindo='0'
          kio=1
          keo=1

          do i=1,nspec_h

            write(cind,'(I4)') i
            call util_string_trim(cind,indi,inde)

            xend=xend-usmaglx_h(i)-usmagspac_h(i)
            x=xend+usmaglx_h(i)/2.0d0

            if (abs(msmag_h(i)).eq.1) then

              nendmag=nendmag+1
              write(cmag,'(I4)') nendmag
              call util_string_trim(cmag,kim,kem)

              cline=
     &          "$xSpec" // cind(indi:inde) // " = $xSpec" // cindo(kio:keo) //
     &          " - $LxEndMag" // cmag(kim:kem) //
     &          " - $SpacerEndMag" // cmag(kim:kem)
              write(lunclc,'(a)')trim(cline)

              cline=
     &          "$xEndMag" // cmag(kim:kem) // " = $xSpec" // cind(indi:inde) //
     &          " + $LxEndMag" // cmag(kim:kem) // " / 2. "
              write(lunclc,'(a)')trim(cline)

              cline="$yEndMag" // cmag(kim:kem) // " = $yMag + $YoffsetEndMag" // cmag(kim:kem)
              write(lunclc,'(a)')trim(cline)

            else

              nendpol=nendpol+1
              write(cpol,'(I4)') nendpol
              call util_string_trim(cpol,kip,kep)

              cline=
     &          "$xSpec" // cind(indi:inde) // " = $xSpec" // cindo(kio:keo) //
     &          " - $LxEndPol" // cpol(kip:kep) // " - $SpacerEndPol" // cpol(kip:kep)
              write(lunclc,'(a)')trim(cline)

              cline=
     &          "$xEndPol" // cpol(kip:kep) // " = $xSpec" // cind(indi:inde) //
     &          " + $LxEndPol" // cpol(kip:kep) // " / 2. "
              write(lunclc,'(a)')trim(cline)

              cline="$yEndPol" // cpol(kip:kep) // " = $yPol + $YoffsetEndPol" // cpol(kip:kep)
              write(lunclc,'(a)')trim(cline)

            endif

            write(lunclc,'(a)') " "

            cindo=cind
            kio=indi
            keo=inde

          enddo !nspec_h

          if (kbundumag.eq.4) then

            if (ixsym_h.ne.0) then
              write(lunclc,'(a)') "$ixsym=1 ! overwrites value of undumag.nam"
            else
              write(lunclc,'(a)') "$ixsym=0 ! overwrites value of undumag.nam"
            endif

            write(lunclc,'(a)')
     &        "$kxsym=int[$ixsym/($ixsym-0.0001)]"
            write(lunclc,'(a)')
     &        "$IronMatSym=$IronMat*abs[$kxsym]"

            if (iysym_h.ne.0) then
              write(lunclc,'(a)') "$iysym=1 ! overwrites value of undumag.nam"
            else
              write(lunclc,'(a)') "$iysym=0 ! overwrites value of undumag.nam"
            endif

            if (izsym_h.ne.0) then
              write(lunclc,'(a)') "$izsym=1 ! overwrites value of undumag.nam"
            else
              write(lunclc,'(a)') "$izsym=0 ! overwrites value of undumag.nam"
            endif

          endif

          write(lunclc,'(a)') " "
          write(lunclc,'(a)') "*-------------------------------------"
          write(lunclc,'(a)') " "
          write(lunclc,'(a)') "! Section of magnetic items"
          write(lunclc,'(a)') " "
          write(lunclc,'(a)') " "

          write(lunclc,'(a)') "!   /---\"
          write(lunclc,'(a)') "!   |   |"
          write(lunclc,'(a)') "!   | > |    First main magnet"
          write(lunclc,'(a)') "!   |   |"
          write(lunclc,'(a)') "!   -----"
          write(lunclc,'(a)') " "

          write(lunclc,'(a)') "& Magnet"
          write(lunclc,'(a)') "BlockChamf mag1 Mag1 $Mcol                  !key, name, mother, color "
          write(lunclc,'(a)') "$xMag1 $yMag $zMag                          !position of magnet"
          write(lunclc,'(a)') "$Br 1.0 0.0 0.0 $MagMat                     !length bc and components of mag. vector, material index"
          write(lunclc,'(a)') "$LxMag $LyMag $LzMag $ChamfM                !dimensions"
          write(lunclc,'(a)') "$nMagDivX $nMagDivY $nMagDivZHalf 1. 1.     !segmentation"
          write(lunclc,'(a)') " "

          write(lunclc,'(a)') " "
          write(lunclc,'(a)') "!   /---\"
          write(lunclc,'(a)') "!   |   |"
          write(lunclc,'(a)') "!   |   |    First full pole"
          write(lunclc,'(a)') "!   |   |"
          write(lunclc,'(a)') "!   -----"
          write(lunclc,'(a)') " "

          write(lunclc,'(a)') "& Pole"
          write(lunclc,'(a)') "BlockChamf pol1 Pol1 $Pcol                                !key, name, mother, color "
          write(lunclc,'(a)') "$xPol1 $yPol $zPol                                        !position of magnet"
          write(lunclc,'(a)') "$IronMat                                                  !material index"
          write(lunclc,'(a)') "$LxPol $LyPol $LzPol $ChamfP                              !dimensions"
          write(lunclc,'(a)') "$nPolDivX $nPolDivY $nPolDivZHalf $FracDivFeY $FracDivFeZ !segmentation"
          write(lunclc,'(a)') " "

          write(lunclc,'(a)') "!   /---\"
          write(lunclc,'(a)') "!   |   |"
          write(lunclc,'(a)') "!   | < |    Second main magnet"
          write(lunclc,'(a)') "!   |   |"
          write(lunclc,'(a)') "!   -----"
          write(lunclc,'(a)') " "

          write(lunclc,'(a)') "& Magnet"
          write(lunclc,'(a)') "BlockChamf mag2 Mag2 $Mcol                  !key, name, mother, color "
          write(lunclc,'(a)') "$xMag2 $yMag $zMag                          !position of magnet"
          write(lunclc,'(a)') "$Br -1.0 0.0 0.0 $MagMat                    !length bc and components of mag. vector, material index"
          write(lunclc,'(a)') "$LxMag $LyMag $LzMag $ChamfM                !dimensions"
          write(lunclc,'(a)') "$nMagDivX $nMagDivY $nMagDivZHalf 1. 1.     !segmentation"
          write(lunclc,'(a)') " "

          write(lunclc,'(a)') "!   /---\"
          write(lunclc,'(a)') "!   |   |"
          write(lunclc,'(a)') "!   |   |    Second full pole"
          write(lunclc,'(a)') "!   |   |"
          write(lunclc,'(a)') "!   -----"
          write(lunclc,'(a)') " "

          write(lunclc,'(a)') "& Pole"
          write(lunclc,'(a)') "BlockChamf pol2 Pol2 $Pcol                                !key, name, mother, color "
          write(lunclc,'(a)') "$xPol2 $yPol $zPol                                        !position of magnet"
          write(lunclc,'(a)') "$IronMat                                                  !material index"
          write(lunclc,'(a)') "$LxPol $LyPol $LzPol $ChamfP                              !dimensions"
          write(lunclc,'(a)') "$nPolDivX $nPolDivY $nPolDivZHalf $FracDivFeY $FracDivFeZ !segmentation"
          write(lunclc,'(a)') " "
          write(lunclc,'(a)') " "

          write(lunclc,'(a)') "*------------------------------------------------------------------------"
          write(lunclc,'(a)') "& Module"
          write(lunclc,'(a)') "0. 0. 0.            !offset of module"
          write(lunclc,'(a)') "$nPeriods           !number of arrays within module"
          write(lunclc,'(a)') "$PerLen 1. 0. 0. 0. !spacing and direction of arrangement, rotation angle"
          write(lunclc,'(a)') " 1. 1. 1.           !scaling of magnetization vector"
          write(lunclc,'(a)') "*------------------------------------------------------------------------"
          write(lunclc,'(a)') " "
          write(lunclc,'(a)') " "

          write(lunclc,'(a)') "!   /---|"
          write(lunclc,'(a)') "!   |   |"
          write(lunclc,'(a)') "!   |   |    First half pole in the center (see $ixsym)"
          write(lunclc,'(a)') "!   |   |"
          write(lunclc,'(a)') "!   -----"
          write(lunclc,'(a)') " "
          write(lunclc,'(a)') "& Special_Pole"
          write(lunclc,'(a)') "BlockUsChamf hpol1 HPol1 $Pcol                                !key, name, mother, color "
          write(lunclc,'(a)') "$xHPol1 $yPol $zPol                                           !position of magnet"
          write(lunclc,'(a)') "$IronMatSym                                                   !material index"
          write(lunclc,'(a)') "$LxHalfPol $LyPol $LzPol $ChamfP                              !dimensions"
          write(lunclc,'(a)') "$nPolDivXHalf $nPolDivY $nPolDivZHalf $FracDivFeY $FracDivFeZ !segmentation"
          write(lunclc,'(a)') " "

          nendmag=0
          nendpol=0

          do i=1,nspec_h

            if (abs(msmag_h(i)).eq.1) then

              nendmag=nendmag+1
              write(cend,*)nendmag
              call util_string_trim(cend,k1end,k2end)

              write(lunclc,'(a)') " "
              write(lunclc,'(a)') "!   /---\"
              write(lunclc,'(a)') "!   |   |"
              if (nendmag/2*2.eq.nendmag) then
                write(lunclc,'(a)') "!   | < |    End magnet " // cend(k1end:k2end)
              else
                write(lunclc,'(a)') "!   | > |    End magnet " // cend(k1end:k2end)
              endif
              write(lunclc,'(a)') "!   |   |"
              write(lunclc,'(a)') "!   -----"

              write(lunclc,'(a)') " "
              write(lunclc,'(a)') "& Special_Magnet"
              write(lunclc,'(a)') "BlockChamf" //
     &          " emag" // cend(k1end:k2end) //
     &          " Emag" // cend(k1end:k2end) // " $Mcol" //
     &          " !key, name, mother, color "

              write(lunclc,'(a)') "$xEndMag" // cend(k1end:k2end) //
     &          " $yEndMag" // cend(k1end:k2end) //
     &          " $zMag" // "    !position of magnet"

              if (nendmag/2*2.eq.nendmag) then
                write(lunclc,'(a)') "$Br -1. 0. 0. $MagMat                            !length bc and components of mag. vector, material index"
              else
                write(lunclc,'(a)') "$Br 1. 0. 0. $MagMat                             !length bc and components of mag. vector, material index"
              endif

              write(lunclc,'(a)') "$LxEndMag" // cend(k1end:k2end) //
     &          " $LyEndMag" // cend(k1end:k2end) //
     &          " $LzMag $ChamfM    !dimensions"

              write(lunclc,'(a)') "$nMagDivX $nMagDivY $nMagDivZHalf 1. 1.          !segmentation"

            else

              nendpol=nendpol+1
              write(cend,*)nendpol
              call util_string_trim(cend,k1end,k2end)

              write(lunclc,'(a)') " "
              write(lunclc,'(a)') " "
              write(lunclc,'(a)') "!   /---\"
              write(lunclc,'(a)') "!   |   |"
              write(lunclc,'(a)') "!   |   |    End pole " // cend(k1end:k2end)
              write(lunclc,'(a)') "!   |   |"
              write(lunclc,'(a)') "!   -----"

              write(lunclc,'(a)') " "
              write(lunclc,'(a)') "& Special_Pole"
              write(lunclc,'(a)') "BlockChamf" //
     &          " epol" // cend(k1end:k2end) //
     &          " Epol" // cend(k1end:k2end) // " $Pcol" //
     &          "                    !key, name, mother, color "

              write(lunclc,'(a)') "$xEndPol" // cend(k1end:k2end) //
     &          " $yEndPol" // cend(k1end:k2end) //
     &          " $zPol" // "                                                                !position of pole"

              write(lunclc,'(a)')
     &          "$IronMat                                                                    !material index"

              write(lunclc,'(a)') "$LxEndPol" // cend(k1end:k2end) //
     &          " $LyEndPol" // cend(k1end:k2end) //
     &          " $LzPol $ChamfP                                                             !dimensions"

              write(lunclc,'(a)') "$nPolDivX $nPolDivY $nPolDivZHalf $FracDivFeY $FracDivFeZ !segmentation"

            endif

            write(lunclc,'(a)') " "

          enddo !nspec_h


          write(lunclc,'(a)') "*------------------------------------------------------------------------"
          write(lunclc,'(a)') "& Materials"
          write(lunclc,'(a)') "2 ! number of materials"
          write(lunclc,'(a)') "1 1 1 undumag_mu.dat               ! number, type, mode, and filename"
          write(lunclc,'(a)') "2 2 3 " // trim(CHUNDUIRON_H) // " ! number, type, mode, and filename"
          write(lunclc,'(a)') "*------------------------------------------------------------------------"

88        flush(lunclc)
          close(lunclc)

        endif !(kbundumag.eq.4)

      endif !(kbundumag.lt.0) then

c      stop "Ende in run_undumag"
      if (kbunduverb_c.eq.0) then
        cline=trim(chundutmp) // " > undumag.log"
      else
        cline=trim(chundutmp)
      endif

      print*,''
      print*, "     --- Spawning UNDUMAG run ---"
      print*,''
      print*,"      ",trim(cline)
      print*,''

      istat=system(trim(cline))

      if (istat.ne.0) then
        print*,"*** Error in run_undumag: Bad return status, check undumag.log ***"
        stop "*** Program WAVE aborted ***"
      endif

      open(newunit=lund,file="undumag.stat")
      read(lund,*)istat
      if (istat.ne.0) then
        print*,''
        print*,"     *** UNDUMAG has crashed ***"
        print*,''
        stop "*** Program WAVE aborted ***"
      else
        print*,''
        print*,"     --- UNDUMAG has finished ---"
        print*,''
      endif
      close(lund)

      kbundumag=1

      return

 90   continue
      write(lungfo,*)"*** Error in run_undumag: File"
      write(lungfo,*)trim(chundunam)
      write(lungfo,*)"not found, Check the namelist UNDUMAGN in wave.in ***"
      write(lungfo,*)"*** Program WAVE aborted ***"
      write(6,*)"*** Error in run_undumag: File"
      write(6,*)trim(chundunam)
      write(6,*)"not found, Check the namelist UNDUMAGN in wave.in ***"
      write(6,*)"*** Program WAVE aborted ***"
      stop
 91   continue
      write(lungfo,*)"*** Error in run_undumag: File"
      write(lungfo,*)trim(chunduclc)
      write(lungfo,*)"not found, Check the namelist UNDUMAGN in wave.in ***"
      write(lungfo,*)"*** Program WAVE aborted ***"
      write(6,*)" Error in run_undumag: File"
      write(6,*)trim(chunduclc)
      write(6,*)"not found, Check the namelist UNDUMAGN in wave.in ***"
      write(6,*)"*** Program WAVE aborted ***"
      stop
 92   continue
      write(lungfo,*)"*** Error in run_undumag: File"
      write(lungfo,*)trim(chundutmp)
      write(lungfo,*)"not found, Check the namelist UNDUMAGN in wave.in ***"
      write(lungfo,*)"*** Program WAVE aborted ***"
      write(6,*)"*** Error in run_undumag: File"
      write(6,*)trim(chundutmp)
      write(6,*)"not found, Check the namelist UNDUMAGN in wave.in ***"
      write(6,*)"*** Program WAVE aborted ***"
      stop

      end
