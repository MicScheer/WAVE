*CMZ :          24/01/2025  17.09.59  by  Michael Scheer
*CMZ :  1.03/02 22/09/2016  17.00.47  by  Michael Scheer
*CMZ :  1.03/01 10/10/2014  08.51.19  by  Michael Scheer
*CMZ :  1.02/00 01/10/2014  11.55.05  by  Michael Scheer
*CMZ :  1.01/02 27/09/2014  16.22.46  by  Michael Scheer
*CMZ :  1.01/00 24/09/2014  14.15.23  by  Michael Scheer
*CMZ :  1.00/01 23/09/2014  22.42.58  by  Michael Scheer
*CMZ :  0.01/03 23/09/2014  13.59.29  by  Michael Scheer
*CMZ :  0.01/02 22/09/2014  17.08.09  by  Michael Scheer
*CMZ :  0.01/01 26/08/2014  14.55.55  by  Michael Scheer
*CMZ :  0.01/00 25/08/2014  15.21.06  by  Michael Scheer
*CMZ :  0.00/06 22/08/2014  15.08.05  by  Michael Scheer
*CMZ :  0.00/05 18/08/2014  14.47.51  by  Michael Scheer
*CMZ :  0.00/04 18/08/2014  10.08.48  by  Michael Scheer
*CMZ :  0.00/03 06/08/2014  09.01.11  by  Michael Scheer
*CMZ :  0.00/02 10/07/2014  16.15.45  by  Michael Scheer
*CMZ :  0.00/01 07/07/2014  09.10.33  by  Michael Scheer
*CMZ : 00.00/02 04/07/2014  21.01.32  by  Michael Scheer
*-- Author :    Michael Scheer   27/06/2014
      subroutine mshplt_init(idev,xpapsiz,ypapsiz,ibxl,ibyb,ibxr,ibyt,
     &  file,viewer,viewerkill,rescale)

      implicit none

      real xpapsiz,ypapsiz,rescale
      integer lun,i,ierr,ibxl,ibxr,ibyb,ibyt,idev
      character(*) file,viewer,viewerkill

*KEEP,MSHPLTINCL.
      include 'mshplt.cmn'
*KEEP,MPLOT.

      real rescale_mshplt

      character(2048)
     &  fileeps_mshplt,
     &  viewer_mshplt,viewer_kill_mshplt

      common/mplotc/
     &  rescale_mshplt,
     &  fileeps_mshplt,
     &  viewer_mshplt,viewer_kill_mshplt
*KEND.

      logical lisopen

      lunbase_ps=100000
1     inquire(unit=lunbase_ps,opened=lisopen)
      if (lisopen) then
        lunbase_ps=lunbase_ps+1000
        if (lunbase_ps.gt.10000000) then
          write(6,*)'*** Error in mshplt_init: No free LUN base found ***'
          stop
        endif
        goto 1
      endif

      lun=lunbase_ps+1

      if (xpapsiz.lt.0..or.ypapsiz.lt.0) then
        ihigzmode_ps=1
      else
        ihigzmode_ps=0
      endif

      if (rescale.le.0) then
        rescale_mshplt=1.
      else
        rescale_mshplt=rescale
      endif

      irunviewer_ps=1

      if (len_trim(viewer).eq.0.or.idev.eq.0) then
        irunviewer_ps=0
      else
        viewer_ps=viewer(1:len_trim(viewer))
        viewer_kill_mshplt=viewerkill(1:len_trim(viewerkill))
      endif

      if (irunviewer_ps.ne.0.and.idev.ne.0) then
        iviewinter_ps=1
      else
        iviewinter_ps=1
      endif

      istat_ps=0
      nfile_ps=0
      inolabs_ps=0

      radtodeg_ps=180./(4.*atan(1.))

      xsizorig_ps=abs(xpapsiz) !cm
      ysizorig_ps=abs(ypapsiz) !cm
      if (xsizorig_ps.le.0) xsizorig_ps=15.
      if (ysizorig_ps.le.0) ysizorig_ps=15.
      xleftorig_ps=4.
      ybottomorig_ps=4.
      xrightorig_ps=xleftorig_ps+xsizorig_ps !cm
      ytoporig_ps=ybottomorig_ps+ysizorig_ps !cm

      kbbxl_ps=ibxl
      kbbxr_ps=ibxr
      kbbyb_ps=ibyb
      kbbyt_ps=ibyt

      if (
     &    kbbxl_ps.lt.0.or.kbbxr_ps.lt.0.or.
     &    kbbyb_ps.lt.0.or.kbbyt_ps.lt.0) then
        kbbxl_ps=int(50. - 50.*(xsizorig_ps-10.)/10.)
        kbbxr_ps=int(450. + 350.*(xsizorig_ps-10.)/10.)
        kbbyb_ps=kbbxl_ps
        kbbyt_ps=kbbxr_ps
      endif

      if (
     &    kbbxl_ps.lt.0.or.kbbxr_ps.lt.0.or.
     &    kbbyb_ps.lt.0.or.kbbyt_ps.lt.0) then
        kbbxl_ps=0
        kbbxr_ps=400
        kbbyb_ps=0
        kbbyt_ps=400
      endif

      open(unit=99,file='.mshplt.bb',status='unknown')
      write(99,*) kbbxl_ps,kbbxr_ps,kbbyb_ps,kbbyt_ps
      close(99)

      mode3d_ps=0

      idrawgtit_ps=0 !to control drawing of global title

      xsiz_ps=xsizorig_ps !cm
      ysiz_ps=ysizorig_ps !cm

      if (isscale_ps.ne.1) scale_ps=pttocm_ps

      xleft_ps=xleftorig_ps !cm
      xright_ps=xleft_ps+xsizorig_ps !cm
      ybottom_ps=ybottomorig_ps !cm
      ytop_ps=ybottom_ps+ysizorig_ps !cm

      x_ps=xleft_ps
      y_ps=ybottom_ps

      ang_ps=0. !degree
      tang_ps=0. !degree
      pi_ps=4.*atan(1.)

      wxmin_ps=0.
      wxmax_ps=xsiz_ps
      wymin_ps=0.
      wymax_ps=ysiz_ps

      scalex_ps=xsiz_ps/(wxmax_ps-wxmin_ps)
      scaley_ps=ysiz_ps/(wymax_ps-wymin_ps)

      ilabmod_ps=0

      scaletxt_ps=ysizorig_ps/20.*1.5
      if (ihigzmode_ps.ne.0) then
        scaletxt_ps=ysizorig_ps/20.
      endif
      call mshplt_set_character_height(0.3) !cm
      call mshplt_set_index_height(0.3*0.8) !cm
      ticsiz_ps=0.2*scaletxt_ps !cm

      offgtitx_ps=-2. !cm
      offgtity_ps=4.*chhe_ps

      if (isoffdate_ps.ne.1) then
        offdatex_ps=-10*chhe_ps
        offdatey_ps=1.5*chhe_ps
      endif

      rlinewidth_ps=0.03*scaletxt_ps !cm

      xlaboff_ps=2.*chhe_ps !offset of x-axis label cm
      ylaboff_ps=2.*chhe_ps !offset of y-axis label cm
      zlaboff_ps=2.5*chhe_ps !offset of z-axis label cm

      xoffexp_ps=-2.0*chhe_ps! long. offset of power term
      yoffexp_ps=ylaboff_ps+1.5*chhe_ps ! transv. offset of power

      xtitoff_ps=4.0*chhe_ps !offset of x-axis title cm
      ytitoff_ps=4.0*chhe_ps !offset of y-axis title cm
      ztitoff_ps=5.*chhe_ps !offset of z-axis title cm

      ygti_ps=1.5*chhe_ps !cm
      gsiz_ps=2.*chhe_ps !cm

      inewpage_ps=2

      ibuffpos_ps=0

      do i=1,9
          write(chch_ps(i),'(a,i1,a)')'\',i,'  '
      enddo
      do i=10,99
        write(chch_ps(i),'(a,i2,a)')'\',i,' '
      enddo
      do i=100,999
        write(chch_ps(i),'(a,i3)')'\',i
      enddo

      ! Font Symbol
      rmtyp_ps(1,0)=456.!dot
      rmtyp_ps(2,0)=-0.04
      rmtyp_ps(1,1)=267. !bullet
      rmtyp_ps(2,1)=-0.33
      rmtyp_ps(1,2)=250. !volle Raute
      rmtyp_ps(2,2)=-0.25
      rmtyp_ps(1,3)=264. !diag. cross
      rmtyp_ps(2,3)=-0.265
      rmtyp_ps(1,4)=340. !hohle Raute
      rmtyp_ps(2,4)=-0.37
      rmtyp_ps(1,5)=304. !Circle with diag. cross
      rmtyp_ps(2,5)=-0.35
      rmtyp_ps(1,6)=305. !Circle with cross
      rmtyp_ps(2,6)=-0.35
      rmtyp_ps(1,6)=452. !Asterix
      rmtyp_ps(2,6)=-0.33
      rmtyp_ps(1,7)=453. !Plus
      rmtyp_ps(2,7)=-0.27
      rmtyp_ps(1,8)=43. !Hash
      rmtyp_ps(2,8)=-0.27
      rmtyp_ps(1,9)=-1 !circle, drawn by mshplt_circle
      rmtyp_ps(2,9)=0.

      !variables for hplot zones

      if (ihigzmode_ps.ne.0) then
        nxzone_ps=1
        nyzone_ps=1
        nzone_ps=nxzone_ps*nyzone_ps
      endif

      xmgl_ps=0.5/scaletxt_ps !cm
      xmgr_ps=0.5/scaletxt_ps !cm
      ymgl_ps=2.0/scaletxt_ps !cm
      ymgu_ps=2.0/scaletxt_ps !cm
      xwin_ps=2.5/scaletxt_ps !cm
      ywin_ps=2.5/scaletxt_ps !cm

      xtit_ps=''
      ytit_ps=''
      ztit_ps=''
      gtit_ps=''

      tsiz_ps=chhe_ps

      if (isDate_ps.ne.1) kDate_ps=1
      if (isBox_ps.ne.1) kBox_ps=0

      if (nfile_ps.eq.0) then
        call mshplt_file_open(lun,file,'',ierr)
        if (ierr.ne.0) goto 9999
      endif

      call mshplt_set_marker_type(1)
      call mshplt_set_marker_size(0.4) !cm

      call mshplt_set_theta_phi(60.,30.)

      call mshplt_set_color(1,0,0,0)
      call mshplt_set_fill_color(3,0,0,0)
      call mshplt_set_line_color(1,0,0,0)
      call mshplt_set_text_color(1,0,0,0)
      call mshplt_set_frame_color(1,0,0,0)
      call mshplt_set_marker_color(1,0,0,0)

      call mshplt_set_label_size(0.3)

      write(lun_ps,'(a)')'% end of mshplt_init'

      return

 9999 istat_ps=-1

      return
      end
*CMZ :          01/08/2018  08.38.20  by  Michael Scheer
*CMZ :  1.03/02 22/09/2016  13.19.58  by  Michael Scheer
*CMZ :  1.02/00 03/10/2014  12.49.08  by  Michael Scheer
*CMZ :  1.01/02 25/09/2014  09.32.57  by  Michael Scheer
*CMZ :  1.01/00 24/09/2014  14.33.11  by  Michael Scheer
*CMZ :  1.00/01 24/09/2014  12.18.04  by  Michael Scheer
*CMZ :  0.01/02 18/09/2014  13.01.02  by  Michael Scheer
*CMZ :  0.01/01 27/08/2014  15.30.42  by  Michael Scheer
*CMZ :  0.01/00 25/08/2014  10.34.58  by  Michael Scheer
*CMZ :  0.00/04 14/08/2014  13.19.56  by  Michael Scheer
*CMZ :  0.00/03 06/08/2014  09.00.02  by  Michael Scheer
*CMZ :  0.00/02 09/07/2014  13.35.02  by  Michael Scheer
*CMZ :  0.00/01 07/07/2014  09.10.33  by  Michael Scheer
*CMZ : 00.00/02 30/06/2014  10.35.42  by  Michael Scheer
*-- Author :    Michael Scheer   27/06/2014
      subroutine mshplt_end

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEEP,mplot.

      real rescale_mshplt

      character(2048)
     &  fileeps_mshplt,
     &  viewer_mshplt,viewer_kill_mshplt

      common/mplotc/
     &  rescale_mshplt,
     &  fileeps_mshplt,
     &  viewer_mshplt,viewer_kill_mshplt
*KEND.

      integer nitemp
      parameter(nitemp=100)

      integer ifile,lun,lunps,ipos(2,nitemp),nwords,istat,leneps
      real scalex,scaley

      character(2048) fileps

      call mshplt_flush_buff
      write(lun_ps,*)'% begin of mshplt_end'
      write(lun_ps,*)'showpage'

c      lunps=lunall_ps(nfile_ps)+1
c      call util_get_free_lun(lunps)

      do ifile=1,nfile_ps
        if (rescale_mshplt.ne.1..and.rescale_mshplt.gt.0.0) then
          leneps=len_trim(fileall_ps(ifile))
          fileps=fileall_ps(ifile)(1:leneps-4)//
     &      '_scaled'//fileall_ps(ifile)(leneps-3:leneps)
          open(newunit=lunps,file=fileps(1:len_trim(fileps)),status='unknown')
          lun=lunall_ps(ifile)
          if (lun.gt.0) then
            write(lunall_ps(ifile),*)'showpage'
            rewind(lun)
1           read(lun,'(a)',end=9)cline_ps
            call util_string_split(cline_ps,nitemp,nwords,ipos,istat)
            if (cline_ps(ipos(1,1):ipos(2,1)).eq.'%%Title') then
              write(lunps,'(a)')'%% Title '//fileps(1:len_trim(fileps))
              goto 1
            endif
            if (nwords.ge.3) then
              if (cline_ps(ipos(1,3):ipos(2,3)).eq.'scale') then
                read(cline_ps,*)scalex,scaley
                write(cline_ps,*)scalex*rescale_mshplt,scaley*rescale_mshplt,
     &            ' scale'
                write(lunps,'(a)')cline_ps(1:len_trim(cline_ps))
                goto 11
              endif
            endif
            write(lunps,'(a)')cline_ps(1:len_trim(cline_ps))
            goto 1
11          read(lun,'(a)',end=9)cline_ps
            write(lunps,'(a)')cline_ps(1:len_trim(cline_ps))
            goto 11
9           close(lunps)
c            close(lunall_ps(ifile))
            lunall_ps(ifile)=-lunall_ps(ifile)
          endif
        else
          close(lunall_ps(ifile))
          lunall_ps(ifile)=-lunall_ps(ifile)
        endif
      enddo

      call sleep(1)

      if(irunviewer_ps.ne.0) then

        if(iviewinter_ps.ne.0) then
          write(6,*)'Hit any key to terminate viewer:'
          read(5,'(a)')
          call system(viewer_kill_mshplt(1:len_trim(viewer_kill_mshplt)))
        endif

      endif

      nfile_ps=0
      kzone_ps=0
      itouched_ps=0

      return
      end
*CMZ :          04/03/2023  17.35.33  by  Michael Scheer
*CMZ :  1.03/02 25/04/2016  12.19.28  by  Michael Scheer
*CMZ :  1.03/01 10/10/2014  13.29.05  by  Michael Scheer
*CMZ :  1.02/00 03/10/2014  11.06.23  by  Michael Scheer
*CMZ :  1.01/00 24/09/2014  16.33.15  by  Michael Scheer
*CMZ :  0.01/02 22/09/2014  20.47.10  by  Michael Scheer
*CMZ :  0.01/01 27/08/2014  08.54.49  by  Michael Scheer
*CMZ :  0.00/06 22/08/2014  15.47.40  by  Michael Scheer
*CMZ :  0.00/05 18/08/2014  13.15.26  by  Michael Scheer
*CMZ :  0.00/04 14/08/2014  12.58.41  by  Michael Scheer
*CMZ :  0.00/03 05/08/2014  15.57.52  by  Michael Scheer
*CMZ :  0.00/02 09/07/2014  14.42.14  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_frame(xmini,xmaxi,ymini,ymaxi,xtit,ytit,chopt)

c Draw frame of axis, ticksize is ticsiz_ps
c chopt controls with axis are to be drawn:
c        L: Left axis with labels
c        l: Left axis without labels
c        R: Right axis with labels
c        r: Right axis without labels
c        B: Bottom axis with labels
c        b: Bottom axis without labels
c        T: Top axis with labels
c        t: Top axis without labels
c
c Option chopt:
c       L: Left axis with labels
c       l: Left axis without labels
c       B: Bottom axis with labels
c       b: Bottom axis without labels
c       R: Right axis with labels
c       r: right axis without labels
c       T: Top axis with labels
c       t: top axis without labels
c       c: Clipping mode is on

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      real xmin,xmax,ymin,ymax,wxmn,wxmx,wymn,wymx,x(4),y(4),
     &  xmini,xmaxi,ymini,ymaxi,
     &  xminin,xmaxin,yminin,ymaxin,
     &  theo,phio
      real sizlab(100),xlabrel(100),ylaboffset(100),anglabrel(100),ticlen(100),
     &  ticposrel(100),titoff,titangrel,ticangrel(100),titsiz,titposrel

      integer il,ir,ib,it,i,iclip,ix,iy,ifirst
      integer kred,kgreen,kblue,kcolor,nlab,ntic
      integer kredo,kgreeno,kblueo,kcoloro

      character(*) xtit,ytit,chopt
      character(2048) copt
      character c1
      character(12) chlab(100)
      call mshplt_flush_buff

      write(lun_ps,'(a)')'% begin of mshplt_frame'

      if (idrawgtit_ps.eq.1) call mshplt_draw_title

      theo=theta_ps
      phio=phi_ps

      if (len_trim(chopt).eq.0) then
        copt='LBrtc'
      else
        copt=chopt
      endif

      isameframe_ps=0
      iclip=0
      il=0
      ir=0
      ib=0
      it=0
      inolabs_ps=1

      kcolor=kFramecolor_ps
      kgreen=kFramegreen_ps
      kred=kFramered_ps
      kblue=kFrameblue_ps

      call mshplt_set_frame_color(kcolor,kgreen,kred,kblue)
      call mshplt_get_line_color(kcoloro,kgreeno,kredo,kblueo)
      call mshplt_set_line_color(kcolor,kgreen,kred,kblue)

      do i=1,len_trim(copt)

        c1=copt(i:i)

        !Left axis
        if (c1.eq.'L') then
          il=1
          if (chhe_ps.eq.0.0) il=-1
        else if (c1.eq.'l') then
          il=-1

        !Right axis
        else if (c1.eq.'R') then
          ir=1
          if (chhe_ps.eq.0.0) ir=-1
        else if (c1.eq.'r') then
          ir=-1

        !Bottom axis
        else if (c1.eq.'B') then
          ib=1
          if (chhe_ps.eq.0.0) ib=-1
        else if (c1.eq.'b') then
          ib=-1

        !Top axis
        else if (c1.eq.'T') then
          it=1
          if (chhe_ps.eq.0.0) it=-1
        else if (c1.eq.'t') then
          it=-1
        else if (c1.eq.'s'.or.c1.eq.'S') then
          isameframe_ps=1
        else if (c1.eq.'c') then
          iclip=1
        endif
      enddo

      if (ihigzmode_ps.ne.0.and.isameframe_ps.eq.0.and.kzone_ps.gt.0) then
        ix=nxzone_ps
        iy=nyzone_ps
        if (kzone_ps.lt.nzone_ps) then
          ifirst=kzone_ps+1
          call mshplt_zone(ix,iy,ifirst,'s')
          else
            if (itouched_ps.ne.0) call mshplt_newpage
          call mshplt_zone(ix,iy,1,' ')
        endif
      endif

      !convert to pad system
      wxmn=xleft_ps
      wxmx=xright_ps
      wymn=ybottom_ps
      wymx=ytop_ps

      xminin=xmini
      xmaxin=xmaxi
      yminin=ymini
      ymaxin=ymaxi

      if (xminin.ne.xminin) then
          xminin=-1.0e30
      endif

      if (xminin.ne.0.0) then
        if (1.0/xminin.eq.0.0) then
          if (xminin.lt.0.0) then
            xminin=-1.0e30
          else
            xminin=1.0e30
          endif
        endif
      endif

      if (xmaxin.ne.xmaxin) then
          xmaxin=1.0e30
      endif

      if (xmaxin.ne.0.0) then
        if (1.0/xmaxin.eq.0.0) then
          if (xmaxin.lt.0.0) then
            xmaxin=-1.0e30
          else
            xmaxin=1.0e30
          endif
        endif
      endif

      if (xminin.eq.xmaxin) then
        xmin=xminin-0.5
        xmax=xmaxin+0.5
      endif

      if (yminin.ne.yminin) then
          yminin=-1.0e30
      endif
      if (yminin.ne.0.0) then
        if (1.0/yminin.eq.0.0) then
          if (yminin.lt.0.0) then
            yminin=-1.0e30
          else
            yminin=1.0e30
          endif
        endif
      endif

      if (ymaxin.ne.ymaxin) then
          ymaxin=1.0e30
      endif
      if (ymaxin.ne.0.0) then
        if (1.0/ymaxin.eq.0.0) then
          if (ymaxin.lt.0.0) then
            ymaxin=-1.0e30
          else
            ymaxin=1.0e30
          endif
        endif
      endif

      xmin=xminin
      xmax=xmaxin

      ymin=yminin
      ymax=ymaxin

      if (yminin.eq.ymaxin) then
        ymin=yminin-0.5
        ymax=ymaxin+0.5
      endif

      if (log10y_ps.ne.0) then
        if (yminin.le.0.0) then
          ymin=-30.
        else
          ymin=float(int(alog10(yminin)))
          if (ymin.le.0.0.and.
     &      (10.**ymin-yminin)/10.**ymin.gt.1.0e-5) ymin=ymin-1.
        endif
        if (ymaxin.le.yminin) then
          ymax=ymin+1
        else
          ymax=float(int(alog10(ymax)))
          if (ymax.ge.0.0.and.
     &      (ymaxin-10.**ymax)/10.**ymax.gt.1.0e-5)
     &      ymax=ymax+1.
        endif
      endif

      wymin_ps=ymin
      wymax_ps=ymax
      scaley_ps=ysiz_ps/(wymax_ps-wymin_ps)

      if (log10x_ps.ne.0) then
        if (xminin.le.0.0) then
          xmin=-30.
        else
          xmin=float(int(alog10(xminin)))
          if (xmin.le.0.0.and.
     &      (10.**xmin-xminin)/10.**xmin.gt.1.0e-5) xmin=xmin-1.
        endif
        if (xmaxin.le.xminin) then
          xmax=xmin+1
        else
          xmax=float(int(alog10(xmax)))
          if (xmax.ge.0.0.and.
     &      (xmaxin-10.**xmax)/10.**xmax.gt.1.0e-5)
     &      xmax=xmax+1.
        endif
      endif

      wxmin_ps=xmin
      wxmax_ps=xmax
      scalex_ps=xsiz_ps/(wxmax_ps-wxmin_ps)

      write(lun_ps,'(a)')'% frame: xleft_ps xright_ps ybottom_ps ytop_ps:'
      write(lun_ps,*)'%',xleft_ps,xright_ps,ybottom_ps,ytop_ps
      write(lun_ps,'(a)')cline_ps(2:len_trim(cline_ps))
      write(lun_ps,'(a)')'% size: xsiz_ps ysiz_ps:'
      write(lun_ps,*)'%',xsiz_ps,ysiz_ps
      write(lun_ps,'(a)')cline_ps(2:len_trim(cline_ps))

      write(lun_ps,*)'% coord.-system: wxmin_ps wxmax_ps:',
     &  wxmin_ps,wxmax_ps
      write(lun_ps,*)'% coord.-system: wymin_ps wymax_ps:',
     &  wymin_ps,wymax_ps
      write(lun_ps,'(a)')cline_ps(2:len_trim(cline_ps))

      if (il.ne.0) then

        if (il.lt.0) then
          inolabs_ps=1
        else
          inolabs_ps=0
        endif

        if (log10y_ps.eq.0) then
          nlab=-100
          sizlab=aLabHeight_ps*(1-inolabs_ps)
          ylaboffset=-ylaboff_ps
          anglabrel=-90.
          ticlen=-ticsiz_ps
          ticangrel=90.
          titsiz=tsiz_ps*(1-inolabs_ps)
          titoff=-ytitoff_ps
          titangrel=0.
          titposrel=0.5
          call mshplt_axis_taylor(
     &      xmin,xmin,ymin,ymax,ymin,ymax,
     &      nlab,chlab,sizlab,xlabrel,ylaboffset,
     &      xoffexp_ps,yoffexp_ps,
     &      anglabrel,
     &      ntic,ticlen,ticposrel,ticangrel,
     &      titsiz,titposrel,titoff,titangrel,ytit
     &      )
c          call mshplt_world_axis(wxmn,wxmn,wymn,wymx,ymin,ymax,ytit,
c     &      -1,0,0.,ytitoff_ps*scaletxt_ps,ylaboff_ps*scaletxt_ps)
        else
          call mshplt_world_log_axis(wxmn,wxmn,wymn,wymx,
     &      10.**ymin,10.**ymax,ytit,
     &      -1,-(1-iNoLabs_ps),0.,ytitoff_ps*scaletxt_ps,
     &      ylaboff_ps*scaletxt_ps)
        endif

      endif !il

      if (ir.ne.0) then

        if (ir.lt.0) then
          inolabs_ps=1
        else
          inolabs_ps=0
        endif

        if (log10y_ps.eq.0) then
          nlab=-100
          sizlab=aLabHeight_ps*(1-inolabs_ps)
          ylaboffset=ylaboff_ps
          anglabrel=-90.
          ticlen=ticsiz_ps
          ticangrel=90.
          titsiz=tsiz_ps*(1-inolabs_ps)
          titoff=ytitoff_ps
          titangrel=0.
          titposrel=0.5
          call mshplt_axis_taylor(
     &      xmax,xmax,ymin,ymax,yminin,ymaxin,
     &      nlab,chlab,sizlab,xlabrel,ylaboffset,
     &      xoffexp_ps,yoffexp_ps,
     &      anglabrel,
     &      ntic,ticlen,ticposrel,ticangrel,
     &      titsiz,titposrel,titoff,titangrel,ytit
     &      )
c          call mshplt_world_axis(wxmx,wxmx,wymn,wymx,ymin,ymax,ytit,
c     &      1,0,0.,ytitoff_ps*scaletxt_ps,ylaboff_ps*scaletxt_ps)
        else
          call mshplt_world_log_axis(
     &      wxmx,wxmx,wymn,wymx,10.**ymin,10.**ymax,ytit,
     &      1,(1-iNoLabs_ps),0.,
     &      ytitoff_ps*scaletxt_ps,ylaboff_ps*scaletxt_ps)
        endif

      endif !ir

      if (ib.ne.0) then

        if (ib.lt.0) then
          inolabs_ps=1
        else
          inolabs_ps=0
        endif

        if (log10x_ps.eq.0) then
          nlab=-100
          sizlab=aLabHeight_ps*(1-inolabs_ps)
          ylaboffset=xlaboff_ps
          anglabrel=0.
          ticlen=ticsiz_ps
          ticangrel=90.
          titsiz=tsiz_ps*(1-inolabs_ps)
          titoff=ytitoff_ps
          titangrel=0.
          titposrel=0.5
          call mshplt_axis_taylor(
     &      xmin,xmax,ymin,ymin,xminin,xmaxin,
     &      nlab,chlab,sizlab,xlabrel,ylaboffset,
     &      xoffexp_ps,yoffexp_ps,
     &      anglabrel,
     &      ntic,ticlen,ticposrel,ticangrel,
     &      titsiz,titposrel,titoff,titangrel,xtit
     &      )
c          call mshplt_world_axis(wxmn,wxmx,wymn,wymn,xmin,xmax,xtit,
c     &      1,1,0.,-xtitoff_ps*scaletxt_ps,-xlaboff_ps*scaletxt_ps)
        else
          call mshplt_world_log_axis(
     &      wxmn,wxmx,wymn,wymn,10.**xmin,10.**xmax,xtit,
     &      -1,-(1-iNoLabs_ps),
     &      0.,xtitoff_ps*scaletxt_ps,xlaboff_ps*scaletxt_ps)
        endif
      endif !ib

      if (it.ne.0) then
        if (it.lt.0) then
          inolabs_ps=1
        else
          inolabs_ps=0
        endif

        if (log10x_ps.eq.0) then
          nlab=-100
          sizlab=aLabHeight_ps*(1-inolabs_ps)
          ylaboffset=-xlaboff_ps
          anglabrel=0.
          ticlen=-ticsiz_ps
          ticangrel=90.
          titsiz=tsiz_ps*(1-inolabs_ps)
          titoff=-xtitoff_ps
          titangrel=0.
          titposrel=0.5
          call mshplt_axis_taylor(
     &      xmin,xmax,ymax,ymax,xminin,xmaxin,
     &      nlab,chlab,sizlab,xlabrel,ylaboffset,
     &      xoffexp_ps,yoffexp_ps,
     &      anglabrel,
     &      ntic,ticlen,ticposrel,ticangrel,
     &      titsiz,titposrel,titoff,titangrel,xtit
     &      )
c          call mshplt_world_axis(wxmn,wxmx,wymx,wymx,xmin,xmax,xtit,
c     &      -1,0,0.,-xtitoff_ps*scaletxt_ps,xlaboff_ps*scaletxt_ps)
        else
          call mshplt_world_log_axis(
     &      wxmn,wxmx,wymx,wymx,10.**xmin,10.**xmax,xtit,
     &      1,(1-iNoLabs_ps),0.,
     &      (xtitoff_ps-chhe_ps)*scaletxt_ps,xlaboff_ps*scaletxt_ps)
        endif
      endif

      x(1)=xmin
      y(1)=ymin
      x(2)=xmax
      y(2)=y(1)
      x(3)=x(2)
      y(3)=ymax
      x(4)=x(1)
      y(4)=y(3)
      if (iclip.eq.1) then
        call mshplt_clip(4,x,y)
      else
        call mshplt_reset_clipping
      endif

      if (kzone_ps.lt.0) then
        kzone_ps=ifirst_ps
      endif

c      iColor_ps=kcolor
c      ired_ps=kred
c      igreen_ps=kgreen
c      iblue_ps=kblue

      call mshplt_set_line_color(kcoloro,kgreeno,kredo,kblueo)
      call mshplt_set_theta_phi(theo,phio) ! due to axis routines

      write(lun_ps,'(a)')'% end of mshplt_frame'

      return
      end
*CMZ :  0.01/02 05/09/2014  15.41.43  by  Michael Scheer
*CMZ :  0.00/04 11/08/2014  17.53.18  by  Michael Scheer
*CMZ :  0.00/02 07/07/2014  12.22.09  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_body

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      write(cline_ps,*)
      call mshplt_fill_buff(cline_ps)

      return
      end
*CMZ :  0.01/02 05/09/2014  15.41.43  by  Michael Scheer
*CMZ :  0.00/03 04/08/2014  13.25.01  by  Michael Scheer
*CMZ :  0.00/02 09/07/2014  14.48.35  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_line(x1,y1,x2,y2)

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      real x(2),y(2),x1,y1,x2,y2

      x(1)=x1
      y(1)=y1
      x(2)=x2
      y(2)=y2

      call mshplt_pline(2,x,y)

      return
      end
*CMZ :  0.01/03 23/09/2014  09.36.37  by  Michael Scheer
*CMZ :  0.01/02 05/09/2014  15.41.43  by  Michael Scheer
*CMZ :  0.00/04 14/08/2014  12.58.41  by  Michael Scheer
*CMZ :  0.00/02 09/07/2014  15.30.29  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_spline(n,x,y)

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      real x(*),y(*),x1,y1,x2,y2
      integer n,i,nsegments
      character(2048) cline

      double precision, dimension (:), allocatable :: sd,xd,yd,x2d,y2d
      double precision xsd,ysd,dsd,ssd

      data nsegments/1000/

      if (n.le.1) return

      if (n.eq.2) then

        if (log10x_ps.eq.0) then
          x1=xleft_ps+scalex_ps*(x(1)-wxmin_ps)
        else
          x1=xleft_ps+scalex_ps*(alog10(x(1))-wxmin_ps)
        endif

        if (log10y_ps.eq.0) then
          y1=ybottom_ps+scaley_ps*(y(1)-wymin_ps)
        else
          y1=ybottom_ps+scaley_ps*(alog10(y(1))-wymin_ps)
        endif

        if (log10x_ps.eq.0) then
          x2=xleft_ps+scalex_ps*(x(2)-wxmin_ps)
        else
          x2=xleft_ps+scalex_ps*(alog10(x(2))-wxmin_ps)
        endif

        if (log10y_ps.eq.0) then
          y2=ybottom_ps+scaley_ps*(y(2)-wymin_ps)
        else
          y2=ybottom_ps+scaley_ps*(alog10(y(2))-wymin_ps)
        endif

        write(cline,*)x1,y1,x2,y2,' lineto'
        call mshplt_fill_buff(cline)
        return
      endif

      allocate(sd(n))
      allocate(xd(n))
      allocate(yd(n))
      allocate(x2d(n))
      allocate(y2d(n))

      do i=1,n
        sd(i)=i
        if (log10x_ps.eq.0) then
          xd(i)=xleft_ps+scalex_ps*(x(i)-wxmin_ps)
        else
          xd(i)=xleft_ps+scalex_ps*(alog10(x(i))-wxmin_ps)
        endif
        if (log10y_ps.eq.0) then
          yd(i)=ybottom_ps+scaley_ps*(y(i)-wymin_ps)
        else
          yd(i)=ybottom_ps+scaley_ps*(alog10(y(i))-wymin_ps)
        endif
      enddo
      call util_spline_interpolation_f90(n,sd,xd,sd(1),xd(1),x2d,-1) !init spline
      call util_spline_interpolation_f90(n,sd,yd,sd(1),yd(1),y2d,-1) !init spline

      dsd=dble(n-1)/dble(nsegments-1)
      ssd=sd(1)
      x1=sngl(xd(1))
      y1=sngl(yd(1))
      write(cline,*)x1,y1,' moveto'
      call mshplt_fill_buff(cline)
      do i=1,nsegments
        call util_spline_interpolation_f90(n,sd,xd,ssd,xsd,x2d,0)
        x2=sngl(xsd)
        call util_spline_interpolation_f90(n,sd,yd,ssd,ysd,y2d,0)
        y2=sngl(ysd)
        write(cline,*)x1,y1,x2,y2,' lineto'
        call mshplt_fill_buff(cline)
        x1=x2
        y1=y2
        ssd=min(ssd+dsd,sd(n))
      enddo

      write(cline,*)'stroke'
      call mshplt_fill_buff(cline)

      deallocate(sd)
      deallocate(xd)
      deallocate(yd)
      deallocate(x2d)
      deallocate(y2d)

      itouched_ps=1

      return
      end
*CMZ :  0.01/02 22/09/2014  14.18.27  by  Michael Scheer
*CMZ :  0.01/00 22/08/2014  18.11.33  by  Michael Scheer
*CMZ :  0.00/02 09/07/2014  13.34.47  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_fill_buff(cline)

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      character(*) cline

      if(ibuffpos_ps.lt.nbuffsize_ps) then
        ibuffpos_ps=ibuffpos_ps+1
        chbuff_ps(ibuffpos_ps)=cline(1:len_trim(cline))
      else
        call mshplt_flush_buff
        if (len_trim(cline).gt.0) then
          ibuffpos_ps=ibuffpos_ps+1
          chbuff_ps(ibuffpos_ps)=cline(1:len_trim(cline))
        endif
      endif

c      kzone_ps=abs(kzone_ps) !to indicate that pad is touched

      return
      end
*CMZ :          10/10/2024  10.57.55  by  Michael Scheer
*CMZ :  0.01/02 11/09/2014  12.59.54  by  Michael Scheer
*CMZ :  0.00/02 09/07/2014  13.35.55  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_flush_buff

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      integer i

c      write(lun_ps,'(a)')'% begin of mshplt_flush_buff'

      do i=1,ibuffpos_ps
        write(lun_ps,'(a)') chbuff_ps(i)(1:len_trim(chbuff_ps(i)))
      enddo
      flush(lun_ps)
      ibuffpos_ps=0

c      write(lun_ps,'(a)')'% end of mshplt_flush_buff'

      return
      end
*CMZ :          06/08/2018  15.49.27  by  Michael Scheer
*CMZ :  1.03/01 08/10/2014  14.16.19  by  Michael Scheer
*CMZ :  1.03/00 06/10/2014  16.06.01  by  Michael Scheer
*CMZ :  1.02/01 05/10/2014  10.00.08  by  Michael Scheer
*CMZ :  1.02/00 03/10/2014  11.12.46  by  Michael Scheer
*CMZ :  1.01/02 26/09/2014  13.31.57  by  Michael Scheer
*CMZ :  1.01/00 25/09/2014  08.49.12  by  Michael Scheer
*CMZ :  1.00/01 24/09/2014  13.15.20  by  Michael Scheer
*CMZ :  0.01/03 23/09/2014  13.56.57  by  Michael Scheer
*CMZ :  0.01/02 05/09/2014  15.41.43  by  Michael Scheer
*CMZ :  0.00/06 21/08/2014  16.57.22  by  Michael Scheer
*CMZ :  0.00/04 14/08/2014  12.58.41  by  Michael Scheer
*CMZ :  0.00/03 06/08/2014  11.22.28  by  Michael Scheer
*CMZ :  0.00/02 09/07/2014  16.03.13  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_frame3d(xminin,xmaxin,yminin,ymaxin,zminin,zmaxin,
     &  xtit,ytit,ztit,chopt)

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.


      real xmin,xmax,ymin,ymax,zmin,zmax,
     &  xminin,xmaxin,yminin,ymaxin,zminin,zmaxin,
     &  xcorn(8),ycorn(8),zcorn(8),ang,xyplen,
     &  ticheight,chhe,titang,rellabang

      integer ilinestyleo,ilinecoloro,iaxis,ibox,ic,iticside,ilabside
      integer kred,kgreen,kblue,kcolor,ifirst
      integer kredo,kgreeno,kblueo

      character(*) xtit,ytit,ztit,chopt
      character(2048) xtitd,ytitd,ztitd,choptd

      data xcorn/-0.5,0.5,0.5,-0.5,-0.5,0.5,0.5,-0.5/
      data ycorn/-0.5,-0.5,0.5,0.5,-0.5,-0.5,0.5,0.5/
      data zcorn/-0.5,-0.5,-0.5,-0.5,0.5,0.5,0.5,0.5/

      xmax=xmaxin
      xmin=xminin
      ymax=ymaxin
      ymin=yminin
      zmax=zmaxin
      zmin=zminin

      if (log10x_ps.ne.0) then

        if (xminin.le.0.0) then
          xmin=-30.
        else
          xmin=float(int(alog10(xminin)))
          if (xmin.le.0.0.and.
     &      (10.**xmin-xminin)/10.**xmin.gt.1.0e-5) xmin=xmin-1.
        endif

        if (xmaxin.le.xminin) then
          xmax=xmin+1
        else
          xmax=float(int(alog10(xmax)))
          if (xmax.ge.0.0.and.
     &      (xmaxin-10.**xmax)/10.**xmax.gt.1.0e-5)
     &      xmax=xmax+1.
        endif
      endif

      if (log10y_ps.ne.0) then

        if (yminin.le.0.0) then
          ymin=-30.
        else
          ymin=float(int(alog10(yminin)))
          if (ymin.le.0.0.and.
     &      (10.**ymin-yminin)/10.**ymin.gt.1.0e-5) ymin=ymin-1.
        endif

        if (ymaxin.le.yminin) then
          ymax=ymin+1
        else
          ymax=float(int(alog10(ymax)))
          if (ymax.ge.0.0.and.
     &      (ymaxin-10.**ymax)/10.**ymax.gt.1.0e-5)
     &      ymax=ymax+1.
        endif
      endif

      if (log10z_ps.ne.0) then

        if (zminin.le.0.0) then
          zmin=-30.
        else
          zmin=float(int(alog10(zminin)))
          if (zmin.le.0.0.and.
     &      (10.**zmin-zminin)/10.**zmin.gt.1.0e-5) zmin=zmin-1.
        endif

        if (zmaxin.le.zminin) then
          zmax=zmin+1
        else
          zmax=float(int(alog10(zmax)))
          if (zmax.ge.0.0.and.
     &      (zmaxin-10.**zmax)/10.**zmax.gt.1.0e-5)
     &      zmax=zmax+1.
        endif
      endif

      ticheight=ticsiz_ps
      chhe=chhe_ps

      write(lun_ps,'(a)')'% begin of mshplt_frame3d'

      if (ihigzmode_ps.ne.0.and.isameframe_ps.eq.0.and.kzone_ps.gt.0) then
        if (kzone_ps.lt.nzone_ps) then
          ifirst=kzone_ps+1
          call mshplt_zone(nxzone_ps,nyzone_ps,ifirst,'s')
          else
            if (itouched_ps.ne.0) call mshplt_newpage
          call mshplt_zone(nxzone_ps,nyzone_ps,1,' ')
        endif
      endif

      if (idrawgtit_ps.eq.1) call mshplt_draw_title

      kcolor=kFramecolor_ps
      kgreen=kFramegreen_ps
      kred=kFramered_ps
      kblue=kFrameblue_ps

      call mshplt_set_frame_color(kcolor,kgreen,kred,kblue)

      if (nzone_ps.le.0) then
        xcorn(1)=-0.5
        xcorn(2)=0.5
        xcorn(3)=0.5
        xcorn(4)=-0.5
        xcorn(5)=-0.5
        xcorn(6)=0.5
        xcorn(7)=0.5
        xcorn(8)=-0.5
        ycorn(1)=-0.5
        ycorn(2)=-0.5
        ycorn(3)=0.5
        ycorn(4)=0.5
        ycorn(5)=-0.5
        ycorn(6)=-0.5
        ycorn(7)=0.5
        ycorn(8)=0.5
        zcorn(1)=-0.5
        zcorn(2)=-0.5
        zcorn(3)=-0.5
        zcorn(4)=-0.5
        zcorn(5)=0.5
        zcorn(6)=0.5
        zcorn(7)=0.5
        zcorn(8)=0.5
      else
        xcorn(1)=-0.625
        xcorn(2)=0.625
        xcorn(3)=0.625
        xcorn(4)=-0.625
        xcorn(5)=-0.625
        xcorn(6)=0.625
        xcorn(7)=0.625
        xcorn(8)=-0.625
        ycorn(1)=-0.625
        ycorn(2)=-0.625
        ycorn(3)=0.625
        ycorn(4)=0.625
        ycorn(5)=-0.625
        ycorn(6)=-0.625
        ycorn(7)=0.625
        ycorn(8)=0.625
        zcorn(1)=-0.625
        zcorn(2)=-0.625
        zcorn(3)=-0.625
        zcorn(4)=-0.625
        zcorn(5)=0.625
        zcorn(6)=0.625
        zcorn(7)=0.625
        zcorn(8)=0.625
      endif

      xtitd=xtit
      ytitd=ytit
      ztitd=ztit
      choptd=chopt

      if(chopt.eq.'') choptd='AB'

      iaxis=0
      ibox=0
      ifbox_ps=0
      do ic=1,len_trim(choptd)
        if (choptd(ic:ic).eq.'A'.or.choptd(ic:ic).eq.'a') iaxis=1
        if (choptd(ic:ic).eq.'B'.or.choptd(ic:ic).eq.'b') then
          ibox=1
          ifbox_ps=1
        endif
      enddo

      ilinestyleo=ilinestyle_ps
      call mshplt_get_line_color(ilinecoloro,kredo,kblueo,kgreeno)

      call mshplt_3dto2d(8,xcorn,ycorn,zcorn,xpcorn_ps,ypcorn_ps)

      wxmin_ps=min(xpcorn_ps(1),xpcorn_ps(2),xpcorn_ps(3),xpcorn_ps(4),
     &  xpcorn_ps(5),xpcorn_ps(6),xpcorn_ps(7),xpcorn_ps(8))
      wxmax_ps=max(xpcorn_ps(1),xpcorn_ps(2),xpcorn_ps(3),xpcorn_ps(4),
     &  xpcorn_ps(5),xpcorn_ps(6),xpcorn_ps(7),xpcorn_ps(8))

      wymin_ps=min(ypcorn_ps(1),ypcorn_ps(2),ypcorn_ps(3),ypcorn_ps(4),
     &  ypcorn_ps(5),ypcorn_ps(6),ypcorn_ps(7),ypcorn_ps(8))
      wymax_ps=max(ypcorn_ps(1),ypcorn_ps(2),ypcorn_ps(3),ypcorn_ps(4),
     &  ypcorn_ps(5),ypcorn_ps(6),ypcorn_ps(7),ypcorn_ps(8))

      scalex_ps=xsiz_ps/(wxmax_ps-wxmin_ps)
      scaley_ps=ysiz_ps/(wymax_ps-wymin_ps)

      if (log10x_ps.eq.0) then
        xmin3d_ps=xmin
        xmax3d_ps=xmax
      else
        if (xmin.le.0.0) xmin=1.0e-30
        if (xmax.le.xmin) xmax=1.0e30
        xmin3d_ps=alog10(xmin)
        xmax3d_ps=alog10(xmax)
      endif

      if (log10y_ps.eq.0) then
        ymin3d_ps=ymin
        ymax3d_ps=ymax
      else
        if (ymin.le.0.0) ymin=1.0e-30
        if (ymax.le.ymin) ymax=1.0e30
        ymin3d_ps=alog10(ymin)
        ymax3d_ps=alog10(ymax)
      endif

      if (log10z_ps.eq.0) then
        zmin3d_ps=zmin
        zmax3d_ps=zmax
      else
        if (zmin.le.0.0) zmin=1.0e-30
        if (zmax.le.zmin) zmax=1.0e30
        zmin3d_ps=alog10(zmin)
        zmax3d_ps=alog10(zmax)
      endif

      call mshplt_set_line_style(1)

      if (iaxis.eq.1) then

        if (theta_ps.le.90.0) then
          if (phi_ps.le.90.0) then

            if (log10x_ps.eq.0) then
              call mshplt_view_axis(
     &          xpcorn_ps(1),xpcorn_ps(2),ypcorn_ps(1),ypcorn_ps(2),
     &          xmin,xmax,xtitd(1:len_trim(xtitd)),-1,0.,0.,
     &          xtitoff_ps,
     &          xlaboff_ps) !x-axis
            else
              xyplen=sqrt(
     &          (xpcorn_ps(2)-xpcorn_ps(1))**2+
     &          (ypcorn_ps(2)-ypcorn_ps(1))**2
     &          )
              ang=atan2(
     &          (ypcorn_ps(2)-ypcorn_ps(1)),
     &          (xpcorn_ps(2)-xpcorn_ps(1))
     &          )*radtodeg_ps
              iticside=-1
              ilabside=1
              titang=0.
              rellabang=-ang
              call mshplt_view_log_axis(
     &          xpcorn_ps(1),ypcorn_ps(1),xyplen,
     &          xminin,xmaxin,ang,
     &          chhe,ticheight,xtit,
     &          iticside,ilabside,
     &          titang,ytitoff_ps,
     &          rellabang,xlaboff_ps) !x-axis
            endif

            if (log10y_ps.eq.0) then
              call mshplt_view_axis(
     &          xpcorn_ps(1),xpcorn_ps(4),ypcorn_ps(1),ypcorn_ps(4),
     &          ymin,ymax,ytitd(1:len_trim(ytitd)),1,180.,180.,
     &          ytitoff_ps,
     &          -ylaboff_ps) !y-axis
            else
              xyplen=sqrt(
     &          (xpcorn_ps(4)-xpcorn_ps(1))**2+
     &          (ypcorn_ps(4)-ypcorn_ps(1))**2
     &          )
              ang=atan2(
     &          (ypcorn_ps(4)-ypcorn_ps(1)),
     &          (xpcorn_ps(4)-xpcorn_ps(1))
     &          )*radtodeg_ps
              iticside=1
              ilabside=-1
              titang=180.
              rellabang=-ang
              call mshplt_view_log_axis(
     &          xpcorn_ps(1),ypcorn_ps(1),xyplen,
     &          yminin,ymaxin,ang,
     &          chhe,ticheight,ytit,
     &          iticside,ilabside,
     &          titang,ztitoff_ps,
     &          rellabang,-2.*zlaboff_ps) !y-axis
            endif

            if (log10z_ps.eq.0) then
              call mshplt_view_axis(
     &          xpcorn_ps(4),xpcorn_ps(8),ypcorn_ps(4),ypcorn_ps(8),
     &          zmin,zmax,ztitd(1:len_trim(ztitd)),1,-90.,0.,
     &          -ztitoff_ps,-zlaboff_ps) !z-axis
            else
              xyplen=sqrt(
     &          (xpcorn_ps(8)-xpcorn_ps(4))**2+
     &          (ypcorn_ps(8)-ypcorn_ps(4))**2
     &          )
              ang=atan2(
     &          (ypcorn_ps(8)-ypcorn_ps(4)),
     &          (xpcorn_ps(8)-xpcorn_ps(4))
     &          )*radtodeg_ps
              iticside=1
              ilabside=-1
              titang=0.
              rellabang=-90.
              call mshplt_view_log_axis(
     &          xpcorn_ps(4),ypcorn_ps(4),xyplen,
     &          zminin,zmaxin,ang,
     &          chhe,ticheight,ztit,
     &          iticside,ilabside,
     &          titang,ztitoff_ps,
     &          rellabang,-zlaboff_ps*1.75) !z-axis
            endif

          else if (phi_ps.le.180.0) then

            if (log10x_ps.eq.0) then
              call mshplt_view_axis(
     &          xpcorn_ps(4),xpcorn_ps(3),ypcorn_ps(4),ypcorn_ps(3),
     &          xmin,xmax,xtitd(1:len_trim(xtitd)),1,180.,180.,
     &          xtitoff_ps,
     &          -xlaboff_ps) !x-axis
            else
              xyplen=sqrt(
     &          (xpcorn_ps(3)-xpcorn_ps(4))**2+
     &          (ypcorn_ps(3)-ypcorn_ps(4))**2
     &          )
              ang=atan2(
     &          (ypcorn_ps(3)-ypcorn_ps(4)),
     &          (xpcorn_ps(3)-xpcorn_ps(4))
     &          )*radtodeg_ps
              iticside=1
              ilabside=-1
              titang=180.
              rellabang=-ang
              call mshplt_view_log_axis(
     &          xpcorn_ps(4),ypcorn_ps(4),xyplen,
     &          xminin,xmaxin,ang,
     &          chhe,ticheight,xtit,
     &          iticside,ilabside,
     &          titang,ztitoff_ps,
     &          rellabang,-2.*xlaboff_ps) !x-axis
            endif

            if (log10y_ps.eq.0) then
              call mshplt_view_axis(
     &          xpcorn_ps(1),xpcorn_ps(4),ypcorn_ps(1),ypcorn_ps(4),
     &          ymin,ymax,ytitd(1:len_trim(ytitd)),1,180.,180.,
     &          ytitoff_ps,
     &          -ylaboff_ps) !y-axis
            else
              xyplen=sqrt(
     &          (xpcorn_ps(4)-xpcorn_ps(1))**2+
     &          (ypcorn_ps(4)-ypcorn_ps(1))**2
     &          )
              ang=atan2(
     &          (ypcorn_ps(4)-ypcorn_ps(1)),
     &          (xpcorn_ps(4)-xpcorn_ps(1))
     &          )*radtodeg_ps
              iticside=1
              ilabside=1
              titang=180.
              rellabang=-ang
              call mshplt_view_log_axis(
     &          xpcorn_ps(1),ypcorn_ps(1),xyplen,
     &          yminin,ymaxin,ang,
     &          chhe,ticheight,ytit,
     &          iticside,ilabside,
     &          titang,ztitoff_ps,
     &          rellabang,-1.5*zlaboff_ps) !y-axis
            endif

            if (log10z_ps.eq.0) then
              call mshplt_view_axis(
     &          xpcorn_ps(3),xpcorn_ps(7),ypcorn_ps(3),ypcorn_ps(7),
     &          zmin,zmax,ztitd(1:len_trim(ztitd)),1,-90.,0.,
     &          -ztitoff_ps,-zlaboff_ps) !z-axis
            else
              xyplen=sqrt(
     &          (xpcorn_ps(7)-xpcorn_ps(3))**2+
     &          (ypcorn_ps(7)-ypcorn_ps(3))**2
     &          )
              ang=atan2(
     &          (ypcorn_ps(7)-ypcorn_ps(3)),
     &          (xpcorn_ps(7)-xpcorn_ps(3))
     &          )*radtodeg_ps
              iticside=1
              ilabside=-1
              titang=0.
              rellabang=-90.
              call mshplt_view_log_axis(
     &          xpcorn_ps(3),ypcorn_ps(3),xyplen,
     &          zminin,zmaxin,ang,
     &          chhe,ticheight,ztit,
     &          iticside,ilabside,
     &          titang,ztitoff_ps,
     &          rellabang,zlaboff_ps*0.6) !z-axis
            endif

          else if (phi_ps.le.270.0) then

            if (log10x_ps.eq.0) then
              call mshplt_view_axis(
     &          xpcorn_ps(4),xpcorn_ps(3),ypcorn_ps(4),ypcorn_ps(3),
     &          xmin,xmax,xtitd(1:len_trim(xtitd)),1,180.,180.,
     &          xtitoff_ps,
     &          -xlaboff_ps) !x-axis
            else
              xyplen=sqrt(
     &          (xpcorn_ps(3)-xpcorn_ps(4))**2+
     &          (ypcorn_ps(3)-ypcorn_ps(4))**2
     &          )
              ang=atan2(
     &          (ypcorn_ps(3)-ypcorn_ps(4)),
     &          (xpcorn_ps(3)-xpcorn_ps(4))
     &          )*radtodeg_ps
              iticside=1
              ilabside=-1
              titang=180.
              rellabang=-ang
              call mshplt_view_log_axis(
     &          xpcorn_ps(4),ypcorn_ps(4),xyplen,
     &          xminin,xmaxin,ang,
     &          chhe,ticheight,xtit,
     &          iticside,ilabside,
     &          titang,ztitoff_ps*0.7,
     &          rellabang,-2.*xlaboff_ps) !x-axis
            endif

            if (log10y_ps.eq.0) then
              call mshplt_view_axis(
     &          xpcorn_ps(2),xpcorn_ps(3),ypcorn_ps(2),ypcorn_ps(3),
     &          ymin,ymax,ytitd(1:len_trim(ytitd)),-1,0.,0.,
     &          1.1*ytitoff_ps,
     &          1.25*ylaboff_ps) !y-axis
            else
              xyplen=sqrt(
     &          (xpcorn_ps(3)-xpcorn_ps(2))**2+
     &          (ypcorn_ps(3)-ypcorn_ps(2))**2
     &          )
              ang=atan2(
     &          (ypcorn_ps(3)-ypcorn_ps(2)),
     &          (xpcorn_ps(3)-xpcorn_ps(2))
     &          )*radtodeg_ps
              iticside=-1
              ilabside=-1
              titang=0.
              rellabang=-ang
              call mshplt_view_log_axis(
     &          xpcorn_ps(2),ypcorn_ps(2),xyplen,
     &          yminin,ymaxin,ang,
     &          chhe,ticheight,ytit,
     &          iticside,ilabside,
     &          titang,ztitoff_ps,
     &          rellabang,1.25*zlaboff_ps) !y-axis
            endif

            if (log10z_ps.eq.0) then
              call mshplt_view_axis(
     &          xpcorn_ps(2),xpcorn_ps(6),ypcorn_ps(2),ypcorn_ps(6),
     &          zmin,zmax,ztitd(1:len_trim(ztitd)),1,-90.,0.,
     &          -ztitoff_ps,-zlaboff_ps) !z-axis
            else
              xyplen=sqrt(
     &          (xpcorn_ps(6)-xpcorn_ps(2))**2+
     &          (ypcorn_ps(6)-ypcorn_ps(2))**2
     &          )
              ang=atan2(
     &          (ypcorn_ps(6)-ypcorn_ps(2)),
     &          (xpcorn_ps(6)-xpcorn_ps(2))
     &          )*radtodeg_ps
              iticside=1
              ilabside=-1
              titang=0.
              rellabang=-90.
              call mshplt_view_log_axis(
     &          xpcorn_ps(2),ypcorn_ps(2),xyplen,
     &          zminin,zmaxin,ang,
     &          chhe,ticheight,ztit,
     &          iticside,ilabside,
     &          titang,ztitoff_ps,
     &          rellabang,-3.*zlaboff_ps*0.6) !z-axis
            endif

          else if (phi_ps.le.360.0) then

            if (log10x_ps.eq.0) then
              call mshplt_view_axis(
     &          xpcorn_ps(1),xpcorn_ps(2),ypcorn_ps(1),ypcorn_ps(2),
     &          xmin,xmax,xtitd(1:len_trim(xtitd)),-1,0.,0.,
     &          1.1*xtitoff_ps,
     &          1.2*xlaboff_ps) !x-axis
            else
              xyplen=sqrt(
     &          (xpcorn_ps(1)-xpcorn_ps(2))**2+
     &          (ypcorn_ps(1)-ypcorn_ps(2))**2
     &          )
              ang=atan2(
     &          (ypcorn_ps(2)-ypcorn_ps(1)),
     &          (xpcorn_ps(2)-xpcorn_ps(1))
     &          )*radtodeg_ps
              iticside=-1
              ilabside=-1
              titang=0.
              rellabang=-ang
              call mshplt_view_log_axis(
     &          xpcorn_ps(1),ypcorn_ps(1),xyplen,
     &          xminin,xmaxin,ang,
     &          chhe,ticheight,xtit,
     &          iticside,ilabside,
     &          titang,ztitoff_ps*0.7,
     &          rellabang,xlaboff_ps) !x-axis
            endif

            if (log10y_ps.eq.0) then
              call mshplt_view_axis(
     &          xpcorn_ps(2),xpcorn_ps(3),ypcorn_ps(2),ypcorn_ps(3),
     &          ymin,ymax,ytitd(1:len_trim(ytitd)),-1,0.,0.,
     &          1.2*ytitoff_ps,
     &          1.1*ylaboff_ps) !y-axis
            else
              xyplen=sqrt(
     &          (xpcorn_ps(3)-xpcorn_ps(2))**2+
     &          (ypcorn_ps(3)-ypcorn_ps(2))**2
     &          )
              ang=atan2(
     &          (ypcorn_ps(3)-ypcorn_ps(2)),
     &          (xpcorn_ps(3)-xpcorn_ps(2))
     &          )*radtodeg_ps
              iticside=-1
              ilabside=1
              titang=0.
              rellabang=-ang
              call mshplt_view_log_axis(
     &          xpcorn_ps(2),ypcorn_ps(2),xyplen,
     &          yminin,ymaxin,ang,
     &          chhe,ticheight,ytit,
     &          iticside,ilabside,
     &          titang,ztitoff_ps*0.8,
     &          rellabang,zlaboff_ps) !y-axis
            endif

            if (log10z_ps.eq.0) then
              call mshplt_view_axis(
     &          xpcorn_ps(1),xpcorn_ps(5),ypcorn_ps(1),ypcorn_ps(5),
     &          zmin,zmax,ztitd(1:len_trim(ztitd)),1,-90.,0.,
     &          -ztitoff_ps,-zlaboff_ps) !z-axis
            else
              xyplen=sqrt(
     &          (xpcorn_ps(5)-xpcorn_ps(1))**2+
     &          (ypcorn_ps(5)-ypcorn_ps(1))**2
     &          )
              ang=atan2(
     &          (ypcorn_ps(5)-ypcorn_ps(1)),
     &          (xpcorn_ps(5)-xpcorn_ps(1))
     &          )*radtodeg_ps
              iticside=1
              ilabside=-1
              titang=0.
              rellabang=-90.
              call mshplt_view_log_axis(
     &          xpcorn_ps(1),ypcorn_ps(1),xyplen,
     &          zminin,zmaxin,ang,
     &          chhe,ticheight,ztit,
     &          iticside,ilabside,
     &          titang,ztitoff_ps,
     &          rellabang,-3.*zlaboff_ps*0.6) !z-axis
            endif

          endif !if (phi_ps.le.90.0) then

        else if (theta_ps.le.180.0) then

          if (phi_ps.le.90.0) then

            if (log10x_ps.eq.0) then
              call mshplt_view_axis(
     &          xpcorn_ps(8),xpcorn_ps(7),ypcorn_ps(8),ypcorn_ps(7),
     &          xmin,xmax,xtitd(1:len_trim(xtitd)),-1,0.,0.,
     &          xtitoff_ps,
     &          xlaboff_ps) !x-axis
            else
              xyplen=sqrt(
     &          (xpcorn_ps(7)-xpcorn_ps(8))**2+
     &          (ypcorn_ps(7)-ypcorn_ps(8))**2
     &          )
              ang=atan2(
     &          (ypcorn_ps(7)-ypcorn_ps(8)),
     &          (xpcorn_ps(7)-xpcorn_ps(8))
     &          )*radtodeg_ps
              iticside=-1
              ilabside=1
              titang=0.
              rellabang=-ang
              call mshplt_view_log_axis(
     &          xpcorn_ps(8),ypcorn_ps(8),xyplen,
     &          xminin,xmaxin,ang,
     &          chhe,ticheight,xtit,
     &          iticside,ilabside,
     &          titang,ytitoff_ps,
     &          rellabang,xlaboff_ps) !x-axis
            endif

            if (log10y_ps.eq.0) then
              call mshplt_view_axis(
     &          xpcorn_ps(6),xpcorn_ps(7),ypcorn_ps(6),ypcorn_ps(7),
     &          ymin,ymax,ytitd(1:len_trim(ytitd)),1,180.,180.,
     &          ytitoff_ps,
     &          -ylaboff_ps) !y-axis
            else
              xyplen=sqrt(
     &          (xpcorn_ps(7)-xpcorn_ps(6))**2+
     &          (ypcorn_ps(7)-ypcorn_ps(6))**2
     &          )
              ang=atan2(
     &          (ypcorn_ps(7)-ypcorn_ps(6)),
     &          (xpcorn_ps(7)-xpcorn_ps(6))
     &          )*radtodeg_ps
              iticside=1
              ilabside=1
              titang=180.
              rellabang=-ang
              call mshplt_view_log_axis(
     &          xpcorn_ps(6),ypcorn_ps(6),xyplen,
     &          yminin,ymaxin,ang,
     &          chhe,ticheight,ytit,
     &          iticside,ilabside,
     &          titang,ztitoff_ps,
     &          rellabang,-1.5*zlaboff_ps) !y-axis
            endif

            if (log10z_ps.eq.0) then
              call mshplt_view_axis(
     &          xpcorn_ps(4),xpcorn_ps(8),ypcorn_ps(4),ypcorn_ps(8),
     &          zmin,zmax,ztitd(1:len_trim(ztitd)),1,+90.,0.,
     &          ztitoff_ps,zlaboff_ps) !z-axis
            else
              xyplen=sqrt(
     &          (xpcorn_ps(8)-xpcorn_ps(4))**2+
     &          (ypcorn_ps(8)-ypcorn_ps(4))**2
     &          )
              ang=atan2(
     &          (ypcorn_ps(8)-ypcorn_ps(4)),
     &          (xpcorn_ps(8)-xpcorn_ps(4))
     &          )*radtodeg_ps
              iticside=-1
              ilabside=-1
              titang=0.
              rellabang=-ang
              call mshplt_view_log_axis(
     &          xpcorn_ps(4),ypcorn_ps(4),xyplen,
     &          zminin,zmaxin,ang,
     &          chhe,ticheight,ztit,
     &          iticside,ilabside,
     &          titang,ztitoff_ps,
     &          rellabang,-1.75*zlaboff_ps) !z-axis
            endif

          else if (phi_ps.le.180.0) then

            if (log10x_ps.eq.0) then
              call mshplt_view_axis(
     &          xpcorn_ps(5),xpcorn_ps(6),ypcorn_ps(5),ypcorn_ps(6),
     &          xmin,xmax,xtitd(1:len_trim(xtitd)),1,180.,180.,
     &          1.1*xtitoff_ps,
     &          -xlaboff_ps) !x-axis
            else
              xyplen=sqrt(
     &          (xpcorn_ps(6)-xpcorn_ps(5))**2+
     &          (ypcorn_ps(6)-ypcorn_ps(5))**2
     &          )
              ang=atan2(
     &          (ypcorn_ps(6)-ypcorn_ps(5)),
     &          (xpcorn_ps(6)-xpcorn_ps(5))
     &          )*radtodeg_ps
              iticside=1
              ilabside=1
              titang=180.
              rellabang=-ang
              call mshplt_view_log_axis(
     &          xpcorn_ps(5),ypcorn_ps(5),xyplen,
     &          xminin,xmaxin,ang,
     &          chhe,ticheight,xtit,
     &          iticside,ilabside,
     &          titang,ztitoff_ps,
     &          rellabang,-1.5*xlaboff_ps) !x-axis
            endif

            if (log10y_ps.eq.0) then
              call mshplt_view_axis(
     &          xpcorn_ps(6),xpcorn_ps(7),ypcorn_ps(6),ypcorn_ps(7),
     &          ymin,ymax,ytitd(1:len_trim(ytitd)),1,180.,180.,
     &          ytitoff_ps,
     &          -ylaboff_ps) !y-axis
            else
              xyplen=sqrt(
     &          (xpcorn_ps(7)-xpcorn_ps(6))**2+
     &          (ypcorn_ps(7)-ypcorn_ps(6))**2
     &          )
              ang=atan2(
     &          (ypcorn_ps(7)-ypcorn_ps(6)),
     &          (xpcorn_ps(7)-xpcorn_ps(6))
     &          )*radtodeg_ps
              iticside=1
              ilabside=1
              titang=180.
              rellabang=-ang
              call mshplt_view_log_axis(
     &          xpcorn_ps(6),ypcorn_ps(6),xyplen,
     &          yminin,ymaxin,ang,
     &          chhe,ticheight,ytit,
     &          iticside,ilabside,
     &          titang,ztitoff_ps,
     &          rellabang,-2.*zlaboff_ps) !y-axis
            endif

            if (log10z_ps.eq.0) then
              call mshplt_view_axis(
     &          xpcorn_ps(3),xpcorn_ps(7),ypcorn_ps(3),ypcorn_ps(7),
     &          zmin,zmax,ztitd(1:len_trim(ztitd)),-1,+90.,180.,
     &          -ztitoff_ps,zlaboff_ps) !z-axis
            else
              xyplen=sqrt(
     &          (xpcorn_ps(7)-xpcorn_ps(3))**2+
     &          (ypcorn_ps(7)-ypcorn_ps(3))**2
     &          )
              ang=atan2(
     &          (ypcorn_ps(7)-ypcorn_ps(3)),
     &          (xpcorn_ps(7)-xpcorn_ps(3))
     &          )*radtodeg_ps
              iticside=-1
              ilabside=-1
              titang=0.
              rellabang=+90.
              call mshplt_view_log_axis(
     &          xpcorn_ps(3),ypcorn_ps(3),xyplen,
     &          zminin,zmaxin,ang,
     &          chhe,ticheight,ztit,
     &          iticside,ilabside,
     &          titang,ztitoff_ps,
     &          rellabang,-1.75*zlaboff_ps) !z-axis
            endif

          else if (phi_ps.le.270.0) then

            if (log10x_ps.eq.0) then
              call mshplt_view_axis(
     &          xpcorn_ps(5),xpcorn_ps(6),ypcorn_ps(5),ypcorn_ps(6),
     &          xmin,xmax,xtitd(1:len_trim(xtitd)),1,180.,180.,
     &          1.2*xtitoff_ps,
     &          -1.25*xlaboff_ps) !x-axis
            else
              xyplen=sqrt(
     &          (xpcorn_ps(6)-xpcorn_ps(5))**2+
     &          (ypcorn_ps(6)-ypcorn_ps(5))**2
     &          )
              ang=atan2(
     &          (ypcorn_ps(6)-ypcorn_ps(5)),
     &          (xpcorn_ps(6)-xpcorn_ps(5))
     &          )*radtodeg_ps
              iticside=1
              ilabside=1
              titang=180.
              rellabang=-ang
              call mshplt_view_log_axis(
     &          xpcorn_ps(5),ypcorn_ps(5),xyplen,
     &          xminin,xmaxin,ang,
     &          chhe,ticheight,xtit,
     &          iticside,ilabside,
     &          titang,ztitoff_ps*0.7,
     &          rellabang,-2.*xlaboff_ps) !x-axis
            endif

            if (log10y_ps.eq.0) then
              call mshplt_view_axis(
     &          xpcorn_ps(5),xpcorn_ps(8),ypcorn_ps(5),ypcorn_ps(8),
     &          ymin,ymax,ytitd(1:len_trim(ytitd)),-1,0.,0.,
     &          1.1*ytitoff_ps,
     &          1.25*ylaboff_ps) !y-axis
            else
              xyplen=sqrt(
     &          (xpcorn_ps(8)-xpcorn_ps(5))**2+
     &          (ypcorn_ps(8)-ypcorn_ps(5))**2
     &          )
              ang=atan2(
     &          (ypcorn_ps(8)-ypcorn_ps(5)),
     &          (xpcorn_ps(8)-xpcorn_ps(5))
     &          )*radtodeg_ps
              iticside=-1
              ilabside=1
              titang=0.
              rellabang=-ang
              call mshplt_view_log_axis(
     &          xpcorn_ps(5),ypcorn_ps(5),xyplen,
     &          yminin,ymaxin,ang,
     &          chhe,ticheight,ytit,
     &          iticside,ilabside,
     &          titang,ztitoff_ps,
     &          rellabang,zlaboff_ps*0.5) !y-axis
            endif

            if (log10z_ps.eq.0) then
              call mshplt_view_axis(
     &          xpcorn_ps(2),xpcorn_ps(6),ypcorn_ps(2),ypcorn_ps(6),
     &          zmin,zmax,ztitd(1:len_trim(ztitd)),-1,+90.,0.,
     &          ztitoff_ps,zlaboff_ps) !z-axis
            else
              xyplen=sqrt(
     &          (xpcorn_ps(6)-xpcorn_ps(2))**2+
     &          (ypcorn_ps(6)-ypcorn_ps(2))**2
     &          )
              ang=atan2(
     &          (ypcorn_ps(6)-ypcorn_ps(2)),
     &          (xpcorn_ps(6)-xpcorn_ps(2))
     &          )*radtodeg_ps
              iticside=-1
              ilabside=-1
              titang=0.
              rellabang=90.
              call mshplt_view_log_axis(
     &          xpcorn_ps(2),ypcorn_ps(2),xyplen,
     &          zminin,zmaxin,ang,
     &          chhe,ticheight,ztit,
     &          iticside,ilabside,
     &          titang,ztitoff_ps,
     &          rellabang,-1.75*zlaboff_ps) !z-axis
            endif

          else if (phi_ps.le.360.0) then

            if (log10x_ps.eq.0) then
              call mshplt_view_axis(
     &          xpcorn_ps(8),xpcorn_ps(7),ypcorn_ps(8),ypcorn_ps(7),
     &          xmin,xmax,xtitd(1:len_trim(xtitd)),-1,0.,0.,
     &          1.1*xtitoff_ps,
     &          1.2*xlaboff_ps) !x-axis
            else
              xyplen=sqrt(
     &          (xpcorn_ps(7)-xpcorn_ps(8))**2+
     &          (ypcorn_ps(7)-ypcorn_ps(8))**2
     &          )
              ang=atan2(
     &          (ypcorn_ps(7)-ypcorn_ps(8)),
     &          (xpcorn_ps(7)-xpcorn_ps(8))
     &          )*radtodeg_ps
              iticside=-1
              ilabside=1
              titang=0.
              rellabang=-ang
              call mshplt_view_log_axis(
     &          xpcorn_ps(8),ypcorn_ps(8),xyplen,
     &          xminin,xmaxin,ang,
     &          chhe,ticheight,xtit,
     &          iticside,ilabside,
     &          titang,ztitoff_ps*0.7,
     &          rellabang,xlaboff_ps*0.5) !x-axis
            endif

            if (log10y_ps.eq.0) then
              call mshplt_view_axis(
     &          xpcorn_ps(5),xpcorn_ps(8),ypcorn_ps(5),ypcorn_ps(8),
     &          ymin,ymax,ytitd(1:len_trim(ytitd)),-1,0.,0.,
     &          1.2*ytitoff_ps,
     &          1.1*ylaboff_ps) !y-axis
            else
              xyplen=sqrt(
     &          (xpcorn_ps(8)-xpcorn_ps(5))**2+
     &          (ypcorn_ps(8)-ypcorn_ps(5))**2
     &          )
              ang=atan2(
     &          (ypcorn_ps(8)-ypcorn_ps(5)),
     &          (xpcorn_ps(8)-xpcorn_ps(5))
     &          )*radtodeg_ps
              iticside=-1
              ilabside=1
              titang=0.
              rellabang=-ang
              call mshplt_view_log_axis(
     &          xpcorn_ps(5),ypcorn_ps(5),xyplen,
     &          yminin,ymaxin,ang,
     &          chhe,ticheight,ytit,
     &          iticside,ilabside,
     &          titang,ztitoff_ps,
     &          rellabang,zlaboff_ps) !y-axis
            endif

            if (log10z_ps.eq.0) then
              call mshplt_view_axis(
     &          xpcorn_ps(1),xpcorn_ps(5),ypcorn_ps(1),ypcorn_ps(5),
     &          zmin,zmax,ztitd(1:len_trim(ztitd)),-1,90.,0.,
     &          ztitoff_ps,zlaboff_ps) !z-axis
            else
              xyplen=sqrt(
     &          (xpcorn_ps(5)-xpcorn_ps(1))**2+
     &          (ypcorn_ps(5)-ypcorn_ps(1))**2
     &          )
              ang=atan2(
     &          (ypcorn_ps(5)-ypcorn_ps(1)),
     &          (xpcorn_ps(5)-xpcorn_ps(1))
     &          )*radtodeg_ps
              iticside=-1
              ilabside=-1
              titang=0.
              rellabang=90.
              call mshplt_view_log_axis(
     &          xpcorn_ps(1),ypcorn_ps(1),xyplen,
     &          zminin,zmaxin,ang,
     &          chhe,ticheight,ztit,
     &          iticside,ilabside,
     &          titang,ztitoff_ps,
     &          rellabang,-1.75*zlaboff_ps) !z-axis
            endif

          endif !if (phi_ps.le.90.0) then
        endif !if (theta_ps.le.90.0) then

      endif !iaxis.eq.1

      if (ibox.ne.0) then
        call mshplt_line_raw(xpcorn_ps(1),ypcorn_ps(1),xpcorn_ps(2),ypcorn_ps(2))
        call mshplt_line_raw(xpcorn_ps(1),ypcorn_ps(1),xpcorn_ps(4),ypcorn_ps(4))
        call mshplt_line_raw(xpcorn_ps(4),ypcorn_ps(4),xpcorn_ps(8),ypcorn_ps(8))
        call mshplt_line_raw(xpcorn_ps(4),ypcorn_ps(4),xpcorn_ps(3),ypcorn_ps(3))
        call mshplt_line_raw(xpcorn_ps(2),ypcorn_ps(2),xpcorn_ps(3),ypcorn_ps(3))
        call mshplt_line_raw(xpcorn_ps(1),ypcorn_ps(1),xpcorn_ps(5),ypcorn_ps(5))
        call mshplt_line_raw(xpcorn_ps(2),ypcorn_ps(2),xpcorn_ps(6),ypcorn_ps(6))
        call mshplt_line_raw(xpcorn_ps(3),ypcorn_ps(3),xpcorn_ps(7),ypcorn_ps(7))
        call mshplt_pline_raw_closed(4,xpcorn_ps(5),ypcorn_ps(5))
      endif

      call mshplt_get_line_color(ilinecoloro,kredo,kblueo,kgreeno)
      call mshplt_set_line_style(ilinestyleo)

      if (kzone_ps.lt.0) then
        kzone_ps=ifirst_ps
      endif

      write(lun_ps,'(a)')'% end of mshplt_frame3d'

      return
      end
*CMZ :  0.01/03 23/09/2014  09.38.57  by  Michael Scheer
*CMZ :  0.01/00 24/08/2014  13.17.26  by  Michael Scheer
*CMZ :  0.00/06 20/08/2014  12.49.55  by  Michael Scheer
*CMZ :  0.00/04 14/08/2014  12.58.41  by  Michael Scheer
*CMZ :  0.00/02 09/07/2014  16.17.24  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_text(x,y,text)

      implicit none

*KEEP,MSHPLT.
      real
     &  pttomm_ps,pttocm_ps, !convert from pt to mm or cm respectively
     &  scale_ps, isscale_ps, !current scale to convert from pt
     &  scalex_ps, !current scale to convert x to pad system
     &  scaley_ps, !current scale to convert y to pad system
     &  xold_ps,yold_ps, !old pen position on pad
     &  x_ps,y_ps, !current pen position on pad
     &  rlinewidth_ps !line width on pad
     & ,chhe_ps, chhe_index_ps ! character height
     & ,ticsiz_ps ! tic height
     & ,ang_ps,tang_ps ! current angle and text angle
     &  ,offgtitx_ps !offset of global title on pad
     &  ,offgtity_ps !offset of global title on pad
     &  ,offdatex_ps !offset of date on pad
     &  ,offdatey_ps !offset of date on pad
     &  ,aLabHeight_ps !character height of axis label
     &  ,xoffexp_ps !xoffset of power term
     &  ,yoffexp_ps !yoffset of power term
     &  ,xlaboff_ps !offset of x-axis label on pad
     &  ,xtitoff_ps !offset of x-axis title on pad
     &  ,ylaboff_ps ! offset of y-axis label on pad
     &  ,ytitoff_ps !offset of y-axis title on pad
     &  ,ztitoff_ps !offset of z-axis title on pad
     &  ,zlaboff_ps ! offset of z-axis label on pad
     &  ,xsiz_ps,ysiz_ps !size on pad of plotting canvas
     &  ,wxmin_ps,wxmax_ps,wymin_ps,wymax_ps !world
     &  ,xleft_ps,xright_ps,ybottom_ps,ytop_ps ! pad on pad
     &  ,xsizorig_ps,ysizorig_ps !size on pad of plotting canvas for paper scale
     &  ,xleftorig_ps,xrightorig_ps,ybottomorig_ps,ytoporig_ps ! orig pad on paper
     &  ,rmtyp_ps(2,0:100) !postscript encoding of characters
     &  ,rmsiz_ps ! marker size
     &  ,phi_ps,theta_ps ! 3d viewing angles
     &  ,pi_ps !Pi
     &  ,xmin3d_ps,xmax3d_ps
     &  ,ymin3d_ps,ymax3d_ps
     &  ,zmin3d_ps,zmax3d_ps
     &  ,xmgl_ps,xmgr_ps,ymgu_ps,ymgl_ps !old hplot variables for margins
     &  ,xwin_ps,ywin_ps ! horizontal and vertical distance of pads
     &  ,ygti_ps !y-offset of global title for hplot
     &  ,gsiz_ps !character size of global title for hplot
     &  ,xcen_ps,ycen_ps !center on canvas
     &  ,scaletxt_ps ! Scaling for text entities
     &  ,ihigzmode_ps !HIGZ-compatibility mode
     &  ,tsiz_ps ! axis-title-size in HIGZ-mode
     &  ,coln_ps,cnLine_ps,cnFill_ps,cnMarker_ps,cnText_ps,cnFrame_ps !color norms
     &  ,clips_ps(2,1000) !clipping range
     &  ,radtodeg_ps
     &  ,xpcorn_ps(8),ypcorn_ps(8) !3d frame box
     &  ,ifbox_ps,xfb_ps(9),yfb_ps(9)

      parameter(pttocm_ps=72/2.54,pttomm_ps=72/25.4)

      integer nmaxfile_ps
      parameter(nmaxfile_ps=10000)

      integer
     &  lun_ps,istat_ps
     &  ,ilabmod_ps !mode off plotting label
     &  ,inolabs_ps !no axis labeling
     &  ,irunviewer_ps !start live viewer
     &  ,iviewinter_ps !interactive stopping of live viewer
     &  ,mtyp_ps ! current marker type
     &  ,icolor_ps, ired_ps,igreen_ps,iblue_ps !color, if <1, rgb values are used
     &  ,ilinestyle_ps !line style: 1: solid, 2: dashed, 3: dotted, 4: dashed-dotted
     &  ,kzone_ps, itouched_ps !current zone and flag of usage
     &  ,nzone_ps,nxzone_ps,nyzone_ps !zones of canvas
     &  ,ifirst_ps !first zone for plotting
     &  ,isamecanvas_ps,isameframe_ps
     &  ,nfile_ps !file counter
     &  ,lunall_ps(nmaxfile_ps) !LUNs of plotfiles
     &  ,kbbxl_ps,kbbxr_ps,kbbyb_ps,kbbyt_ps !bounding box (kbbxl_ps,kbbxr_ps,...)
     &  ,inewpage_ps !new page
     &  ,kFrameColor_ps,kFrameRed_ps,kFrameBlue_ps,kFrameGreen_ps
     &  ,kFillColor_ps,kFillRed_ps,kFillBlue_ps,kFillGreen_ps
     &  ,kLineColor_ps,kLineRed_ps,kLineBlue_ps,kLineGreen_ps
     &  ,kTextColor_ps,kTextRed_ps,kTextBlue_ps,kTextGreen_ps
     &  ,kMarkerColor_ps,kMarkerRed_ps,kMarkerBlue_ps,kMarkerGreen_ps
     &  ,idrawgtit_ps !draw global title on top of page
     &  ,lunbase_ps
     &  ,log10x_ps,log10y_ps,log10z_ps !flags for logarithmic axis
     &  ,nclips_ps ! number of clipping points
     &  ,mode3d_ps
     &  ,isbox_ps,isdate_ps,isoffdate_ps
     &  ,kBox_ps,kDate_ps ! to plot surrounding box and date on top of plots

      character(2048)
     &  file_ps,viewer_ps,cline_ps,filebase_ps,
     &  xtit_ps,ytit_ps,ztit_ps,gtit_ps !axis-title and global title
     &  ,fileall_ps(nmaxfile_ps) !plotfiles

      character(4) chmarker_ps !current marker

      integer nbuffsize_ps,ibuffpos_ps
      parameter(nbuffsize_ps=100000)
      character(2048) chbuff_ps(nbuffsize_ps)
      character(4) chch_ps(1000)

      common/mshpltc/
     &  scale_ps, isscale_ps,scalex_ps,scaley_ps,
     &  xold_ps,yold_ps,x_ps,y_ps,rlinewidth_ps,chhe_ps, chhe_index_ps,
     &  ticsiz_ps,
     &  ang_ps,tang_ps
     & ,xlaboff_ps,ylaboff_ps
     & ,xtitoff_ps,ytitoff_ps
     &  ,ztitoff_ps !offset of z-axis title on pad
     &  ,zlaboff_ps ! offset of z-axis label on pad
     & ,xsiz_ps,ysiz_ps
     &  ,wxmin_ps,wxmax_ps,wymin_ps,wymax_ps !world
     &  ,xleft_ps,xright_ps,ybottom_ps,ytop_ps ! pad on pad
     &  ,ilabmod_ps,rmtyp_ps,rmsiz_ps
     & ,lun_ps,istat_ps,inolabs_ps,irunviewer_ps,iviewinter_ps,mtyp_ps
     &  ,file_ps,viewer_ps,filebase_ps
     &  ,ibuffpos_ps,chbuff_ps,chch_ps,chmarker_ps,cline_ps
     &  ,phi_ps,theta_ps ! 3d viewing angles
     &  ,pi_ps !Pi
     &  ,icolor_ps, ired_ps,igreen_ps,iblue_ps !color, if <1, rgb values are used
     &  ,ilinestyle_ps !line style: 1: solid, 2: dashed, 3: dotted, 4: dashed-dotted
     &  ,xmin3d_ps,xmax3d_ps
     &  ,ymin3d_ps,ymax3d_ps
     &  ,zmin3d_ps,zmax3d_ps
     &  ,offgtitx_ps !offset of global title on pad
     &  ,offgtity_ps !offset of global title on pad
     &  ,offdatex_ps !offset of date on pad
     &  ,offdatey_ps !offset of date on pad
     &  ,xmgl_ps,xmgr_ps,ymgu_ps,ymgl_ps !old hplot variables for margins
     &  ,xwin_ps,ywin_ps ! horizontal and vertical distance of pads
     &  ,xsizorig_ps,ysizorig_ps !size on pad of plotting canvas for paper scale
     &  ,xleftorig_ps,xrightorig_ps,ybottomorig_ps,ytoporig_ps ! orig pad on paper
     &  ,kzone_ps, itouched_ps !current zone and flag of usage
     &  ,nzone_ps,nxzone_ps,nyzone_ps !zones of canvas
     &  ,ifirst_ps !first zone for plotting
     &  ,isamecanvas_ps,isameframe_ps
     &  ,xtit_ps,ytit_ps,ztit_ps,gtit_ps !axis-title and global title
     &  ,nfile_ps !file counter
     &  ,lunall_ps !LUNs of plotfiles
     &  ,fileall_ps !plotfiles
     &  ,kbbxl_ps,kbbxr_ps,kbbyb_ps,kbbyt_ps !bounding box (kbbxl_ps,kbbxr_ps,...)
     &  ,inewpage_ps !new page
     &  ,kBox_ps,kDate_ps ! to plot surrounding box and date on top of plots
     &  ,ygti_ps !y-offset of global title for hplot
     &  ,gsiz_ps !character size of global title for hplot
     &  ,kFrameColor_ps,kFrameRed_ps,kFrameBlue_ps,kFrameGreen_ps
     &  ,kLineColor_ps,kLineRed_ps,kLineBlue_ps,kLineGreen_ps
     &  ,kFillColor_ps,kFillRed_ps,kFillBlue_ps,kFillGreen_ps
     &  ,kTextColor_ps,kTextRed_ps,kTextBlue_ps,kTextGreen_ps
     &  ,kMarkerColor_ps,kMarkerRed_ps,kMarkerBlue_ps,kMarkerGreen_ps
     &  ,xcen_ps,ycen_ps !center on canvas
     &  ,scaletxt_ps ! Scaling for text entities
     &  ,ihigzmode_ps !HIGZ-compatibility mode
     &  ,idrawgtit_ps !draw global title on top of page
     &  ,tsiz_ps ! axis-title-size in HIGZ-mode
     &  ,lunbase_ps
     &  ,coln_ps,cnLine_ps,cnFill_ps,cnMarker_ps,cnText_ps,cnFrame_ps !color norms
     &  ,log10x_ps,log10y_ps,log10z_ps !flags for logarithmic axis
     &  ,aLabHeight_ps !character height of axis label
     &  ,clips_ps,nclips_ps
     &  ,radtodeg_ps
     &  ,mode3d_ps
     &  ,isbox_ps,isdate_ps,isoffdate_ps
     &  ,xoffexp_ps !xoffset of power term
     &  ,yoffexp_ps !yoffset of power term
     &  ,xpcorn_ps,ypcorn_ps
     &  ,ifbox_ps,xfb_ps,yfb_ps
*KEND.

      real x1,y1,x,y
      integer ic,ir,ig,ib
      character(*) text
      character(2048) cline

      itouched_ps=1

      if (icolor_ps.ne.kTextColor_ps) then
        call mshplt_get_text_color(ic,ir,ig,ib)
        call mshplt_set_text_color(ic,ir,ig,ib)
      endif

      if (log10x_ps.eq.0) then
        x1=xleft_ps+scalex_ps*(x-wxmin_ps)
      else
        x1=xleft_ps+scalex_ps*(alog10(x)-wxmin_ps)
      endif

      if (log10y_ps.eq.0) then
        y1=ybottom_ps+scaley_ps*(y-wymin_ps)
      else
        y1=ybottom_ps+scaley_ps*(alog10(y)-wymin_ps)
      endif

      write(cline,*)
     &  x1,y1,' moveto ',
     &  tang_ps,' rotate ',
     &  '(',text(1:len_trim(text)),') show',
     &  -tang_ps,' rotate'
      call mshplt_fill_buff(cline)

      return
      end
*CMZ :  1.03/01 09/10/2014  14.41.17  by  Michael Scheer
*CMZ :  1.00/01 23/09/2014  19.03.09  by  Michael Scheer
*CMZ :  0.01/03 23/09/2014  09.38.20  by  Michael Scheer
*CMZ :  0.01/02 05/09/2014  15.41.43  by  Michael Scheer
*CMZ :  0.01/00 24/08/2014  14.25.20  by  Michael Scheer
*CMZ :  0.00/04 14/08/2014  12.58.41  by  Michael Scheer
*CMZ :  0.00/02 10/07/2014  15.51.37  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_marker(n,x,y)

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      real x(*),y(*),x1,y1,yoffset,r
      integer n,i,mtyp
      character(2048) cline
      character(4) chmark

      if (mtyp_ps.eq.-9999) return

      itouched_ps=1

      if (mtyp_ps.eq.9) then
        do i=1,n
          r=rmsiz_ps*scaletxt_ps
          call mshplt_circle(x(i),y(i),r)
        enddo
        return
      endif

      mtyp=nint(rmtyp_ps(1,mtyp_ps))
      chmark=chch_ps(mtyp)
      yoffset=rmtyp_ps(2,mtyp_ps)*rmsiz_ps*scaletxt_ps

      write(cline,'(a,f12.5,a)')'/Symbol findfont ',rmsiz_ps*scaletxt_ps,
     &  ' scalefont setfont'
      call mshplt_fill_buff(cline)

      do i=1,n

        if (log10x_ps.eq.0) then
          x1=xleft_ps+scalex_ps*(x(i)-wxmin_ps)
        else
          x1=xleft_ps+scalex_ps*(alog10(x(i))-wxmin_ps)
        endif
        if (log10y_ps.eq.0) then
          y1=ybottom_ps+scaley_ps*(y(i)-wymin_ps)
        else
          y1=ybottom_ps+scaley_ps*(alog10(y(i))-wymin_ps)
        endif

        write(cline,'(2f12.5,a,f12.5,a)')x1,y1,
     &    ' moveto ('//chmark(1:len_trim(chmark))
     &    //') dup stringwidth pop 2 div neg 0 rmoveto 0 ',
     &    yoffset,' rmoveto show'

        call mshplt_fill_buff(cline)

      enddo

      write(cline,'(a,f12.5,a)')'/Helvetica findfont ',chhe_ps,
     &  ' scalefont setfont'
      call mshplt_fill_buff(cline)

      return
      end
*CMZ :  1.02/00 30/09/2014  21.25.28  by  Michael Scheer
*CMZ :  1.01/02 25/09/2014  16.02.25  by  Michael Scheer
*CMZ :  0.01/02 05/09/2014  15.41.43  by  Michael Scheer
*CMZ :  0.00/06 20/08/2014  12.06.59  by  Michael Scheer
*CMZ :  0.00/03 04/08/2014  13.33.24  by  Michael Scheer
*CMZ :  0.00/02 09/07/2014  16.03.13  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_3dto2d(n,x,y,z,xp,yp)

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      real x(*),y(*),z(*),xp(*),yp(*)

      real t(3,3),tx(3,3),tz(3,3),
     &  cosphi,sinphi,costhe,sinthe,phio,theo,phi,the
      integer i,j,k,m,n

      data phio/-9999./
      data theo/-9999./

      save

      if (phi_ps.ne.phio.or.theta_ps.ne.theo) then

        phi=phi_ps/radtodeg_ps
        the=theta_ps/radtodeg_ps

        sinphi=sin(phi)
        cosphi=cos(phi)
c        sinthe=sin(the)
c        costhe=cos(the)
        if (theta_ps.gt.90.) then
          sinthe=cos(the) ! we use hplot-convention here
          costhe=-sin(the)
        else
          sinthe=cos(the) ! we use hplot-convention here
          costhe=sin(the)
        endif

        !rotate around z-axis

        tz(1,1)=cosphi
        tz(1,2)=-sinphi
        tz(1,3)=0.0

        tz(2,1)=sinphi
        tz(2,2)=cosphi
        tz(2,3)=0.0

        tz(3,1)=0.0
        tz(3,2)=0.0
        tz(3,3)=1.0

        !rotate around x-axis

        tx(1,1)=1.0
        tx(1,2)=0.0
        tx(1,3)=0.0

        tx(2,1)=0.0
        tx(2,2)=costhe
        tx(2,3)=sinthe

        tx(3,1)=0.0
        tx(3,2)=-sinthe
        tx(3,3)=costhe

        t=0.0
        do i=1,3
          do j=1,3
            do k=1,3
              t(i,j)=t(i,j)+tx(i,k)*tz(k,j)
            enddo
          enddo
        enddo

        phio=phi_ps
        theo=theta_ps

      endif

      do m=1,n
        xp(m)=t(1,1)*x(m)+t(1,2)*y(m)+t(1,3)*z(m)
        yp(m)=t(2,1)*x(m)+t(2,2)*y(m)+t(2,3)*z(m)
      enddo

      return
      end
*CMZ :  1.02/00 01/10/2014  08.23.56  by  Michael Scheer
*CMZ :  0.01/02 05/09/2014  15.41.43  by  Michael Scheer
*CMZ :  0.00/03 11/07/2014  12.27.06  by  Michael Scheer
*CMZ :  0.00/02 07/07/2014  12.22.09  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_set_theta_phi(theta,phi)

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      real theta,phi,p,t

      t=mod(theta,180.)
      if (t.lt.0.0) t=t+180.

      p=mod(phi,360.)
      if (p.lt.0.0) p=p+360.

      phi_ps=p
      theta_ps=t

      return
      end
*CMZ :  1.01/02 26/09/2014  13.59.31  by  Michael Scheer
*CMZ :  0.01/03 22/09/2014  21.04.09  by  Michael Scheer
*CMZ :  0.01/02 11/09/2014  12.03.40  by  Michael Scheer
*CMZ :  0.00/04 14/08/2014  12.58.41  by  Michael Scheer
*CMZ :  0.00/03 06/08/2014  09.32.19  by  Michael Scheer
*CMZ :  0.00/02 09/07/2014  14.48.35  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_pline(n,x,y)

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      real x(*),y(*),x1,y1,x2,y2
      integer n,i,ic,ir,ig,ib,kc,kr,kg,kb
      character(2048) cline

      if (n.le.1) return

      call mshplt_get_line_color(ic,ir,ig,ib)
      call mshplt_get_color(kc,kr,kg,kb)
      if (ic.ne.kc.or.ir.ne.kr.or.kg.ne.ig.or.ib.ne.kb) then
        call mshplt_set_line_color(ic,ir,ig,ib)
      endif

      if (log10x_ps.eq.0.or.mode3d_ps.ne.0) then
        x1=xleft_ps+scalex_ps*(x(1)-wxmin_ps)
      else
        x1=xleft_ps+scalex_ps*(alog10(x(1))-wxmin_ps)
      endif

      if (log10y_ps.eq.0.or.mode3d_ps.ne.0) then
        y1=ybottom_ps+scaley_ps*(y(1)-wymin_ps)
      else
        y1=ybottom_ps+scaley_ps*(alog10(y(1))-wymin_ps)
      endif

      write(cline,*)x1,y1,' moveto'
      call mshplt_fill_buff(cline)

      do i=2,n
        if (log10x_ps.eq.0.or.mode3d_ps.ne.0) then
          x2=xleft_ps+scalex_ps*(x(i)-wxmin_ps)
        else
          x2=xleft_ps+scalex_ps*(alog10(x(i))-wxmin_ps)
        endif
        if (log10y_ps.eq.0.or.mode3d_ps.ne.0) then
          y2=ybottom_ps+scaley_ps*(y(i)-wymin_ps)
        else
          y2=ybottom_ps+scaley_ps*(alog10(y(i))-wymin_ps)
        endif
        write(cline,*)x2,y2,' lineto'
        call mshplt_fill_buff(cline)
      enddo

      itouched_ps=1

      write(cline,*)'stroke'
      call mshplt_fill_buff(cline)

      return
      end
*CMZ :  0.01/03 22/09/2014  20.54.58  by  Michael Scheer
*CMZ :  0.01/02 05/09/2014  15.41.43  by  Michael Scheer
*CMZ :  0.00/03 06/08/2014  09.18.49  by  Michael Scheer
*CMZ :  0.00/02 09/07/2014  14.48.35  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_pline_closed(n,x,y)

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      integer n
      real x(*),y(*)

      if (n.le.1) return

      call mshplt_pline(n,x,y)
      call mshplt_line_raw(x(n),y(n),x(1),y(1))

      return
      end
*CMZ :  1.01/00 24/09/2014  14.35.22  by  Michael Scheer
*CMZ :  0.01/02 05/09/2014  15.41.43  by  Michael Scheer
*CMZ :  0.00/03 06/08/2014  11.12.41  by  Michael Scheer
*CMZ :  0.00/02 09/07/2014  16.03.13  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_pline3d(n,x,y,z)

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      real x(*),y(*),z(*)
      real, dimension(:), allocatable :: xn,yn,zn,xx,yy
      integer n,i

      allocate(xn(n))
      allocate(yn(n))
      allocate(zn(n))
      allocate(xx(n))
      allocate(yy(n))

      do i=1,n
        xn(i)=-0.5+(x(i)-xmin3d_ps)/(xmax3d_ps-xmin3d_ps)
        yn(i)=-0.5+(y(i)-ymin3d_ps)/(ymax3d_ps-ymin3d_ps)
        if (log10z_ps.eq.0) then
          zn(i)=-0.5+(z(i)-zmin3d_ps)/(zmax3d_ps-zmin3d_ps)
        else
          zn(i)=-0.5+(alog10(z(i))-zmin3d_ps)/(zmax3d_ps-zmin3d_ps)
        endif
      enddo

      call mshplt_3dto2d(n,xn,yn,zn,xx,yy)
      call mshplt_pline(n,xx,yy)

      deallocate(xn)
      deallocate(yn)
      deallocate(zn)
      deallocate(xx)
      deallocate(yy)

      return
      end
*CMZ :  0.01/02 05/09/2014  15.41.43  by  Michael Scheer
*CMZ :  0.00/03 06/08/2014  11.30.38  by  Michael Scheer
*CMZ :  0.00/02 09/07/2014  16.03.13  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_pline3d_closed(n,x,y,z)

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      integer n
      real x(*),y(*),z(*),xx(2),yy(2),zz(2)

      call mshplt_pline3d(n,x,y,z)

      xx(1)=x(n)
      yy(1)=y(n)
      zz(1)=z(n)

      xx(2)=x(1)
      yy(2)=y(1)
      zz(2)=z(1)

      call mshplt_pline3d(2,xx,yy,zz)

      return
      end
*CMZ :  0.01/02 05/09/2014  15.41.43  by  Michael Scheer
*CMZ :  0.00/04 06/08/2014  14.52.15  by  Michael Scheer
*CMZ :  0.00/02 07/07/2014  12.22.09  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_set_viewer_interactive(iinter)

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      integer iinter

      iviewinter_ps=iinter

      return
      end
*CMZ :  0.01/02 05/09/2014  15.41.43  by  Michael Scheer
*CMZ :  0.00/04 07/08/2014  08.43.21  by  Michael Scheer
*CMZ :  0.00/02 07/07/2014  12.22.09  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_set_pad_size(xsiz,ysiz)

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      real xsiz,ysiz

      xsiz_ps=xsiz
      ysiz_ps=ysiz

      xleft_ps=(xleft_ps+xright_ps)/2.-xsiz/2.
      ybottom_ps=(ytop_ps+ybottom_ps)/2.-ysiz/2.

      xright_ps=xleft_ps+xsiz
      ytop_ps=ybottom_ps+ysiz

      return
      end
*CMZ :  0.01/02 05/09/2014  15.41.43  by  Michael Scheer
*CMZ :  0.00/04 07/08/2014  08.43.21  by  Michael Scheer
*CMZ :  0.00/02 07/07/2014  12.22.09  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_get_pad_size(xsiz,ysiz)

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      real xsiz,ysiz

      xsiz=xsiz_ps
      ysiz=ysiz_ps

      return
      end
*CMZ :  0.01/02 05/09/2014  15.41.43  by  Michael Scheer
*CMZ :  0.00/04 07/08/2014  08.45.57  by  Michael Scheer
*CMZ :  0.00/02 07/07/2014  12.22.09  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_get_pad(xleft,xright,ybottom,ytop)

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      real xleft,xright,ybottom,ytop

      xleft=xleft_ps
      xright=xright_ps
      ybottom=ybottom_ps
      ytop=ytop_ps

      return
      end
*CMZ :  0.01/02 18/09/2014  12.57.25  by  Michael Scheer
*CMZ :  0.01/01 27/08/2014  09.42.51  by  Michael Scheer
*CMZ :  0.00/05 18/08/2014  12.37.32  by  Michael Scheer
*CMZ :  0.00/04 07/08/2014  08.45.57  by  Michael Scheer
*CMZ :  0.00/02 07/07/2014  12.22.09  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_set_pad(xleft,xright,ybottom,ytop)

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      real xleft,xright,ybottom,ytop

      write(cline_ps,*)
     &  '% begin of mshplt_set_pad xleft, xright,ybottom, ytop:'
      call mshplt_fill_buff(cline_ps(2:len_trim(cline_ps)))
      write(cline_ps,*)
     &  '%',xleft, xright,ybottom, ytop
      call mshplt_fill_buff(cline_ps(2:len_trim(cline_ps)))

      call mshplt_fill_buff('initclip') !reset clipping

      xleft_ps=xleft
      xright_ps=xright
      ybottom_ps=ybottom
      ytop_ps=ytop

      xsiz_ps=xright-xleft
      ysiz_ps=ytop-ybottom

      write(cline_ps,*)
     &  '% end of mshplt_set_pad xleft, xright,ybottom, ytop:'
      call mshplt_fill_buff(cline_ps(2:len_trim(cline_ps)))
      write(cline_ps,*)
     &  '%',xleft, xright,ybottom, ytop
      call mshplt_fill_buff(cline_ps(2:len_trim(cline_ps)))

      return
      end
*CMZ :  0.01/03 22/09/2014  20.54.58  by  Michael Scheer
*CMZ :  0.01/02 05/09/2014  15.41.43  by  Michael Scheer
*CMZ :  0.00/04 07/08/2014  15.19.17  by  Michael Scheer
*CMZ :  0.00/02 07/07/2014  12.22.09  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_draw_pad_box

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      call mshplt_line_raw(0.,0.,xsiz_ps,0.)
      call mshplt_line_raw(xsiz_ps,0.,xsiz_ps,ysiz_ps)
      call mshplt_line_raw(xsiz_ps,ysiz_ps,0.,ysiz_ps)
      call mshplt_line_raw(0.,ysiz_ps,0.,0.)

      return
      end
*CMZ :  0.01/02 05/09/2014  15.41.43  by  Michael Scheer
*CMZ :  0.01/00 23/08/2014  10.04.36  by  Michael Scheer
*CMZ :  0.00/06 21/08/2014  16.34.08  by  Michael Scheer
*CMZ :  0.00/04 07/08/2014  11.34.50  by  Michael Scheer
*CMZ :  0.00/02 07/07/2014  12.22.09  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_title(title)

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      character(*) title

      gtit_ps=title
      idrawgtit_ps=1

      return
      end
*CMZ :  1.01/00 24/09/2014  14.08.29  by  Michael Scheer
*CMZ :  0.01/03 22/09/2014  21.17.47  by  Michael Scheer
*CMZ :  0.01/02 05/09/2014  15.41.43  by  Michael Scheer
*CMZ :  0.00/06 21/08/2014  16.21.19  by  Michael Scheer
*CMZ :  0.00/04 07/08/2014  11.13.32  by  Michael Scheer
*CMZ :  0.00/02 07/07/2014  12.22.09  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_date(chdate)

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      character(*) chdate
      character(10) dtday,dttime,dtzone
      integer idatetime(8)

      integer lendate

      lendate=len_trim(chdate)
      if (lendate.eq.0) return

      if (chdate.eq.'today') then
        call date_and_time(dtday,dttime,dtzone,idatetime)
        write(cline_ps,*)dtday(1:4),'.',dtday(5:6),'.',dtday(7:8),'  ',
     &    dttime(1:2),':',dttime(3:4)
      else
        cline_ps=chdate(1:lendate)
      endif

      lendate=len_trim(cline_ps)

      call mshplt_text_raw(
     &  xrightorig_ps-chhe_ps*scaletxt_ps*lendate+offdatex_ps*scaletxt_ps,
     &  ytoporig_ps-offdatey_ps*scaletxt_ps,
     &  cline_ps(1:lendate)
     &  )

      return
      end
*CMZ :  0.01/02 05/09/2014  15.41.43  by  Michael Scheer
*CMZ :  0.00/06 21/08/2014  16.58.38  by  Michael Scheer
*CMZ :  0.00/04 07/08/2014  11.18.06  by  Michael Scheer
*CMZ :  0.00/02 07/07/2014  12.22.09  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_set_title_offset(x,y)

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      real x,y

      offgtitx_ps=x
      offgtity_ps=y

      return
      end
*CMZ :  0.01/02 05/09/2014  15.41.43  by  Michael Scheer
*CMZ :  0.00/06 21/08/2014  16.58.38  by  Michael Scheer
*CMZ :  0.00/04 07/08/2014  11.19.17  by  Michael Scheer
*CMZ :  0.00/02 07/07/2014  12.22.09  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_get_title_offset(x,y)

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      real x,y

      x=offgtitx_ps
      y=offgtity_ps

      return
      end
*CMZ :  0.01/02 05/09/2014  15.41.43  by  Michael Scheer
*CMZ :  0.00/04 07/08/2014  11.26.55  by  Michael Scheer
*CMZ :  0.00/02 07/07/2014  12.22.09  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_get_axis_title_offset(x,y)

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      real x,y

      x=xtitoff_ps
      y=ytitoff_ps

      return
      end
*CMZ :  0.01/02 05/09/2014  15.41.43  by  Michael Scheer
*CMZ :  0.00/04 07/08/2014  11.26.55  by  Michael Scheer
*CMZ :  0.00/02 07/07/2014  12.22.09  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_set_axis_title_offset(x,y)

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      real x,y

      xtitoff_ps=x
      ytitoff_ps=y

      return
      end
*CMZ :  0.01/02 05/09/2014  15.41.43  by  Michael Scheer
*CMZ :  0.00/06 22/08/2014  15.23.52  by  Michael Scheer
*CMZ :  0.00/04 07/08/2014  11.49.49  by  Michael Scheer
*CMZ :  0.00/02 10/07/2014  15.55.54  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_set_character_height(siz)

      implicit none

      real siz

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      chhe_ps=siz

      if (ihigzmode_ps.eq.0) then
        write(cline_ps,*)'/Helvetica findfont',chhe_ps*scaletxt_ps,
     &    ' scalefont setfont'
      else
        write(cline_ps,*)'/Helvetica findfont',chhe_ps*scaletxt_ps*1.3,
     &    ' scalefont setfont'
      endif

      call mshplt_fill_buff(cline_ps)

      return
      end
*CMZ :  0.01/02 05/09/2014  15.41.43  by  Michael Scheer
*CMZ :  0.00/04 07/08/2014  11.56.06  by  Michael Scheer
*CMZ :  0.00/02 07/07/2014  12.22.09  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_get_axis_label_offset(x,y)

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      real x,y

      x=xlaboff_ps
      y=ylaboff_ps

      return
      end
*CMZ :  0.01/02 05/09/2014  15.41.43  by  Michael Scheer
*CMZ :  0.00/06 21/08/2014  16.54.44  by  Michael Scheer
*CMZ :  0.00/04 07/08/2014  11.55.55  by  Michael Scheer
*CMZ :  0.00/02 07/07/2014  12.22.09  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_set_axis_label_offset(x,y)

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      real x,y

      xlaboff_ps=x
      ylaboff_ps=y

      return
      end
*CMZ :  1.03/02 25/04/2016  15.48.20  by  Michael Scheer
*CMZ :  1.03/01 10/10/2014  13.21.20  by  Michael Scheer
*CMZ :  0.01/03 23/09/2014  12.06.56  by  Michael Scheer
*CMZ :  0.01/02 05/09/2014  15.41.43  by  Michael Scheer
*CMZ :  0.01/00 23/08/2014  11.47.48  by  Michael Scheer
*CMZ :  0.00/04 14/08/2014  12.58.41  by  Michael Scheer
*CMZ :  0.00/02 09/07/2014  16.17.24  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_text_ndc(x,y,text)

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      real x1,y1,x,y
      integer ic,ir,ig,ib
      character(*) text
      character(2048) cline

      itouched_ps=1

      if (icolor_ps.ne.kTextColor_ps) then
        call mshplt_get_text_color(ic,ir,ig,ib)
        call mshplt_set_text_color(ic,ir,ig,ib)
      endif

      x1=xleft_ps+x*xsiz_ps
      y1=ybottom_ps+y*ysiz_ps

      write(cline,*) x1,y1,
     &  ' moveto ', tang_ps,' rotate (',text(1:len_trim(text)),') show ',
     &  -tang_ps,' rotate'
      call mshplt_fill_buff(cline)

      return
      end
*CMZ :  0.01/02 05/09/2014  15.41.43  by  Michael Scheer
*CMZ :  0.00/06 19/08/2014  15.09.19  by  Michael Scheer
*CMZ :  0.00/04 07/08/2014  15.52.16  by  Michael Scheer
*CMZ :  0.00/03 04/08/2014  15.42.32  by  Michael Scheer
*CMZ :  0.00/02 07/07/2014  12.22.09  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_set_line_color(icolor,ired,igreen,iblue)

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      integer icolor,ired,igreen,iblue

      kLineColor_ps=icolor

      if (icolor.gt.0) then
        if (icolor.eq.1) then
          kLineRed_ps=0
          kLineGreen_ps=0
          kLineBlue_ps=0
        else if (icolor.eq.2) then
          kLineRed_ps=1
          kLineGreen_ps=0
          kLineBlue_ps=0
        else if (icolor.eq.3) then
          kLineRed_ps=0
          kLineGreen_ps=1
          kLineBlue_ps=0
        else if (icolor.eq.4) then
          kLineRed_ps=0
          kLineGreen_ps=0
          kLineBlue_ps=1
        else if (icolor.eq.5) then
          kLineRed_ps=1
          kLineGreen_ps=1
          kLineBlue_ps=0
        else if (icolor.eq.6) then
          kLineRed_ps=1
          kLineGreen_ps=0
          kLineBlue_ps=1
        else if (icolor.eq.7) then
          kLineRed_ps=0
          kLineGreen_ps=1
          kLineBlue_ps=1
        else if (icolor.eq.8) then
          kLineRed_ps=35
          kLineGreen_ps=85
          kLineBlue_ps=33
        endif
      else
        kLineColor_ps=-1
        kLineRed_ps=ired
        kLineGreen_ps=igreen
        kLineBlue_ps=iblue
      endif !icolor.gt.0

      cnLine_ps=sqrt(float(kLineRed_ps**2+kLineGreen_ps**2+kLineBlue_ps**2))
      if (cnLine_ps.le.0.) cnLine_ps=1.

      write(cline_ps,*)
     &  kLineRed_ps/cnLine_ps,
     &  kLineGreen_ps/cnLine_ps,
     &  kLineBlue_ps/cnLine_ps,
     &  ' setrgbcolor'
      call mshplt_fill_buff(cline_ps)

      return
      end
*CMZ :  0.01/02 05/09/2014  15.41.43  by  Michael Scheer
*CMZ :  0.00/04 07/08/2014  15.53.15  by  Michael Scheer
*CMZ :  0.00/03 04/08/2014  16.13.05  by  Michael Scheer
*CMZ :  0.00/02 07/07/2014  12.22.09  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_set_line_style(istyle)

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      integer istyle

      if (istyle.eq.ilinestyle_ps) return

      if (istyle.ge.0.and.istyle.le.4) then

        ilinestyle_ps=istyle

        if (ilinestyle_ps.eq.1) then
          write(cline_ps,*)'[] 0 setdash'
        else if (ilinestyle_ps.eq.2) then
          write(cline_ps,*)'[',6*rlinewidth_ps, 3*rlinewidth_ps,'] 0 setdash'
        else if (ilinestyle_ps.eq.3) then
          write(cline_ps,*)'[',2*rlinewidth_ps, 2*rlinewidth_ps,'] 0 setdash'
        else if (ilinestyle_ps.eq.4) then
          write(cline_ps,*)'[',
     &      8*rlinewidth_ps,
     &      4*rlinewidth_ps,
     &      2*rlinewidth_ps,
     &      4*rlinewidth_ps,
     &      '] 0 setdash'
        endif

        call mshplt_fill_buff(cline_ps)

      endif

      return
      end
*CMZ :  1.03/01 10/10/2014  08.54.46  by  Michael Scheer
*CMZ :  0.01/02 05/09/2014  15.41.43  by  Michael Scheer
*CMZ :  0.01/01 27/08/2014  11.57.12  by  Michael Scheer
*CMZ :  0.00/04 14/08/2014  12.02.07  by  Michael Scheer
*CMZ :  0.00/03 04/08/2014  15.28.51  by  Michael Scheer
*CMZ :  0.00/02 09/07/2014  11.35.10  by  Michael Scheer
*CMZ : 00.00/02 30/06/2014  10.16.31  by  Michael Scheer
*-- Author :    Michael Scheer   27/06/2014
      subroutine mshplt_set_scale(scale)

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      real scale

      scale_ps=scale

      if (isscale_ps.ne.1) then
        isscale_ps=1
      else
        write(cline_ps,*)scale_ps,scale_ps,' scale'
        call mshplt_fill_buff(cline_ps)
      endif

      return
      end
*CMZ :  0.01/02 05/09/2014  15.41.43  by  Michael Scheer
*CMZ :  0.00/04 07/08/2014  16.06.35  by  Michael Scheer
*CMZ :  0.00/02 07/07/2014  12.22.09  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_get_frame(xmin,xmax,ymin,ymax)

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      real xmin,xmax,ymin,ymax

      xmin=wxmin_ps
      xmax=wxmax_ps
      ymin=wymin_ps
      ymax=wymax_ps

      return
      end
*CMZ :  0.01/02 05/09/2014  15.41.43  by  Michael Scheer
*CMZ :  0.00/04 07/08/2014  16.15.06  by  Michael Scheer
*CMZ :  0.00/02 10/07/2014  15.55.54  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_set_tic_size(siz)

      implicit none

      real siz

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      ticsiz_ps=siz

      return
      end
*CMZ :  1.03/01 09/10/2014  14.42.41  by  Michael Scheer
*CMZ :  0.01/02 05/09/2014  15.41.43  by  Michael Scheer
*CMZ :  0.01/00 24/08/2014  14.21.01  by  Michael Scheer
*CMZ :  0.00/04 08/08/2014  16.50.25  by  Michael Scheer
*CMZ :  0.00/02 10/07/2014  16.13.29  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_set_marker_type(ktyp)

      implicit none

      integer mtyp,ktyp

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      mtyp=ktyp

      if (mtyp.eq.-9999) goto 9999

      if (ktyp.eq.20) then
        mtyp=1
      else if (ktyp.eq.24) then
        mtyp=9
      else if (ktyp.eq.31) then
        mtyp=5
      endif

      if (mtyp.lt.0.or.mtyp.gt.9) then
        write(6,*) '*** Error in mshplt_set_marker_type: Range of valid marker types is 0 -> 9'
        return
      endif

      if (mtyp.eq.9) then
        chmarker_ps=chch_ps(nint(rmtyp_ps(1,1)))
      endif

9999  mtyp_ps=mtyp
      return
      end
*CMZ :  1.02/00 03/10/2014  10.00.50  by  Michael Scheer
*CMZ :  0.01/02 05/09/2014  15.41.43  by  Michael Scheer
*CMZ :  0.00/04 08/08/2014  15.00.08  by  Michael Scheer
*CMZ :  0.00/02 10/07/2014  15.21.58  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_set_marker_size(siz)

      implicit none

      real siz

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      rmsiz_ps=siz

      return
      end
*CMZ :  1.02/01 03/10/2014  14.42.26  by  Michael Scheer
*CMZ :  0.01/02 05/09/2014  15.41.43  by  Michael Scheer
*CMZ :  0.00/06 19/08/2014  15.10.09  by  Michael Scheer
*CMZ :  0.00/04 08/08/2014  16.20.54  by  Michael Scheer
*CMZ :  0.00/03 04/08/2014  15.42.32  by  Michael Scheer
*CMZ :  0.00/02 07/07/2014  12.22.09  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_set_color(icolor,ired,igreen,iblue)

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      integer icolor,ired,igreen,iblue

      if (icolor.gt.0) then
        if (icolor.eq.1) then
          ired_ps=0
          igreen_ps=0
          iblue_ps=0
        else if (icolor.eq.2) then
          ired_ps=1
          igreen_ps=0
          iblue_ps=0
        else if (icolor.eq.3) then
          ired_ps=0
          igreen_ps=1
          iblue_ps=0
        else if (icolor.eq.4) then
          ired_ps=0
          igreen_ps=0
          iblue_ps=1
        else if (icolor.eq.5) then
          ired_ps=1
          igreen_ps=1
          iblue_ps=0
        else if (icolor.eq.6) then
          ired_ps=1
          igreen_ps=0
          iblue_ps=1
        else if (icolor.eq.7) then
          ired_ps=0
          igreen_ps=1
          iblue_ps=1
        else if (icolor.eq.8) then
          ired_ps=35
          igreen_ps=85
          iblue_ps=33
        endif
      else
        icolor_ps=-1
        ired_ps=ired
        igreen_ps=igreen
        iblue_ps=iblue
      endif

      coln_ps=sqrt(float(ired_ps**2+igreen_ps**2+iblue_ps**2))
      if (icolor_ps.lt.0.or.coln_ps.le.0.0) coln_ps=1.

      write(cline_ps,*)
     &  ired_ps/coln_ps,
     &  igreen_ps/coln_ps,
     &  iblue_ps/coln_ps,
     &  ' setrgbcolor'
      call mshplt_fill_buff(cline_ps(1:len_trim(cline_ps)))

      return
      end
*CMZ :  0.01/02 05/09/2014  15.41.43  by  Michael Scheer
*CMZ :  0.00/06 19/08/2014  15.01.54  by  Michael Scheer
*CMZ :  0.00/04 08/08/2014  16.34.41  by  Michael Scheer
*CMZ :  0.00/03 04/08/2014  15.42.32  by  Michael Scheer
*CMZ :  0.00/02 07/07/2014  12.22.09  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_get_line_color(icolor,ired,igreen,iblue)

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      integer icolor,ired,igreen,iblue

      icolor=kLineColor_ps
      ired=kLineRed_ps
      igreen=kLineGreen_ps
      iblue=kLineBlue_ps

      return
      end
*CMZ :  0.01/02 05/09/2014  15.41.43  by  Michael Scheer
*CMZ :  0.00/04 08/08/2014  16.34.07  by  Michael Scheer
*CMZ :  0.00/03 04/08/2014  15.42.32  by  Michael Scheer
*CMZ :  0.00/02 07/07/2014  12.22.09  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_get_color(icolor,ired,igreen,iblue)

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      integer icolor,ired,igreen,iblue

      icolor=icolor_ps
      ired=ired_ps
      igreen=igreen_ps
      iblue=iblue_ps

      return
      end
*CMZ :  0.01/03 23/09/2014  09.22.45  by  Michael Scheer
*CMZ :  0.01/02 15/09/2014  15.09.38  by  Michael Scheer
*CMZ :  0.00/04 11/08/2014  18.23.32  by  Michael Scheer
*CMZ :  0.00/02 07/07/2014  12.22.09  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_circle(x,y,r)

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      real x,y,r,x1,y1

        if (log10x_ps.eq.0) then
          x1=xleft_ps+scalex_ps*(x-wxmin_ps)
        else
          x1=xleft_ps+scalex_ps*(alog10(x)-wxmin_ps)
        endif
        if (log10y_ps.eq.0) then
          y1=ybottom_ps+scaley_ps*(y-wymin_ps)
        else
          y1=ybottom_ps+scaley_ps*(alog10(y)-wymin_ps)
        endif

      write(cline_ps,*)'newpath ',x1,y1,r*2.*chhe_ps,' 0 360 arc stroke'
      call mshplt_fill_buff(cline_ps)

      return
      end
*CMZ :  0.01/02 05/09/2014  15.41.43  by  Michael Scheer
*CMZ :  0.00/04 11/08/2014  18.31.12  by  Michael Scheer
*CMZ :  0.00/02 07/07/2014  12.22.09  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_scale_character_height(fac)

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      real fac,chhe

      chhe=chhe_ps
      call mshplt_set_character_height(fac*chhe)

      return
      end
*CMZ :  0.01/02 05/09/2014  15.41.43  by  Michael Scheer
*CMZ :  0.00/04 11/08/2014  18.38.59  by  Michael Scheer
*CMZ :  0.00/02 10/07/2014  15.21.58  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_scale_marker_size(fac)

      implicit none

      real fac

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      rmsiz_ps=rmsiz_ps*fac

      return
      end
*CMZ :  0.01/02 05/09/2014  15.41.43  by  Michael Scheer
*CMZ :  0.00/04 11/08/2014  20.51.17  by  Michael Scheer
*CMZ :  0.00/02 10/07/2014  15.21.58  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_get_marker_size(siz)

      implicit none

      real siz

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      siz=rmsiz_ps

      return
      end
*CMZ :  0.01/03 23/09/2014  12.00.19  by  Michael Scheer
*CMZ :  0.01/02 05/09/2014  15.41.43  by  Michael Scheer
*CMZ :  0.00/04 14/08/2014  12.58.41  by  Michael Scheer
*CMZ :  0.00/03 06/08/2014  09.32.19  by  Michael Scheer
*CMZ :  0.00/02 09/07/2014  14.48.35  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_clip(n,x,y)

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

c Controls clipping:
c      n=0: return -nclips_ps
c      n<0: Activate clipping for already defined clips_ps
c      n>2: Set clipping for (x,y) path
c     else: return

c see also: mshplt_get_clip

      real x(*),y(*),x1,y1,x2,y2
      integer n,i
      character(2048) cline

      if (n.gt.0.and.n.le.2) then
        nclips_ps=0
        return
      endif

      write(cline,*)'initclip' !reset clipping
      call mshplt_fill_buff(cline)

      if (n.eq.0) then
        n=-nclips_ps
        return
      endif

c Just to switch clipping on for already defined x,y
      if (n.lt.0) then

        nclips_ps=-n

        if (n.gt.1000) then
          Print*,'*** Error in mshplt_clip: Dimension exceeded ***'
          return
        endif

        x1=xleft_ps+scalex_ps*(clips_ps(1,1)-wxmin_ps)
        y1=ybottom_ps+scaley_ps*(clips_ps(2,1)-wymin_ps)

        write(cline,*)'newpath ', x1,y1,' moveto'
        call mshplt_fill_buff(cline)

        do i=2,nclips_ps
          x2=xleft_ps+scalex_ps*(clips_ps(1,i)-wxmin_ps)
          y2=ybottom_ps+scaley_ps*(clips_ps(2,i)-wymin_ps)
          write(cline,*)x1,y1,x2,y2,' lineto'
          call mshplt_fill_buff(cline)
          x1=x2
          y1=y2
        enddo

        write(cline,*)'closepath clip newpath'
        call mshplt_fill_buff(cline)

      else ! if (n.lt.0) then

        nclips_ps=n

        if (n.gt.1000) then
          Print*,'*** Error in mshplt_clip: Dimension exceeded ***'
          return
        endif

        nclips_ps=n
        clips_ps(1,1)=x(1)
        clips_ps(2,1)=y(1)

        x1=xleft_ps+scalex_ps*(x(1)-wxmin_ps)
        y1=ybottom_ps+scaley_ps*(y(1)-wymin_ps)

        write(cline,*)'newpath ', x1,y1,' moveto'
        call mshplt_fill_buff(cline)

        do i=2,nclips_ps
          clips_ps(1,i)=x(i)
          clips_ps(2,i)=y(i)
          x2=xleft_ps+scalex_ps*(x(i)-wxmin_ps)
          y2=ybottom_ps+scaley_ps*(y(i)-wymin_ps)
          write(cline,*)x1,y1,x2,y2,' lineto'
          call mshplt_fill_buff(cline)
          x1=x2
          y1=y2
        enddo

        write(cline,*)'closepath clip newpath'
        call mshplt_fill_buff(cline)

      endif ! if (n.lt.0) then

      return
      end
*CMZ :  0.01/02 05/09/2014  15.41.43  by  Michael Scheer
*CMZ :  0.00/04 14/08/2014  12.58.41  by  Michael Scheer
*CMZ :  0.00/02 07/07/2014  12.22.09  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_set_eps_file(file)

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      character(2048) file

      file_ps=file(1:len_trim(file))

      return
      end
*CMZ :          27/10/2017  08.39.09  by  Michael Scheer
*CMZ :  0.01/02 05/09/2014  15.41.43  by  Michael Scheer
*CMZ :  0.01/00 23/08/2014  09.16.36  by  Michael Scheer
*CMZ :  0.00/06 22/08/2014  15.16.33  by  Michael Scheer
*CMZ :  0.00/05 18/08/2014  14.45.59  by  Michael Scheer
*CMZ :  0.00/04 18/08/2014  10.16.50  by  Michael Scheer
*CMZ :  0.00/02 07/07/2014  12.22.09  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_zone(nxzon,nyzon,ifirst,chopt)

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      integer nxzon,nyzon,ifirst,ix,iy,i
      real wid,hig,xl,xr,yl,yt
      character(*) chopt
      character choptd
      integer ic
      equivalence (ic,choptd)

      call mshplt_flush_buff

      choptd=chopt
      if (len_trim(chopt).ne.0.and.ic.ne.0) then
        write(cline_ps,*)'% begin of mshplt_zone',nxzon,nyzon,ifirst,chopt
      else
        write(cline_ps,*)'% begin of mshplt_zone',nxzon,nyzon,ifirst
      endif
      write(lun_ps,'(a)')cline_ps(2:len_trim(cline_ps))

      nxzone_ps=nxzon
      nyzone_ps=nyzon
      nzone_ps=nxzone_ps*nyzone_ps

      if (chopt.eq.'S'.or.chopt.eq.'s') then
        isamecanvas_ps=1
      else
        isamecanvas_ps=0
        isameframe_ps=0
      endif

      if (nzone_ps.gt.0) then
        if (ifirst.gt.nzone_ps.or.ifirst.le.0) then
          print*,'*** WARNING IN MSHPLT_ZONE: IFIRST > NXZONE*NYZONE OR IFIRST <= 0'
          print*,'*** WARNING IN MSHPLT_ZONE: IFRST SET TO ONE'
          ifirst_ps=1
        else
          ifirst_ps=ifirst
        endif
      endif

      wid=(xsizorig_ps-xmgl_ps-xmgr_ps-xwin_ps*(nxzon-1))/nxzon
      hig=(ysizorig_ps-ymgl_ps-ymgu_ps-ywin_ps*(nyzon-1))/nyzon

      ix=0
      iy=1
      do i=1,ifirst
        ix=ix+1
        if (ix.gt.nxzon) then
          ix=1
          iy=iy+1
        endif
      enddo

      xl=xleftorig_ps+xmgl_ps+(ix-1)*(wid+xwin_ps)
      xr=xl+wid
      yl=ybottomorig_ps+ymgl_ps+(nyzon-iy)*(hig+ywin_ps)
      yt=yl+hig
      call mshplt_set_pad(xl,xr,yl,yt)

      if (isamecanvas_ps.eq.0.and.abs(kzone_ps).ge.nzone_ps
     &    .and.itouched_ps.ne.0) then
        call mshplt_newpage
      endif

      kzone_ps=ifirst
      kzone_ps=min(-1,-abs(kzone_ps)) !to indicate call to mshplt_zone

      write(cline_ps,*)'% end of mshplt_zone'
      write(lun_ps,'(a)')cline_ps(2:len_trim(cline_ps))

      return
      end
*CMZ :  1.01/02 27/09/2014  15.10.10  by  Michael Scheer
*CMZ :  0.01/02 05/09/2014  15.41.43  by  Michael Scheer
*CMZ :  0.00/05 18/08/2014  15.03.41  by  Michael Scheer
*CMZ :  0.00/04 14/08/2014  12.47.46  by  Michael Scheer
*CMZ :  0.00/02 07/07/2014  12.22.09  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_newpage

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      integer ierr
      character(2048) file
      data file/'filedummy'/

      if (itouched_ps.eq.0) return

      call mshplt_flush_buff

      inewpage_ps=1
      itouched_ps=0
      call mshplt_file_open(-1,file(1:len_trim(file)),'',ierr)

      return
      end
*CMZ :  1.03/01 07/10/2014  15.16.02  by  Michael Scheer
*CMZ :  1.02/00 02/10/2014  21.58.12  by  Michael Scheer
*CMZ :  1.01/02 27/09/2014  16.14.01  by  Michael Scheer
*CMZ :  0.01/02 12/09/2014  09.10.17  by  Michael Scheer
*CMZ :  0.01/00 23/08/2014  09.58.55  by  Michael Scheer
*CMZ :  0.00/06 22/08/2014  15.02.41  by  Michael Scheer
*CMZ :  0.00/05 18/08/2014  15.01.02  by  Michael Scheer
*CMZ :  0.00/04 16/08/2014  15.29.25  by  Michael Scheer
*CMZ :  0.00/02 07/07/2014  12.22.09  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_file_open(luni,file,chopt,ierr)

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      character(10) dtday,dttime,dtzone
      integer idatetime(8),lun,luni,ierr,ifound,lunmax,i
      integer lchopt,lfile,l1num,l2num,nlun
      character(*) file,chopt
      character(2048) cline,cnum

      character(8192) fileeps

      if (nfile_ps.ge.nmaxfile_ps) then
        ierr=-1000
        return
      endif

      lun=luni

      write(cnum,*)nfile_ps

      lunmax=0
      do i=1,nfile_ps
        if (lunall_ps(nfile_ps).gt.lunmax) lunmax=lunall_ps(nfile_ps)
      enddo

      ifound=0
      do i=1,nfile_ps
        if (lunall_ps(nfile_ps).eq.luni) ifound=1
      enddo

      if (ifound.eq.1) then
        lun=lunmax+1
        print*,'*** Warning in mshplt_file_open: LUN',luni,'is in use!'
        print*,'*** Changed to',lun
      endif

      nfile_ps=nfile_ps+1

      l2num=len_trim(cnum)
      do l1num=1,l2num
        if (cnum(l1num:l1num).ne.' ') goto 9
      enddo

9     continue

      lfile=len_trim(file)

      if (lun.gt.0) then
        nlun=lun
        if (file(lfile-2:lfile).ne.'eps') then
          filebase_ps=file(1:lfile)
          fileeps=file(1:lfile)//'_'//cnum(l1num:l2num)//'.eps'
        else
          filebase_ps=file(1:lfile-4)
          fileeps=file(1:lfile-4)//'_'//cnum(l1num:l2num)//'.eps'
        endif
      else
        nlun=lun_ps+1
        fileeps=filebase_ps(1:len_trim(filebase_ps))//'_'//cnum(l1num:l2num)//'.eps'
      endif

      if (nfile_ps.eq.1) fileeps=file

      lfile=len_trim(fileeps)
      lchopt=len_trim(chopt)

      open(unit=nlun,file=fileeps(1:lfile),status='unknown',iostat=ierr)
      file_ps=fileeps(1:lfile)

      CALL date_and_time(dtday,dttime,dtzone,idatetime)

      lun_ps=nlun

      write(lun_ps,'(a)')'%!PS-Adobe-2.0'
      write(cline,*)'%%BoundingBox: ',kbbxl_ps,kbbyb_ps,kbbxr_ps,kbbyt_ps

      write(lun_ps,'(a)')cline(2:len_trim(cline))
      write(lun_ps,'(a)')'%%Title: '//file_ps(1:len_trim(file_ps))
      write(lun_ps,'(a)')'%%Pages: (atend)'
c      write(lun_ps,'(a)')'%%Creator: mshplt'
*KEEP,MSHPLTVERSION.
      write(lun_ps,'(a)') '%%Creator mshplt 1.03/02'
*KEND.
      write(cline,*)'%%CreationDate: ',
     &  dtday(1:4),'.',dtday(5:6),'.',dtday(7:8),'   ',
     &  dttime(1:2),':',dttime(3:4),':',dttime(5:6)
      write(lun_ps,'(a)')cline(2:len_trim(cline))
      write(lun_ps,'(a)')'%%EndComments'

      write(lun_ps,*)scale_ps,scale_ps,' scale'
      write(lun_ps,*)rlinewidth_ps,' setlinewidth'
      write(lun_ps,*)xleftorig_ps,ybottomorig_ps,' moveto %goto origin'

      call mshplt_draw_title

      if (ihigzmode_ps.ne.0) then
        write(lun_ps,*)'/Helvetica findfont ',chhe_ps*scaletxt_ps*1.5,
     &    ' scalefont setfont'
      else
        write(lun_ps,*)'/Helvetica findfont ',chhe_ps*scaletxt_ps,
     &    ' scalefont setfont'
      endif

      if (inewpage_ps.ne.0) then

        if (kBox_ps.gt.0) then
          write(lun_ps,*)xsizorig_ps,0.0,' rlineto'
          write(lun_ps,*)0.0,ysizorig_ps,' rlineto'
          write(lun_ps,*)-xsizorig_ps,0.0,' rlineto'
          write(lun_ps,*)0.0,-ysizorig_ps,' rlineto'
          write(lun_ps,*)'stroke'
        endif

        if (ihigzmode_ps.ne.0) then

          if (kDate_ps.gt.0) then
            write(lun_ps,*)xrightorig_ps+offdatex_ps*scaletxt_ps*1.5,
     &        ytoporig_ps+offdatey_ps*scaletxt_ps-1.25,
     &        'moveto ( ',
     &        dtday(1:4),'.',dtday(5:6),'.',dtday(7:8),'   ',
     &        dttime(1:2),':',dttime(3:4),
     &        ') show'
          endif

        else

          if (kDate_ps.gt.0) then
            write(lun_ps,*)xrightorig_ps+offdatex_ps,
     &        ytoporig_ps+offdatey_ps,
     &        'moveto ( ',
     &        dtday(1:4),'.',dtday(5:6),'.',dtday(7:8),'   ',
     &        dttime(1:2),':',dttime(3:4),
     &        ') show'
          endif

        endif

        inewpage_ps=0
      endif

      write(lun_ps,*)xleftorig_ps,ybottomorig_ps,' moveto %goto origin'
      lunall_ps(nfile_ps)=nlun
      fileall_ps(nfile_ps)=fileeps(1:lfile)

      if (irunviewer_ps.ne.0) then
        call system(
     &    viewer_ps(1:len_trim(viewer_ps))//' '//
     &    file_ps(1:len_trim(file_ps)))
      endif

      write(lun_ps,'(a)')'% end of mshplt_file_open'
      return
      end
*CMZ :  1.03/02 01/04/2016  12.31.53  by  Michael Scheer
*CMZ :  0.01/02 05/09/2014  15.41.43  by  Michael Scheer
*CMZ :  0.00/06 19/08/2014  13.50.41  by  Michael Scheer
*CMZ :  0.00/04 11/08/2014  17.53.18  by  Michael Scheer
*CMZ :  0.00/02 07/07/2014  12.22.09  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_hplset(chopt,val)

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      character(*) chopt
      real val

      if (chopt.eq.'YGTI'.or.chopt.eq.'ygti') then
        ygti_ps=val
        offgtity_ps=val
      else if (chopt.eq.'GSIZ'.or.chopt.eq.'gsiz') then
        gsiz_ps=val
      else if (chopt.eq.'YMGL'.or.chopt.eq.'ymgl') then
        ymgl_ps=val
      else if (chopt.eq.'YMGU'.or.chopt.eq.'ymgu') then
        ymgu_ps=val
      else if (chopt.eq.'XMGL'.or.chopt.eq.'xmgl') then
        xmgl_ps=val
      else if (chopt.eq.'XMGR'.or.chopt.eq.'xmgr') then
        xmgr_ps=val
      else if (chopt.eq.'YWIN'.or.chopt.eq.'ywin') then
        ywin_ps=val
      else if (chopt.eq.'YWIN'.or.chopt.eq.'xwin') then
        xwin_ps=val
      else
        print*,'*** WARNING in mshplt_hplset: Option ',
     &    chopt(1:len_trim(chopt)),' not available! ***'
      endif

      return
      end
*CMZ :          25/09/2016  13.35.45  by  Michael Scheer
*CMZ :  1.03/01 09/10/2014  15.40.58  by  Michael Scheer
*CMZ :  0.01/02 05/09/2014  15.41.43  by  Michael Scheer
*CMZ :  0.00/06 19/08/2014  15.54.32  by  Michael Scheer
*CMZ :  0.00/04 11/08/2014  17.53.18  by  Michael Scheer
*CMZ :  0.00/02 07/07/2014  12.22.09  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_hplopt(chopt,n)

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      character(*) chopt
      integer n

      if (n.gt.1) then
        print*,'*** WARNING in mshplt_hplopt: Number of options must not exceed one!'
      endif

      if (chopt.eq.'DATE'.or.chopt.eq.'date') then
        kDate_ps=1
        isDate_ps=1
      else if (chopt.eq.'NDAT'.or.chopt.eq.'ndat') then
        kDate_ps=0
        isDate_ps=1
      else if (chopt.eq.'BOX '.or.chopt.eq.'box ') then
        kBox_ps=1
        isBox_ps=1
      else if (chopt.eq.'NBOX'.or.chopt.eq.'nbox') then
        isBox_ps=1
        kBox_ps=0
      else
        print*,'*** WARNING in mshplt_hplopt: Option ',
     &    chopt(1:len_trim(chopt)),' not available! ***'
      endif

      return
      end
*CMZ :          26/05/2017  10.10.47  by  Michael Scheer
*CMZ :  1.03/02 25/04/2016  09.27.22  by  Michael Scheer
*CMZ :  0.01/02 05/09/2014  15.41.43  by  Michael Scheer
*CMZ :  0.00/06 19/08/2014  14.33.17  by  Michael Scheer
*CMZ :  0.00/04 11/08/2014  17.53.18  by  Michael Scheer
*CMZ :  0.00/02 07/07/2014  12.22.09  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_igset(pname,val)

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      character(*) pname
      real val

      if (pname.eq.'PLCI'.or.pname.eq.'plci') then
        call mshplt_set_line_color(int(val),0,0,0)
      else if (pname.eq.'TXCI'.or.pname.eq.'txci') then
        call mshplt_set_text_color(int(val),0,0,0)
      else if (pname.eq.'PMCI'.or.pname.eq.'pmci') then
        call mshplt_set_marker_color(int(val),0,0,0)
      else if (pname.eq.'MTYP'.or.pname.eq.'mtyp') then
        call mshplt_set_marker_type(int(val))
      else if (pname.eq.'MSCF'.or.pname.eq.'mscf') then
        call mshplt_set_marker_size(val)
      else if (pname.eq.'LWID'.or.pname.eq.'lwid') then
        call mshplt_set_line_width(val)
      else if (pname.eq.'CHHE'.or.pname.eq.'chhe') then
        call mshplt_set_character_height(val)
      else
        print*,'*** WARNING in mshplt_igset: Parameter ',
     &    pname(1:len_trim(pname)),' not available! ***'
      endif

      return
      end
*CMZ :  0.01/02 05/09/2014  16.58.57  by  Michael Scheer
*CMZ :  0.00/06 19/08/2014  15.11.12  by  Michael Scheer
*CMZ :  0.00/04 07/08/2014  15.52.16  by  Michael Scheer
*CMZ :  0.00/03 04/08/2014  15.42.32  by  Michael Scheer
*CMZ :  0.00/02 07/07/2014  12.22.09  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_set_marker_color(icolor,ired,igreen,iblue)

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      integer icolor,ired,igreen,iblue

      if (icolor.gt.0) then
        if (icolor.eq.1) then
          kMarkerRed_ps=0
          kMarkerGreen_ps=0
          kMarkerBlue_ps=0
        else if (icolor.eq.2) then
          kMarkerRed_ps=1
          kMarkerGreen_ps=0
          kMarkerBlue_ps=0
        else if (icolor.eq.3) then
          kMarkerRed_ps=0
          kMarkerGreen_ps=1
          kMarkerBlue_ps=0
        else if (icolor.eq.4) then
          kMarkerRed_ps=0
          kMarkerGreen_ps=0
          kMarkerBlue_ps=1
        else if (icolor.eq.5) then
          kMarkerRed_ps=1
          kMarkerGreen_ps=1
          kMarkerBlue_ps=0
        else if (icolor.eq.6) then
          kMarkerRed_ps=1
          kMarkerGreen_ps=0
          kMarkerBlue_ps=1
        else if (icolor.eq.7) then
          kMarkerRed_ps=0
          kMarkerGreen_ps=1
          kMarkerBlue_ps=1
        else if (icolor.eq.8) then
          kMarkerRed_ps=35
          kMarkerGreen_ps=85
          kMarkerBlue_ps=33
        endif
      else
        kMarkerColor_ps=-1
        kMarkerRed_ps=ired
        kMarkerGreen_ps=igreen
        kMarkerBlue_ps=iblue
      endif !icolor.gt.0

      cnMarker_ps=sqrt(float(kMarkerRed_ps**2+kMarkerGreen_ps**2+kMarkerBlue_ps**2))
      if (cnMarker_ps.le.0) cnMarker_ps=1.

      write(cline_ps,*)
     &  kMarkerRed_ps/cnMarker_ps,
     &  kMarkerGreen_ps/cnMarker_ps,
     &  kMarkerBlue_ps/cnMarker_ps,
     &  ' setrgbcolor'
      call mshplt_fill_buff(cline_ps)

      return
      end
*CMZ :  0.01/02 05/09/2014  16.58.57  by  Michael Scheer
*CMZ :  0.00/06 19/08/2014  15.11.41  by  Michael Scheer
*CMZ :  0.00/04 07/08/2014  15.52.16  by  Michael Scheer
*CMZ :  0.00/03 04/08/2014  15.42.32  by  Michael Scheer
*CMZ :  0.00/02 07/07/2014  12.22.09  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_set_text_color(icolor,ired,igreen,iblue)

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      integer icolor,ired,igreen,iblue

      kTextColor_ps=icolor

      if (icolor.gt.0) then
        if (icolor.eq.1) then
          kTextRed_ps=0
          kTextGreen_ps=0
          kTextBlue_ps=0
        else if (icolor.eq.2) then
          kTextRed_ps=1
          kTextGreen_ps=0
          kTextBlue_ps=0
        else if (icolor.eq.3) then
          kTextRed_ps=0
          kTextGreen_ps=1
          kTextBlue_ps=0
        else if (icolor.eq.4) then
          kTextRed_ps=0
          kTextGreen_ps=0
          kTextBlue_ps=1
        else if (icolor.eq.5) then
          kTextRed_ps=1
          kTextGreen_ps=1
          kTextBlue_ps=0
        else if (icolor.eq.6) then
          kTextRed_ps=1
          kTextGreen_ps=0
          kTextBlue_ps=1
        else if (icolor.eq.7) then
          kTextRed_ps=0
          kTextGreen_ps=1
          kTextBlue_ps=1
        else if (icolor.eq.8) then
          kTextRed_ps=35
          kTextGreen_ps=85
          kTextBlue_ps=33
        endif
      else
        kTextColor_ps=-1
        kTextRed_ps=ired
        kTextGreen_ps=igreen
        kTextBlue_ps=iblue
      endif !icolor.gt.0

      cnText_ps=sqrt(float(kTextRed_ps**2+kTextGreen_ps**2+kTextBlue_ps**2))
      if (cnText_ps.le.0.) cnText_ps=1.

      write(cline_ps,*)
     &  kTextRed_ps/cnText_ps,
     &  kTextGreen_ps/cnText_ps,
     &  kTextBlue_ps/cnText_ps,
     &  ' setrgbcolor'
      call mshplt_fill_buff(cline_ps)

      return
      end
*CMZ :  0.01/02 05/09/2014  15.41.43  by  Michael Scheer
*CMZ :  0.00/06 19/08/2014  15.02.58  by  Michael Scheer
*CMZ :  0.00/04 08/08/2014  16.34.41  by  Michael Scheer
*CMZ :  0.00/03 04/08/2014  15.42.32  by  Michael Scheer
*CMZ :  0.00/02 07/07/2014  12.22.09  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_get_marker_color(icolor,ired,igreen,iblue)

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      integer icolor,ired,igreen,iblue

      icolor=kMarkerColor_ps
      ired=kMarkerRed_ps
      igreen=kMarkerGreen_ps
      iblue=kMarkerBlue_ps

      return
      end
*CMZ :  0.01/02 05/09/2014  15.41.43  by  Michael Scheer
*CMZ :  0.00/06 19/08/2014  15.03.38  by  Michael Scheer
*CMZ :  0.00/04 08/08/2014  16.34.41  by  Michael Scheer
*CMZ :  0.00/03 04/08/2014  15.42.32  by  Michael Scheer
*CMZ :  0.00/02 07/07/2014  12.22.09  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_get_text_color(icolor,ired,igreen,iblue)

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      integer icolor,ired,igreen,iblue

      icolor=kTextColor_ps
      ired=kTextRed_ps
      igreen=kTextGreen_ps
      iblue=kTextBlue_ps

      return
      end
*CMZ :  0.01/02 05/09/2014  15.41.43  by  Michael Scheer
*CMZ :  0.01/00 23/08/2014  09.39.32  by  Michael Scheer
*CMZ :  0.00/04 11/08/2014  17.53.18  by  Michael Scheer
*CMZ :  0.00/02 07/07/2014  12.22.09  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_pause

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      character(64) cans

      write(6,*)
      write(6,*)'Hit any key to continue or type "stop" to terminate the programm!'
      read(5,'(a)')cans

      if(cans.eq.'stop'.or.cans.eq.'STOP') then
        call mshplt_end
        stop
      endif

      return
      end
*CMZ :  1.03/02 25/04/2016  12.25.01  by  Michael Scheer
*CMZ :  1.03/01 10/10/2014  13.26.55  by  Michael Scheer
*CMZ :  0.01/02 05/09/2014  15.41.43  by  Michael Scheer
*CMZ :  0.01/00 23/08/2014  10.54.42  by  Michael Scheer
*CMZ :  0.00/06 21/08/2014  16.34.08  by  Michael Scheer
*CMZ :  0.00/04 07/08/2014  11.34.50  by  Michael Scheer
*CMZ :  0.00/02 07/07/2014  12.22.09  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_draw_title

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      real x,y
      integer lentit,icolor,ired,igreen,iblue
      character(2048) gtit

      write(lun_ps,'(a)')"% begin of mshplt_draw_title"

      call mshplt_get_text_color(icolor,ired,igreen,iblue)
      call mshplt_set_text_color(icolor,ired,igreen,iblue)

      call mshplt_flush_buff

      if (idrawgtit_ps.ne.1) then
        write(lun_ps,'(a)')"% end of mshplt_draw_title"
        return
      endif

      gtit=gtit_ps
      lentit=len_trim(gtit)

      x=xleftorig_ps+xsizorig_ps/2.
     &  -gsiz_ps*lentit/5.*scaletxt_ps+offgtitx_ps*scaletxt_ps
      y=ybottomorig_ps+ysizorig_ps-offgtity_ps*scaletxt_ps

      if (ihigzmode_ps.eq.0) then
        write(lun_ps,*)'/Helvetica findfont ',gsiz_ps*scaletxt_ps,
     &    ' scalefont setfont'
      else
        write(lun_ps,*)'/Helvetica findfont ',gsiz_ps*scaletxt_ps*1.5,
     &    ' scalefont setfont'
      endif

      write(lun_ps,*)x,y,' moveto'
      write(lun_ps,*)'(',gtit(1:lentit),') show'
      write(lun_ps,*)xleftorig_ps,ybottomorig_ps,' moveto %goto origin'

      if (ihigzmode_ps.eq.0) then
        write(lun_ps,*)'/Helvetica findfont ',chhe_ps*scaletxt_ps,
     &    ' scalefont setfont'
      else
        write(lun_ps,*)'/Helvetica findfont ',chhe_ps*scaletxt_ps*1.5,
     &    ' scalefont setfont'
      endif

      write(lun_ps,'(a)')"% end of mshplt_draw_title"
      return
      end
*CMZ :  0.01/03 22/09/2014  21.17.47  by  Michael Scheer
*CMZ :  0.01/02 05/09/2014  15.41.43  by  Michael Scheer
*CMZ :  0.01/00 24/08/2014  13.06.42  by  Michael Scheer
*CMZ :  0.00/04 11/08/2014  17.53.18  by  Michael Scheer
*CMZ :  0.00/02 07/07/2014  12.22.09  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_draw_axis_titles(xtit,ytit)

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      real x,y,chhe,tsiz,tang
      character(*) xtit,ytit

      chhe=chhe_ps
      tsiz=tsiz_ps

      call mshplt_set_character_height(tsiz)

      if (ihigzmode_ps.ne.0) then
        x=wxmax_ps
     &    -len_trim(xtit)*tsiz_ps/scalex_ps*scaletxt_ps
        y=wymin_ps-ytitoff_ps/scaley_ps*scaletxt_ps
      else
        x=(wxmax_ps+wxmin_ps)/2.
     &    -len_trim(xtit)*tsiz_ps/2./scalex_ps*scaletxt_ps
        y=wymin_ps-ytitoff_ps/scaley_ps*scaletxt_ps
      endif

      call mshplt_text_raw(x,y,xtit)

      if (ihigzmode_ps.ne.0) then
        y=wymax_ps
     &    -len_trim(ytit)*tsiz_ps/scaley_ps*scaletxt_ps
        x=wxmin_ps-xtitoff_ps/scalex_ps*scaletxt_ps
      else
        y=(wymax_ps+wymin_ps)/2.
     &    -len_trim(ytit)*tsiz_ps/2./scaley_ps*scaletxt_ps
        x=wxmin_ps-xtitoff_ps/scalex_ps*scaletxt_ps
      endif

      tang=tang_ps
      tang_ps=90.

      call mshplt_text_raw(x,y,ytit)

      call mshplt_set_character_height(chhe)

      tang_ps=tang

      return
      end
*CMZ :  1.01/02 25/09/2014  12.33.10  by  Michael Scheer
*CMZ :  0.01/02 05/09/2014  15.41.43  by  Michael Scheer
*CMZ :  0.01/00 24/08/2014  12.42.28  by  Michael Scheer
*CMZ :  0.00/06 22/08/2014  15.23.52  by  Michael Scheer
*CMZ :  0.00/04 07/08/2014  11.49.49  by  Michael Scheer
*CMZ :  0.00/02 10/07/2014  15.55.54  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_set_text_angle(tang)

      implicit none

      real tang

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      tang_ps=tang

c      write(cline_ps,*)tang_ps,' rotate'
c      call mshplt_fill_buff(cline_ps)

      return
      end
*CMZ :  1.02/00 03/10/2014  10.54.54  by  Michael Scheer
*CMZ :  0.01/02 05/09/2014  16.58.57  by  Michael Scheer
*CMZ :  0.00/06 19/08/2014  15.09.19  by  Michael Scheer
*CMZ :  0.00/04 07/08/2014  15.52.16  by  Michael Scheer
*CMZ :  0.00/03 04/08/2014  15.42.32  by  Michael Scheer
*CMZ :  0.00/02 07/07/2014  12.22.09  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_set_frame_color(icolor,ired,igreen,iblue)

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      integer icolor,ired,igreen,iblue

      kFrameColor_ps=icolor

      if (icolor.gt.0) then
        if (icolor.eq.1) then
          kFrameRed_ps=0
          kFrameGreen_ps=0
          kFrameBlue_ps=0
        else if (icolor.eq.2) then
          kFrameRed_ps=1
          kFrameGreen_ps=0
          kFrameBlue_ps=0
        else if (icolor.eq.3) then
          kFrameRed_ps=0
          kFrameGreen_ps=1
          kFrameBlue_ps=0
        else if (icolor.eq.4) then
          kFrameRed_ps=0
          kFrameGreen_ps=0
          kFrameBlue_ps=1
        else if (icolor.eq.5) then
          kFrameRed_ps=1
          kFrameGreen_ps=1
          kFrameBlue_ps=0
        else if (icolor.eq.6) then
          kFrameRed_ps=1
          kFrameGreen_ps=0
          kFrameBlue_ps=1
        else if (icolor.eq.7) then
          kFrameRed_ps=0
          kFrameGreen_ps=1
          kFrameBlue_ps=1
        else if (icolor.eq.8) then
          kFrameRed_ps=35
          kFrameGreen_ps=85
          kFrameBlue_ps=33
        endif
      else
        kFrameColor_ps=-1
        kFrameRed_ps=ired
        kFrameGreen_ps=igreen
        kFrameBlue_ps=iblue
      endif !icolor.gt.0

      cnFrame_ps=sqrt(float(kFrameRed_ps**2+kFrameGreen_ps**2+kFrameBlue_ps**2))
      if (cnFrame_ps.le.0) cnFrame_ps=1.

      write(cline_ps,*)
     &  kFrameRed_ps/cnFrame_ps,
     &  kFrameGreen_ps/cnFrame_ps,
     &  kFrameBlue_ps/cnFrame_ps,
     &  ' setrgbcolor'
      call mshplt_fill_buff(cline_ps)

      call mshplt_flush_buff

      return
      end
*CMZ :  0.01/03 22/09/2014  20.54.58  by  Michael Scheer
*CMZ :  0.01/02 05/09/2014  15.41.43  by  Michael Scheer
*CMZ :  0.00/04 07/08/2014  15.19.17  by  Michael Scheer
*CMZ :  0.00/02 07/07/2014  12.22.09  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_box(x1,y1,x2,y2)

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      real x1,y1,x2,y2

      call mshplt_line_raw(x1,y1,x2,y1)
      call mshplt_line_raw(x2,y1,x2,y2)
      call mshplt_line_raw(x2,y2,x1,y2)
      call mshplt_line_raw(x1,y2,x1,y1)

      return
      end
*CMZ :  0.01/02 05/09/2014  15.41.43  by  Michael Scheer
*CMZ :  0.00/04 11/08/2014  17.53.18  by  Michael Scheer
*CMZ :  0.00/02 07/07/2014  12.22.09  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_set_logx(iset)

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      integer iset

      log10x_ps=iset

      return
      end
*CMZ :  0.01/02 05/09/2014  15.41.43  by  Michael Scheer
*CMZ :  0.00/04 11/08/2014  17.53.18  by  Michael Scheer
*CMZ :  0.00/02 07/07/2014  12.22.09  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_set_logy(iset)

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      integer iset

      log10y_ps=iset

      return
      end
*CMZ :  0.01/02 05/09/2014  15.41.43  by  Michael Scheer
*CMZ :  0.00/04 11/08/2014  17.53.18  by  Michael Scheer
*CMZ :  0.00/02 07/07/2014  12.22.09  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_set_logz(iset)

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      integer iset

      log10z_ps=iset

      return
      end
*CMZ :  0.01/02 05/09/2014  15.41.43  by  Michael Scheer
*CMZ :  0.00/04 11/08/2014  17.53.18  by  Michael Scheer
*CMZ :  0.00/02 07/07/2014  12.22.09  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_get_logx(iget)

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      integer iget

      iget=log10x_ps

      return
      end
*CMZ :  0.01/02 05/09/2014  15.41.43  by  Michael Scheer
*CMZ :  0.00/04 11/08/2014  17.53.18  by  Michael Scheer
*CMZ :  0.00/02 07/07/2014  12.22.09  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_get_logy(iget)

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      integer iget

      iget=log10y_ps

      return
      end
*CMZ :  0.01/02 05/09/2014  15.41.43  by  Michael Scheer
*CMZ :  0.00/04 11/08/2014  17.53.18  by  Michael Scheer
*CMZ :  0.00/02 07/07/2014  12.22.09  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_get_logz(iget)

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      integer iget

      iget=log10z_ps

      return
      end
*CMZ :          02/08/2018  18.13.30  by  Michael Scheer
*CMZ :  1.03/02 23/09/2016  13.46.02  by  Michael Scheer
*CMZ :  1.02/00 01/10/2014  14.43.01  by  Michael Scheer
*CMZ :  1.01/02 27/09/2014  14.55.40  by  Michael Scheer
*CMZ :  1.01/00 24/09/2014  14.41.50  by  Michael Scheer
*CMZ :  1.00/01 24/09/2014  11.26.19  by  Michael Scheer
*CMZ :  0.01/03 23/09/2014  14.12.14  by  Michael Scheer
*CMZ :  0.01/02 22/09/2014  19.04.48  by  Michael Scheer
*CMZ :  0.01/01 27/08/2014  08.54.49  by  Michael Scheer
*CMZ :  0.00/06 22/08/2014  15.47.40  by  Michael Scheer
*CMZ :  0.00/05 18/08/2014  13.15.26  by  Michael Scheer
*CMZ :  0.00/04 14/08/2014  12.58.41  by  Michael Scheer
*CMZ :  0.00/03 05/08/2014  15.57.52  by  Michael Scheer
*CMZ :  0.00/02 09/07/2014  14.42.14  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_axis_taylor(
     &  xmin,xmax,ymin,ymax,smin,smax,
     &  nlab,chlab,rlabsiz,xlabrel,ylaboff,xoffexp,yoffexp,rlabangrel,
     &  ntic,ticsiz,ticposrel,ticangrel,
     &  titsiz,titposrel,titoff,titangrel,chtit
     &  )

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      integer :: ical=0

      integer nlab,ntic,i,nbig,ic,ifirst,nexp,iprimitiv,k,nskip,iexpmode,
     &  nexpmax,kbig(100),ndig,nticmax,inv

      real xmin,xmax,ymin,ymax,axang,x(2),y(2),ango,chheo,
     &  rlabsiz(*),xlabrel(*),ylaboff(*),xoffexp,yoffexp,rlabangrel(*),axlen,
     &  titsiz,titangrel,ticsiz(*),ticposrel(*),ticangrel(*),cosang,sinang,
     &  titposrel,titoff,xlab(100),smin,smax,rmanti,rtotlabx,rtotlaby,
     &  costitang,sintitang,costicang,sinticang,coslabang,sinlabang
     &  ,ticsizo,dxabs,xlabinv(100)

      character(*) chtit,chlab(*)
      character(12) chreal
      character(3) chexp

      ical=ical+1

      ndig=100

      if (ymax.eq.ymin.and.xmax.eq.xmin) then
        print*,'*** Error in mshplt_axis_taylor: Zero length axis ***'
        return
      endif

      call mshplt_flush_buff
      write(lun_ps,'(a)')'% begin of mshplt_axis_taylor'

      ticsizo=ticsiz(1)

      ango=ang_ps
      chheo=chhe_ps
      axang=atan2((ymax-ymin)*scaley_ps,(xmax-xmin)*scalex_ps) !radian
      cosang=cos(axang)
      sinang=sin(axang)
      axang=axang/1.745329251994e-2 !degree, on paper
      axlen=sqrt(((ymax-ymin)*scaley_ps)**2+((xmax-xmin)*scalex_ps)**2) !paper

      ifirst=1
      kbig=0
      iprimitiv=0

      call mshplt_line_raw(xmin,ymin,xmax,ymax)

      if (nlab.lt.0) then

c nlab.lt.0 means automatic mode, where abs(nlab) is the garantied size of the
c considered arrays. Nlab is number of labels and ntic is number of all tics
c Most values are overwritten and returned, so they must not be constants.

        if (len(chlab(1)).lt.12) then
          write(lun_ps,*)
     &      '%*** Error in mshplt_axis_taylor: For automatic mode, the size of'
          write(lun_ps,*)
     &      '%*** must be at least 12, i.e. character(12)'
          write(6,*)
     &      '*** Error in mshplt_axis_taylor: For automatic mode, the size of'
          write(6,*)
     &      '*** must be at least 12, i.e. character(12)'
          goto 9999
        endif

        if (-nlab.gt.100) then
          write(6,*)
     &      '*** Error in mshplt_axis_taylor: Dimension of labels exceeded'
          write(lun_ps,*)
     &      '%*** Error in mshplt_axis_taylor: Dimension of labels exceeded'
          goto 9999
        endif

        nticmax=max(3,int(axlen))

c111     continue
        iprimitiv=1
c Problems with negative axis
        if (smax.lt.0.0d0) then
          call mshplt_axis_division(abs(smax),abs(smin),nticmax,nlab,xlabinv)
          do inv=1,nlab
            xlab(inv)=-xlabinv(nlab-inv+1)
            write(chlab(1),'(2pe12.3)')xlab(inv)
            read(chlab(1),*)xlab(inv)
          enddo
        else
          call mshplt_axis_division(smin,smax,nticmax,nlab,xlab)
          do inv=1,nlab
            write(chlab(1),'(2pe12.3)')xlab(inv)
            read(chlab(1),*)xlab(inv)
          enddo
        endif

        dxabs=abs(xlab(2)-xlab(1))
        iexpmode=0
        write(chlab(1),'(2pe12.3)')dxabs
        read(chlab(1),*)dxabs
        write(chlab(1),'(1pe12.5)')dxabs
        read(chlab(1)(10:12),'(i3)')nexpmax
        if (abs(xlab(nlab)-xlab(1)).gt.9999.
     &    .or.
     &    abs(xlab(nlab)-xlab(1)).lt.0.001) iexpmode=1

        if (iexpmode.eq.0) nexpmax=0

        ntic=nlab
        nbig=0

        if (nlab.eq.1) then
          ticposrel(1)=0.5
          xlabrel(1)=0.5
        else if (nlab.eq.2) then
          ticposrel(1)=0.25
          xlabrel(1)=0.25
          ticposrel(2)=0.75
          xlabrel(2)=0.75
        endif

        do i=1,ntic
          ticsiz(i)=ticsizo
          ticangrel(i)=ticangrel(1)
          rlabangrel(i)=rlabangrel(1)
          rlabsiz(i)=rlabsiz(1)
          if (iexpmode.eq.0) then
            write(chlab(i),'(g12.5)')xlab(i)
          else
            write(chlab(i),'(1pe12.5)')xlab(i)
          endif
        enddo

        if (nlab.ge.2) then

          do i=1,ntic
            ticposrel(i)=(xlab(i)-smin)/(smax-smin)
            xlabrel(i)=ticposrel(i)
          enddo

          write(chreal,'(2pe12.3)')xlab(2)-xlab(1)

          if (chreal(2:3).eq.'10'
     &        .or.
     &        chreal(2:3).eq.'25'
     &        .or.
     &        chreal(2:3).eq.'50'
     &        ) then

            iprimitiv=0

            do i=1,ntic
              chlab(i)=''
              ticsiz(i)=ticsizo
              ticangrel(i)=90.
              ticposrel(i)=(xlab(i)-smin)/(smax-smin)
              rlabangrel(i)=rlabangrel(1)
              rlabsiz(i)=rlabsiz(1)
              xlabrel(i)=ticposrel(i)
c              write(chreal,'(g12.5)')xlab(i)/10**nexpmax
            enddo !nlab

          endif
        endif !nlab > 2

c        if (iprimitiv.ne.0) then

c          write(chlab(i),'(g12.3)')xlab(i)

c        else

          do i=1,ntic

            if (xlab(i).eq.0.0) then
              chlab(i)='0.0'
              cycle
            endif

            write(chlab(i),'(f12.4)')xlab(i)/10.**nexpmax

            read(chlab(i)(1:8),*)rmanti
            nexp=0

            ifirst=0
            do ic=1,len_trim(chlab(i))
              if (chlab(i)(ic:ic).ne.' '.and.ifirst.eq.0) then
                ifirst=ic
                goto 123
              endif
            enddo
123         continue

            chlab(i)=chlab(i)(ifirst:len_trim(chlab(i)))

            do ic=len_trim(chlab(i)),1,-1
              if (chlab(i)(ic:ic).eq.'0') then
                chlab(i)(ic:ic)=''
              else
                goto 124
              endif
            enddo

124         continue

            if(chlab(i)(len_trim(chlab(i)):len_trim(chlab(i))).eq.'.')
     &        chlab(i)(len_trim(chlab(i)):len_trim(chlab(i)))=''

          enddo ! labels

c        endif !iprimitiv
      endif !nlab.lt.0, i.e. automatic mode

      do i=1,ntic

        costicang=cos((ticangrel(i)+axang)*1.745329251994e-2)
        sinticang=sin((ticangrel(i)+axang)*1.745329251994e-2)

        x(1)=xmin+ticposrel(i)*(xmax-xmin)
        y(1)=ymin+ticposrel(i)*(ymax-ymin)
        call mshplt_view_to_world(x(1),y(1),x(1),y(1))

        x(2)=x(1)+ticsiz(i)*costicang
        y(2)=y(1)+ticsiz(i)*sinticang

        call mshplt_world_to_view(x(1),y(1),x(1),y(1))
        call mshplt_world_to_view(x(2),y(2),x(2),y(2))
c20160923        if (x(1).le.xmax) then
          call mshplt_pline_raw(2,x,y)
c20160923        endif
      enddo

      rtotlabx=0.
      do i=1,ntic
        if (len_trim(chlab(i)).ne.0) then
          coslabang=cos((rlabangrel(i)+axang)*1.745329251994e-2)
          sinlabang=sin((rlabangrel(i)+axang)*1.745329251994e-2)
          rtotlabx=rtotlabx+len_trim(chlab(i))*rlabsiz(i)/2.*coslabang
          rtotlaby=rtotlaby+len_trim(chlab(i))*rlabsiz(i)/2.*sinlabang
        endif
      enddo

      if (sqrt(rtotlabx**2+rtotlaby**2).gt.axlen*0.5) then
c        nticmax=nticmax-1
c        goto 111
        if (nbig.ge.2) then
          do k=1,ntic
            if (kbig(k).eq.0) then
              chlab(k)=''
            endif
          enddo
        else !nbig
          if (ntic.gt.2.and.mod(ntic,2).eq.0) then
            nskip=int(sqrt(rtotlabx**2+rtotlaby**2)/axlen+0.5)*2
            if (ntic.gt.4) then
              do i=1,ntic,nskip+1
                do k=1,nskip
                  chlab(i-1+k)=''
                enddo
              enddo
            else
              chlab(2)=''
              chlab(4)=''
            endif
          else
            do i=2,ntic-1
              if (i.ne.ntic/2+1) chlab(i)=''
            enddo
          endif
        endif !nbig
      endif !to many letters

      do i=1,nlab
        if (rlabsiz(i).gt.0.) then

          coslabang=cos((rlabangrel(i)+axang)*1.745329251994e-2)
          sinlabang=sin((rlabangrel(i)+axang)*1.745329251994e-2)

          call mshplt_set_text_angle(axang+rlabangrel(i))
          call mshplt_set_character_height(rlabsiz(i))

          x(1)=xmin+xlabrel(i)*(xmax-xmin)
          y(1)=ymin+xlabrel(i)*(ymax-ymin)

          call mshplt_view_to_world(x(1),y(1),x(1),y(1))

          x(1)=x(1)+ylaboff(i)*sinang
     &      -len_trim(chlab(i))*rlabsiz(i)/4.*coslabang

          y(1)=y(1)-ylaboff(i)*cosang
     &      -sign(1.,90.-theta_ps)*rlabsiz(i)/2.*sinang-len_trim(chlab(i))*rlabsiz(i)/4.*sinlabang

          call mshplt_world_to_view(x(1),y(1),x(1),y(1))
          call mshplt_text_raw(x(1),y(1),chlab(i))

c20170602??          if (i.eq.nlab.and.iprimitiv.eq.0.and.iexpmode.ne.0) then
          if (i.eq.nlab.and.iexpmode.ne.0) then

            x(1)=xmax
            y(1)=ymax

            call mshplt_view_to_world(x(1),y(1),x(1),y(1))

            if (ylaboff(i).gt.0.) then
              x(1)=x(1)+xoffexp*cosang+yoffexp*sinang
              y(1)=y(1)+xoffexp*sinang-yoffexp*cosang
            else
              x(1)=x(1)+xoffexp*cosang-(yoffexp+3.*rlabsiz(i))*sinang
              y(1)=y(1)+xoffexp*sinang+yoffexp*cosang
            endif

            call mshplt_world_to_view(x(1),y(1),x(1),y(1))
            call mshplt_set_character_height(rlabsiz(i)*1.2)
            call mshplt_text_raw(x(1),y(1),'x 10')

            write(chexp,'(i3)')nexpmax

            if (chexp(1:2).eq.'  ') then
              call mshplt_text_upper(chexp(3:3))
            else if (chexp(1:1).eq.' ') then
              call mshplt_text_upper(chexp(2:3))
            else
              call mshplt_text_upper(chexp(1:3))
            endif

          endif
        endif
      enddo

c axis-title

      if (titsiz.gt.0.0) then

        costitang=cos((titangrel+axang)*1.745329251994e-2)
        sintitang=sin((titangrel+axang)*1.745329251994e-2)
        call mshplt_set_text_angle(titangrel+axang)
        call mshplt_set_character_height(titsiz)
        x(1)=xmin+titposrel*(xmax-xmin)
        y(1)=ymin+titposrel*(ymax-ymin)
        call mshplt_view_to_world(x(1),y(1),x(1),y(1))
        x(1)=x(1)+titoff*sintitang
     &  -len_trim(chtit)*titsiz/4.*costitang
        y(1)=y(1)-titoff*costitang
     &    -len_trim(chtit)*titsiz/4.*sintitang
        call mshplt_world_to_view(x(1),y(1),x(1),y(1))
        call mshplt_text_raw(x(1),y(1),chtit)

      endif

9999  continue

      call mshplt_set_character_height(chheo)
      call mshplt_set_text_angle(ango)
      call mshplt_flush_buff

      write(cline_ps,*)'% end of mshplt_axis_taylor'
      write(lun_ps,'(a)')cline_ps(2:len_trim(cline_ps))

      return
      end
*CMZ :  0.01/02 11/09/2014  17.22.55  by  Michael Scheer
*CMZ :  0.00/06 22/08/2014  15.31.43  by  Michael Scheer
*CMZ :  0.00/04 07/08/2014  08.43.21  by  Michael Scheer
*CMZ :  0.00/03 05/08/2014  15.53.35  by  Michael Scheer
*CMZ :  0.00/02 08/07/2014  11.10.50  by  Michael Scheer
*CMZ :  0.00/01 07/07/2014  12.00.23  by  Michael Scheer
*CMZ : 00.00/02 04/07/2014  21.11.57  by  Michael Scheer
*-- Author :    Michael Scheer   11/06/2014
      subroutine mshplt_world_log_axis(
     &  wxmin,wxmax,wymin,wymax,smin,smax,title,
     &  iticside,ilabside,titang,offtit,offlab
     &  )

c Draw an axis from (wxmin,wxmax) to (wymin,wymax) on pad
c with labeling from smin to smax, and title title

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      real wxmin,wxmax,wymin,wymax,smin,smax,xlength,ang,chhe,tiche,
     &  titang,offtit,offlab
      integer iticside,ilabside

      character(*) title

      write(lun_ps,*)'% begin of mshplt_world_log_axis'

      ilabmod_ps=ilabside

      xlength=sqrt((wxmax-wxmin)**2+(wymax-wymin)**2)

      if (xlength.eq.0.0.or.smin.eq.smax) then
        print*,'*** Error in mshplt_world_log_axis: Zero-length axis!'
        return
      endif

      ang=atan2(wymax-wymin,wxmax-wxmin)*57.29578

      if (ihigzmode_ps.eq.0) then
        chhe=chhe_ps*scaletxt_ps
        tiche=ticsiz_ps
      else
        chhe=chhe_ps*scaletxt_ps*1.5
        tiche=ticsiz_ps
      endif

      call mshplt_log_axis(wxmin,wymin,xlength,smin,smax,ang,chhe,
     &  tiche,title,iticside,ilabside,titang,offtit,offlab)

      write(lun_ps,*)'% end of mshplt_world_log_axis'

      return
      end
*CMZ :          20/10/2023  20.40.14  by  Michael Scheer
*CMZ :  1.03/01 08/10/2014  12.21.58  by  Michael Scheer
*CMZ :  1.01/02 26/09/2014  11.35.12  by  Michael Scheer
*CMZ :  1.01/00 24/09/2014  16.31.07  by  Michael Scheer
*CMZ :  1.00/01 23/09/2014  20.14.56  by  Michael Scheer
*CMZ :  0.01/03 22/09/2014  21.17.48  by  Michael Scheer
*CMZ :  0.01/02 22/09/2014  19.54.46  by  Michael Scheer
*-- Author :    Michael Scheer   08/09/2014
      subroutine mshplt_log_axis(
     &  xorig, yorig, xylength, smin, smax, ang,
     &  chhe, ticheight, title,iticside,ilabside,titang,offtit,offlab)

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      real xorig, yorig, xylength, smin, smax, ang,chhe,ticheight
      real titang,offtit,offlab,sminl,smaxl,
     &  x(2),y(2),aLabSide,xdec,ydec,xv,yv,nexp,ango,chheo

      integer iticside,ilabside,j,l,idecmin,idecmax
      character(*) title

      integer ntic,nlab,ndec,i
      real xmin,xmax,ymin,ymax,sinang,cosang,dxy,tside,costit,sintit

      character(8) c8

      call mshplt_flush_buff

      write(lun_ps,'(a)')'% begin of mshplt_log_axis'

      if (iticside.lt.0) then
        tside=1.
      else if (iticside.gt.0) then
        tside=-1.
      else
        tside=0.
      endif

      if (ilabside.lt.0) then
        aLabSide=-1.
      else if (ilabside.gt.0) then
        aLabSide=0.25
      else
        aLabSide=0.
      endif

      if (smin.le.0.0) then
        sminl=-30.
      else
        sminl=alog10(smin)
      endif

      if (smax.le.0.0) then
        smaxl=30.
      else
        smaxl=alog10(smax)
      endif

      cosang=cos(ang*1.745329e-2)
      sinang=sin(ang*1.745329e-2)

      call  mshplt_world_to_view(xorig,yorig,xmin,ymin)

      call  mshplt_world_to_view(
     &  xorig+cosang*xylength,yorig+sinang*xylength,
     &  xmax,ymax)

      call mshplt_line_raw(xmin,ymin,xmax,ymax)

      idecmin=nint(sminl)
      if (sminl.ne.idecmin) then
        print*,'*** Warning in mshplt_log_axis: Routine must be called with integer power of 10'
      endif
      idecmax=nint(smaxl)
      if (smaxl.ne.idecmax) then
        print*,'*** Warning in mshplt_log_axis: Routine must be called with integer power of 10'
      endif

      if (idecmax.lt.30.and.10.0**idecmax.lt.smax.or.idecmax.eq.idecmin) then
        idecmax=idecmax+1
      endif

      ndec=idecmax-idecmin

      nlab=ndec*9
      ntic=nlab

      dxy=xylength/ndec

      write(c8,'(i8)')max(abs(idecmin),abs(idecmax))
      do l=1,8
        if(c8(l:l).ne.' ') then
          nexp=8-l+1
          goto 1
        endif
      enddo

 1    continue

      do i=1,ndec

        if (i.eq.1) then

          xdec=xorig
          ydec=yorig

          x(1)=xdec
          x(2)=x(1)+ticheight*tside*sinang
          y(1)=ydec
          y(2)=y(1)+ticheight*tside*cosang

          call  mshplt_world_to_view(x(1),y(1),x(1),y(1))
          call  mshplt_world_to_view(x(2),y(2),x(2),y(2))

          if (ilabside.ne.0.and.aLabHeight_ps.gt.0.0) then
            x(2)=xdec+aLabSide*sinang*(offlab+(1.+nexp)*aLabHeight_ps/2.)
     &        -aLabHeight_ps*cosang
            y(2)=ydec+aLabSide*offlab*cosang-aLabHeight_ps
     &        +aLabHeight_ps*cosang
            call  mshplt_world_to_view(x(2),y(2),xv,yv)
            call mshplt_set_character_height(aLabHeight_ps*1.2)
            call mshplt_text_raw(xv,yv,"10")

            write(c8,'(i8)')idecmin
            do l=1,8
              if(c8(l:l).ne.' ') goto 12
            enddo
 12         continue

            call mshplt_text_upper(c8(l:8))

          endif
        endif

        xdec=xorig+i*dxy*cosang
        ydec=yorig+i*dxy*sinang

        x(1)=xdec
        x(2)=x(1)+ticheight*tside*sinang
        y(1)=ydec
        y(2)=y(1)+ticheight*tside*cosang

        call  mshplt_world_to_view(x(1),y(1),x(1),y(1))
        call  mshplt_world_to_view(x(2),y(2),x(2),y(2))

        call mshplt_pline_raw(2,x,y)

        if (ilabside.ne.0.and.aLabHeight_ps.gt.0.0) then
          x(2)=xdec+aLabSide*sinang*(offlab+(1.+nexp)*aLabHeight_ps/2.)
     &      -aLabHeight_ps*cosang
          y(2)=ydec+aLabSide*offlab*cosang-aLabHeight_ps
     &      +aLabHeight_ps*cosang
          call  mshplt_world_to_view(x(2),y(2),xv,yv)
          call mshplt_set_character_height(aLabHeight_ps*1.2)
          call mshplt_text_raw(xv,yv,"10")

          write(c8,'(i8)')idecmin+i
          do l=1,8
            if(c8(l:l).ne.' ') goto 11
          enddo
 11       continue

          call mshplt_text_upper(c8(l:8))

        endif

        do j=2,9
          x(1)=xorig+((i-1)*dxy+alog10(float(j))*dxy)*cosang
          x(2)=x(1)+ticheight*sinang*tside*0.8
          y(1)=yorig+((i-1)*dxy+alog10(float(j))*dxy)*sinang
          y(2)=y(1)+ticheight*cosang*tside*0.8
          call  mshplt_world_to_view(x(1),y(1),x(1),y(1))
          call  mshplt_world_to_view(x(2),y(2),x(2),y(2))
          call mshplt_pline_raw(2,x,y)
        enddo
      enddo

      chheo=chhe_ps
      ango=ang_ps

      costit=cos((ang+titang)*1.745329e-2)
      sintit=sin((ang+titang)*1.745329e-2)

      call mshplt_set_text_angle(ang+titang)

      if (aLabSide.gt.0.0) then
        x(1)=xorig
     &    +0.5*xylength*costit
     &    -len_trim(title)*chhe/4.*costit
     &    +(6.+nexp)/2.*aLabSide*(offtit+aLabHeight_ps)*sintit
        y(1)=yorig
     &    +0.5*xylength*sintit
     &    -len_trim(title)*chhe/4.*sintit
     &    +(6.+nexp)/2.*aLabSide*(offtit+aLabHeight_ps)*costit
      else if (aLabSide.lt.0.0) then
        x(1)=xorig
     &    +0.5*xylength*costit
     &    +aLabSide*len_trim(title)*chhe/4.*costit
     &    +aLabSide*(offtit+aLabHeight_ps)*sintit
        y(1)=yorig
     &    +0.5*xylength*sintit
     &    +aLabSide*len_trim(title)*chhe/4.*sintit
     &    +aLabSide*(offtit+aLabHeight_ps)*costit
      else
        goto 9999
      endif

      call mshplt_world_to_view(x(1),y(1),xv,yv)
      call mshplt_set_character_height(chhe)
      call mshplt_text_raw(xv,yv,title(1:len_trim(title)))

 9999 continue

      call mshplt_set_character_height(chheo)
      call mshplt_set_text_angle(ango)

      call mshplt_flush_buff
      write(lun_ps,'(a)')'% end of mshplt_log_axis'

      return
      end
*CMZ :  0.01/02 08/09/2014  17.36.39  by  Michael Scheer
*-- Author :    Michael Scheer   08/09/2014
      subroutine mshplt_number(
     &  x,     ! x coordinate of string ( user units )
     &  y,     ! y coordinate of string ( user units )
     &  chhe,     ! height of string ( user units )
     &  rnum,     ! number to be pplotted
     &  ang,      ! angle in respect to x-axis
     &  ndig)  ! number of digits right of decimal point

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      real rnum,x,y,chhe,ang
      integer ndig,ifirst

      character(16) chnum
      character(8) chform

      chform(1:4)='(f10.'

      if (ndig.ge.0) then
        chform(1:4)='(f10.'
        write(chform(6:6),'(i1)')ndig
        chform(7:7)=')'
        write(chnum,chform)rnum
      else
        chform(1:5)='(I10)'
        write(chnum,chform)ifix(rnum+0.5)
      endif

      do ifirst=1,16
        if (chnum(ifirst:ifirst).ne.' ') goto 9
      enddo

9     continue

      if (ihigzmode_ps.eq.0) then
        write(lun_ps,*)'/Helvetica findfont',chhe*scaletxt_ps,
     &    ' scalefont setfont'
      else
        write(lun_ps,*)'/Helvetica findfont',chhe*scaletxt_ps*1.3,
     &    ' scalefont setfont'
      endif

      write(lun_ps,*)
     &  x,y,' moveto ',
     &  ang,' rotate ',
     &  '(',chnum(ifirst:len_trim(chnum)),') show',
     &  -ang,' rotate'

      return
      end
*CMZ :  0.01/02 11/09/2014  11.44.51  by  Michael Scheer
*CMZ :  0.00/04 14/08/2014  12.58.41  by  Michael Scheer
*CMZ :  0.00/03 06/08/2014  09.32.19  by  Michael Scheer
*CMZ :  0.00/02 09/07/2014  14.48.35  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_view_to_world(xv,yv,xw,yw)

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      real xv,yv,xw,yw,x,y

      x=xv
      y=yv

      xw=xleft_ps+scalex_ps*(x-wxmin_ps)
      yw=ybottom_ps+scaley_ps*(y-wymin_ps)

      return
      end
*CMZ :  0.01/02 11/09/2014  11.53.53  by  Michael Scheer
*CMZ :  0.00/04 14/08/2014  12.58.41  by  Michael Scheer
*CMZ :  0.00/03 06/08/2014  09.32.19  by  Michael Scheer
*CMZ :  0.00/02 09/07/2014  14.48.35  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_world_to_view(xw,yw,xv,yv)

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      real xv,yv,xw,yw,x,y

      x=xw
      y=yw

      xv=wxmin_ps+(x-xleft_ps)/scalex_ps
      yv=wymin_ps+(y-ybottom_ps)/scaley_ps

      return
      end
*CMZ :          31/07/2018  12.07.04  by  Michael Scheer
*CMZ :  1.02/00 03/10/2014  12.49.41  by  Michael Scheer
*CMZ :  1.01/02 25/09/2014  09.32.57  by  Michael Scheer
*CMZ :  1.01/00 24/09/2014  14.33.11  by  Michael Scheer
*CMZ :  1.00/01 24/09/2014  12.18.04  by  Michael Scheer
*CMZ :  0.01/02 18/09/2014  13.01.02  by  Michael Scheer
*CMZ :  0.01/01 27/08/2014  15.30.42  by  Michael Scheer
*CMZ :  0.01/00 25/08/2014  10.34.58  by  Michael Scheer
*CMZ :  0.00/04 14/08/2014  13.19.56  by  Michael Scheer
*CMZ :  0.00/03 06/08/2014  09.00.02  by  Michael Scheer
*CMZ :  0.00/02 09/07/2014  13.35.02  by  Michael Scheer
*CMZ :  0.00/01 07/07/2014  09.10.33  by  Michael Scheer
*CMZ : 00.00/02 30/06/2014  10.35.42  by  Michael Scheer
*-- Author :    Michael Scheer   27/06/2014
      subroutine mshplt_stop

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEEP,mplot.

      real rescale_mshplt

      character(2048)
     &  fileeps_mshplt,
     &  viewer_mshplt,viewer_kill_mshplt

      common/mplotc/
     &  rescale_mshplt,
     &  fileeps_mshplt,
     &  viewer_mshplt,viewer_kill_mshplt
*KEND.

      integer nitemp
      parameter(nitemp=100)

      integer ifile,lun,lunps,ipos(2,nitemp),nwords,istat,leneps
      real scalex,scaley

      character(2048) fileps

      call mshplt_flush_buff
      write(lun_ps,*)'% begin of mshplt_end'
      write(lun_ps,*)'showpage'

c      lunps=lunall_ps(nfile_ps)+1
c      call util_get_free_lun(lunps)

      do ifile=1,nfile_ps
        if (rescale_mshplt.ne.1..and.rescale_mshplt.gt.0.0) then
          leneps=len_trim(fileall_ps(ifile))
          fileps=fileall_ps(ifile)(1:leneps-4)//
     &      '_scaled'//fileall_ps(ifile)(leneps-3:leneps)
          open(newunit=lunps,file=fileps(1:len_trim(fileps)),status='unknown')
          lun=lunall_ps(ifile)
          if (lun.gt.0) then
            write(lunall_ps(ifile),*)'showpage'
            rewind(lun)
1           read(lun,'(a)',end=9)cline_ps
            call util_string_split(cline_ps,nitemp,nwords,ipos,istat)
            if (cline_ps(ipos(1,1):ipos(2,1)).eq.'%%Title') then
              write(lunps,'(a)')'%% Title '//fileps(1:len_trim(fileps))
              goto 1
            endif
            if (nwords.ge.3) then
              if (cline_ps(ipos(1,3):ipos(2,3)).eq.'scale') then
                read(cline_ps,*)scalex,scaley
                write(cline_ps,*)scalex*rescale_mshplt,scaley*rescale_mshplt,
     &            ' scale'
                write(lunps,'(a)')cline_ps(1:len_trim(cline_ps))
                goto 11
              endif
            endif
            write(lunps,'(a)')cline_ps(1:len_trim(cline_ps))
            goto 1
11          read(lun,'(a)',end=9)cline_ps
            write(lunps,'(a)')cline_ps(1:len_trim(cline_ps))
            goto 11
9           close(lunps)
c            close(lunall_ps(ifile))
            lunall_ps(ifile)=-lunall_ps(ifile)
          endif
        else
          close(lunall_ps(ifile))
          lunall_ps(ifile)=-lunall_ps(ifile)
        endif
      enddo

      call sleep(1)

      if(irunviewer_ps.ne.0) then

        if(iviewinter_ps.ne.0) then
          write(6,*)'Hit any key to terminate viewer:'
          read(5,'(a)')
          call system(viewer_kill_mshplt(1:len_trim(viewer_kill_mshplt)))
        endif

      endif

      stop
      end
*CMZ :  0.01/02 11/09/2014  14.09.59  by  Michael Scheer
*CMZ :  0.00/06 22/08/2014  15.23.52  by  Michael Scheer
*CMZ :  0.00/04 07/08/2014  11.49.49  by  Michael Scheer
*CMZ :  0.00/02 10/07/2014  15.55.54  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_set_label_size(siz)

      implicit none

      real siz

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      aLabHeight_ps=siz

      call mshplt_fill_buff('% begin of mshplt_set_label_size')

      if (ihigzmode_ps.eq.0) then
        write(cline_ps,*)'/Helvetica findfont',aLabHeight_ps*scaletxt_ps,
     &    ' scalefont setfont'
      else
        write(cline_ps,*)'/Helvetica findfont',aLabHeight_ps*scaletxt_ps*1.3,
     &    ' scalefont setfont'
      endif

      call mshplt_fill_buff(cline_ps)

      call mshplt_fill_buff('% end of mshplt_set_label_size')

      return
      end
*CMZ :          30/01/2018  14.43.10  by  Michael Scheer
*CMZ :  1.03/01 07/10/2014  15.45.02  by  Michael Scheer
*CMZ :  1.02/00 01/10/2014  08.53.27  by  Michael Scheer
*CMZ :  1.01/02 27/09/2014  12.27.11  by  Michael Scheer
*CMZ :  0.01/02 19/09/2014  17.39.37  by  Michael Scheer
*-- Author :    Michael Scheer   15/09/2014
      subroutine mshplt_axis_division(xmin,xmax,nlabmax,nlab,xlab)

      implicit none

      real xmin,xmax,dx,xhigh,xdiff,rmandiff,xlab(100),rmanmax,x,eps
      integer il,nlab,nexpdiff,nexpmax,nlabmax

      data eps/1.0e-5/

      character(12) chreal

      if (xmax.ne.xmax.or.xmin.ne.xmin) then
        print*,"*** Error in mshplt_axis_division: Xmin or Xmax not a number ***"
        nlab=0
        return
      endif

      if (nlabmax.eq.0) then
        nlab=0
        return
      else if (nlabmax.eq.1) then
        nlab=1
        xlab(1)=(xmax+xmin)/2.
        return
      endif

      if (xmax-xmin.eq.0.0.or.
     &    abs(xmax-xmin).le.1.0e-5*(abs(xmax)+abs(xmin))
     &    ) then
        nlab=2
        xlab(1)=xmin
        xlab(2)=xmax
        return
      endif

      xdiff=max(abs(xmax),abs(xmin))
      xdiff=abs(xmax-xmin)

      write(chreal,'(1pe12.5)')xdiff
      read(chreal(10:12),'(i3)')nexpdiff
      read(chreal(1:8),*)rmandiff

      il=int(rmandiff)

      if (il.eq.1) then
        dx=0.1
        if (int(rmandiff/dx).gt.nlabmax) dx=0.25
      else if (il.eq.2) then
        dx=0.2
        if (int(rmandiff/dx).gt.nlabmax) dx=0.5
      else if (il.eq.3) then
        dx=0.5
        if (int(rmandiff/dx).gt.nlabmax) dx=1.
      else if (il.eq.4) then
        dx=0.5
        if (int(rmandiff/dx).gt.nlabmax) dx=1.
      else if (il.eq.5) then
        dx=1.
        if (int(rmandiff/dx).gt.nlabmax) dx=2.
      else if (il.eq.6) then
        dx=1.
        if (int(rmandiff/dx).gt.nlabmax) dx=2.
      else if (il.eq.7) then
        dx=1.
        if (int(rmandiff/dx).gt.nlabmax) dx=2.5
      else if (il.eq.8) then
        dx=1.
        if (int(rmandiff/dx).gt.nlabmax) dx=2.
      else if (il.eq.9) then
        dx=1.
        if (int(rmandiff/dx).gt.nlabmax) dx=2.5
      else
        dx=0.5
        if (int(rmandiff/dx).gt.nlabmax) dx=1.
      endif

      do while (int(rmandiff/dx).gt.nlabmax)
        if (dx.eq.0.25) then
          dx=0.5
        else if (dx.eq.0.5) then
          dx=1.
        else if (dx.eq.1.) then
          dx=2.5
        else if (dx.eq.2.) then
          dx=5.
        else
          dx=dx*2.
        endif
      enddo

      dx=dx*10.**nexpdiff

      write(chreal,'(1pe12.5)')xmax
      read(chreal(10:12),'(i3)')nexpmax
      read(chreal(1:8),*)rmanmax
      xhigh=int(rmanmax)*10.**nexpmax

      if (xmin*xmax.lt.0.0) xhigh=0.0

      x=xhigh

      if (xhigh.lt.xmin) then
        do while(x.lt.xmin+dx/1000.)
          x=x+dx
        enddo
        xhigh=x
      endif

      x=xhigh

      nlab=0
      do while(x.le.xmax+dx/1000.)
        x=x+dx
        nlab=nlab+1
      enddo

      x=xhigh-dx
      do while(x.ge.xmin-dx/1000.)
        x=x-dx
        nlab=nlab+1
      enddo

      xlab(1)=x+dx
      do il=2,nlab
        xlab(il)=xlab(il-1)+dx
      enddo

      do  il=1,nlab
        if (abs(xlab(il)/dx).lt.eps) xlab(il)=0.0
      enddo

      return
      end
*CMZ :  0.01/02 22/09/2014  13.46.15  by  Michael Scheer
*CMZ :  0.00/06 22/08/2014  15.23.52  by  Michael Scheer
*CMZ :  0.00/04 07/08/2014  11.49.49  by  Michael Scheer
*CMZ :  0.00/02 10/07/2014  15.55.54  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_set_index_height(siz)

      implicit none

      real siz

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      chhe_index_ps=siz

      if (ihigzmode_ps.eq.0) then
        write(cline_ps,*)'/Helvetica findfont',chhe_index_ps*scaletxt_ps,
     &    ' scalefont setfont'
      else
        write(cline_ps,*)'/Helvetica findfont',chhe_index_ps*scaletxt_ps*1.3,
     &    ' scalefont setfont'
      endif

      call mshplt_fill_buff(cline_ps)

      return
      end
*CMZ :  1.01/02 25/09/2014  14.49.04  by  Michael Scheer
*CMZ :  0.01/03 23/09/2014  11.31.17  by  Michael Scheer
*CMZ :  0.01/02 22/09/2014  17.02.18  by  Michael Scheer
*CMZ :  0.00/04 11/08/2014  17.53.18  by  Michael Scheer
*CMZ :  0.00/02 07/07/2014  12.22.09  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_text_relative_position(dx,dy,chtext,chhe)

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      real dx,dy,chhe,chheo
      character(*) chtext

      chheo=chhe_ps

      call mshplt_flush_buff
      call mshplt_fill_buff('gsave currentpoint pop 0 translate')
      call mshplt_set_character_height(chhe)
      write(cline_ps,*)dx,dy,' rmoveto ',tang_ps,' rotate 0 0 ('
     &  //chtext(1:len_trim(chtext))//
     &  ') show ',-tang_ps,' rotate grestore'
      call mshplt_fill_buff(cline_ps)
      call mshplt_set_character_height(chheo)

      return
      end
*CMZ :  0.01/02 22/09/2014  17.05.14  by  Michael Scheer
*CMZ :  0.00/04 11/08/2014  17.53.18  by  Michael Scheer
*CMZ :  0.00/02 07/07/2014  12.22.09  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_text_upper(chtext)

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      character(*) chtext
      real chhe

      chhe=chhe_index_ps
      call mshplt_text_relative_position(chhe/4.,chhe,chtext,chhe)

      return
      end
*CMZ :  0.01/02 22/09/2014  17.03.34  by  Michael Scheer
*CMZ :  0.00/06 22/08/2014  15.23.52  by  Michael Scheer
*CMZ :  0.00/04 07/08/2014  11.49.49  by  Michael Scheer
*CMZ :  0.00/02 10/07/2014  15.55.54  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_get_character_height(siz)

      implicit none

      real siz

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      siz=chhe_ps

      return
      end
*CMZ :  0.01/02 22/09/2014  17.10.33  by  Michael Scheer
*CMZ :  0.00/04 11/08/2014  17.53.18  by  Michael Scheer
*CMZ :  0.00/02 07/07/2014  12.22.09  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_text_lower(chtext)

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      character(*) chtext
      real chhe

      chhe=chhe_index_ps
      call mshplt_text_relative_position(chhe/4.,-chhe/2.,chtext,chhe)

      return
      end
*CMZ :  0.01/03 23/09/2014  12.13.27  by  Michael Scheer
*CMZ :  0.01/02 22/09/2014  17.12.02  by  Michael Scheer
*CMZ :  0.00/04 11/08/2014  17.53.18  by  Michael Scheer
*CMZ :  0.00/02 07/07/2014  12.22.09  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_text_upper_lower(chup,chlow)

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      character(*) chup,chlow

      call mshplt_text_upper(chup)
      call mshplt_text_lower(chlow)

      return
      end
*CMZ :  0.01/03 22/09/2014  20.55.58  by  Michael Scheer
*CMZ :  0.01/02 05/09/2014  15.41.43  by  Michael Scheer
*CMZ :  0.00/03 04/08/2014  13.25.01  by  Michael Scheer
*CMZ :  0.00/02 09/07/2014  14.48.35  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_line_raw(x1,y1,x2,y2)

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      real x(2),y(2),x1,y1,x2,y2

      x(1)=x1
      y(1)=y1
      x(2)=x2
      y(2)=y2

      call mshplt_pline_raw(2,x,y)

      return
      end
*CMZ :  0.01/03 22/09/2014  20.55.58  by  Michael Scheer
*CMZ :  0.01/02 11/09/2014  12.03.40  by  Michael Scheer
*CMZ :  0.00/04 14/08/2014  12.58.41  by  Michael Scheer
*CMZ :  0.00/03 06/08/2014  09.32.19  by  Michael Scheer
*CMZ :  0.00/02 09/07/2014  14.48.35  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_pline_raw(n,x,y)

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      real x(*),y(*),x1,y1,x2,y2
      integer n,i,ic,ir,ig,ib,kc,kr,kg,kb
      character(2048) cline

      if (n.le.1) return

      call mshplt_get_line_color(ic,ir,ig,ib)
      call mshplt_get_color(kc,kr,kg,kb)
      if (ic.ne.kc.or.ir.ne.kr.or.kg.ne.ig.or.ib.ne.kb) then
        call mshplt_set_line_color(ic,ir,ig,ib)
      endif

      x1=xleft_ps+scalex_ps*(x(1)-wxmin_ps)
      y1=ybottom_ps+scaley_ps*(y(1)-wymin_ps)

      write(cline,*)x1,y1,' moveto'
      call mshplt_fill_buff(cline)

      do i=2,n
        x2=xleft_ps+scalex_ps*(x(i)-wxmin_ps)
        y2=ybottom_ps+scaley_ps*(y(i)-wymin_ps)
        write(cline,*)x2,y2,' lineto'
        call mshplt_fill_buff(cline)
      enddo

      itouched_ps=1

      write(cline,*)'stroke'
      call mshplt_fill_buff(cline)

      return
      end
*CMZ :  0.01/03 22/09/2014  20.55.58  by  Michael Scheer
*CMZ :  0.01/02 05/09/2014  15.41.43  by  Michael Scheer
*CMZ :  0.00/03 06/08/2014  09.18.49  by  Michael Scheer
*CMZ :  0.00/02 09/07/2014  14.48.35  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_pline_raw_closed(n,x,y)

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      integer n
      real x(*),y(*)

      if (n.le.1) return

      call mshplt_pline_raw(n,x,y)
      call mshplt_line_raw(x(n),y(n),x(1),y(1))

      return
      end
*CMZ :  0.01/03 23/09/2014  12.07.12  by  Michael Scheer
*CMZ :  0.01/00 24/08/2014  13.17.26  by  Michael Scheer
*CMZ :  0.00/06 20/08/2014  12.49.55  by  Michael Scheer
*CMZ :  0.00/04 14/08/2014  12.58.41  by  Michael Scheer
*CMZ :  0.00/02 09/07/2014  16.17.24  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_text_raw(x,y,text)

      implicit none

*KEEP,mshplt.
      real
     &  pttomm_ps,pttocm_ps, !convert from pt to mm or cm respectively
     &  scale_ps, isscale_ps, !current scale to convert from pt
     &  scalex_ps, !current scale to convert x to pad system
     &  scaley_ps, !current scale to convert y to pad system
     &  xold_ps,yold_ps, !old pen position on pad
     &  x_ps,y_ps, !current pen position on pad
     &  rlinewidth_ps !line width on pad
     & ,chhe_ps, chhe_index_ps ! character height
     & ,ticsiz_ps ! tic height
     & ,ang_ps,tang_ps ! current angle and text angle
     &  ,offgtitx_ps !offset of global title on pad
     &  ,offgtity_ps !offset of global title on pad
     &  ,offdatex_ps !offset of date on pad
     &  ,offdatey_ps !offset of date on pad
     &  ,aLabHeight_ps !character height of axis label
     &  ,xoffexp_ps !xoffset of power term
     &  ,yoffexp_ps !yoffset of power term
     &  ,xlaboff_ps !offset of x-axis label on pad
     &  ,xtitoff_ps !offset of x-axis title on pad
     &  ,ylaboff_ps ! offset of y-axis label on pad
     &  ,ytitoff_ps !offset of y-axis title on pad
     &  ,ztitoff_ps !offset of z-axis title on pad
     &  ,zlaboff_ps ! offset of z-axis label on pad
     &  ,xsiz_ps,ysiz_ps !size on pad of plotting canvas
     &  ,wxmin_ps,wxmax_ps,wymin_ps,wymax_ps !world
     &  ,xleft_ps,xright_ps,ybottom_ps,ytop_ps ! pad on pad
     &  ,xsizorig_ps,ysizorig_ps !size on pad of plotting canvas for paper scale
     &  ,xleftorig_ps,xrightorig_ps,ybottomorig_ps,ytoporig_ps ! orig pad on paper
     &  ,rmtyp_ps(2,0:100) !postscript encoding of characters
     &  ,rmsiz_ps ! marker size
     &  ,phi_ps,theta_ps ! 3d viewing angles
     &  ,pi_ps !Pi
     &  ,xmin3d_ps,xmax3d_ps
     &  ,ymin3d_ps,ymax3d_ps
     &  ,zmin3d_ps,zmax3d_ps
     &  ,xmgl_ps,xmgr_ps,ymgu_ps,ymgl_ps !old hplot variables for margins
     &  ,xwin_ps,ywin_ps ! horizontal and vertical distance of pads
     &  ,ygti_ps !y-offset of global title for hplot
     &  ,gsiz_ps !character size of global title for hplot
     &  ,xcen_ps,ycen_ps !center on canvas
     &  ,scaletxt_ps ! Scaling for text entities
     &  ,ihigzmode_ps !HIGZ-compatibility mode
     &  ,tsiz_ps ! axis-title-size in HIGZ-mode
     &  ,coln_ps,cnLine_ps,cnFill_ps,cnMarker_ps,cnText_ps,cnFrame_ps !color norms
     &  ,clips_ps(2,1000) !clipping range
     &  ,radtodeg_ps
     &  ,xpcorn_ps(8),ypcorn_ps(8) !3d frame box
     &  ,ifbox_ps,xfb_ps(9),yfb_ps(9)

      parameter(pttocm_ps=72/2.54,pttomm_ps=72/25.4)

      integer nmaxfile_ps
      parameter(nmaxfile_ps=10000)

      integer
     &  lun_ps,istat_ps
     &  ,ilabmod_ps !mode off plotting label
     &  ,inolabs_ps !no axis labeling
     &  ,irunviewer_ps !start live viewer
     &  ,iviewinter_ps !interactive stopping of live viewer
     &  ,mtyp_ps ! current marker type
     &  ,icolor_ps, ired_ps,igreen_ps,iblue_ps !color, if <1, rgb values are used
     &  ,ilinestyle_ps !line style: 1: solid, 2: dashed, 3: dotted, 4: dashed-dotted
     &  ,kzone_ps, itouched_ps !current zone and flag of usage
     &  ,nzone_ps,nxzone_ps,nyzone_ps !zones of canvas
     &  ,ifirst_ps !first zone for plotting
     &  ,isamecanvas_ps,isameframe_ps
     &  ,nfile_ps !file counter
     &  ,lunall_ps(nmaxfile_ps) !LUNs of plotfiles
     &  ,kbbxl_ps,kbbxr_ps,kbbyb_ps,kbbyt_ps !bounding box (kbbxl_ps,kbbxr_ps,...)
     &  ,inewpage_ps !new page
     &  ,kFrameColor_ps,kFrameRed_ps,kFrameBlue_ps,kFrameGreen_ps
     &  ,kFillColor_ps,kFillRed_ps,kFillBlue_ps,kFillGreen_ps
     &  ,kLineColor_ps,kLineRed_ps,kLineBlue_ps,kLineGreen_ps
     &  ,kTextColor_ps,kTextRed_ps,kTextBlue_ps,kTextGreen_ps
     &  ,kMarkerColor_ps,kMarkerRed_ps,kMarkerBlue_ps,kMarkerGreen_ps
     &  ,idrawgtit_ps !draw global title on top of page
     &  ,lunbase_ps
     &  ,log10x_ps,log10y_ps,log10z_ps !flags for logarithmic axis
     &  ,nclips_ps ! number of clipping points
     &  ,mode3d_ps
     &  ,isbox_ps,isdate_ps,isoffdate_ps
     &  ,kBox_ps,kDate_ps ! to plot surrounding box and date on top of plots

      character(2048)
     &  file_ps,viewer_ps,cline_ps,filebase_ps,
     &  xtit_ps,ytit_ps,ztit_ps,gtit_ps !axis-title and global title
     &  ,fileall_ps(nmaxfile_ps) !plotfiles

      character(4) chmarker_ps !current marker

      integer nbuffsize_ps,ibuffpos_ps
      parameter(nbuffsize_ps=100000)
      character(2048) chbuff_ps(nbuffsize_ps)
      character(4) chch_ps(1000)

      common/mshpltc/
     &  scale_ps, isscale_ps,scalex_ps,scaley_ps,
     &  xold_ps,yold_ps,x_ps,y_ps,rlinewidth_ps,chhe_ps, chhe_index_ps,
     &  ticsiz_ps,
     &  ang_ps,tang_ps
     & ,xlaboff_ps,ylaboff_ps
     & ,xtitoff_ps,ytitoff_ps
     &  ,ztitoff_ps !offset of z-axis title on pad
     &  ,zlaboff_ps ! offset of z-axis label on pad
     & ,xsiz_ps,ysiz_ps
     &  ,wxmin_ps,wxmax_ps,wymin_ps,wymax_ps !world
     &  ,xleft_ps,xright_ps,ybottom_ps,ytop_ps ! pad on pad
     &  ,ilabmod_ps,rmtyp_ps,rmsiz_ps
     & ,lun_ps,istat_ps,inolabs_ps,irunviewer_ps,iviewinter_ps,mtyp_ps
     &  ,file_ps,viewer_ps,filebase_ps
     &  ,ibuffpos_ps,chbuff_ps,chch_ps,chmarker_ps,cline_ps
     &  ,phi_ps,theta_ps ! 3d viewing angles
     &  ,pi_ps !Pi
     &  ,icolor_ps, ired_ps,igreen_ps,iblue_ps !color, if <1, rgb values are used
     &  ,ilinestyle_ps !line style: 1: solid, 2: dashed, 3: dotted, 4: dashed-dotted
     &  ,xmin3d_ps,xmax3d_ps
     &  ,ymin3d_ps,ymax3d_ps
     &  ,zmin3d_ps,zmax3d_ps
     &  ,offgtitx_ps !offset of global title on pad
     &  ,offgtity_ps !offset of global title on pad
     &  ,offdatex_ps !offset of date on pad
     &  ,offdatey_ps !offset of date on pad
     &  ,xmgl_ps,xmgr_ps,ymgu_ps,ymgl_ps !old hplot variables for margins
     &  ,xwin_ps,ywin_ps ! horizontal and vertical distance of pads
     &  ,xsizorig_ps,ysizorig_ps !size on pad of plotting canvas for paper scale
     &  ,xleftorig_ps,xrightorig_ps,ybottomorig_ps,ytoporig_ps ! orig pad on paper
     &  ,kzone_ps, itouched_ps !current zone and flag of usage
     &  ,nzone_ps,nxzone_ps,nyzone_ps !zones of canvas
     &  ,ifirst_ps !first zone for plotting
     &  ,isamecanvas_ps,isameframe_ps
     &  ,xtit_ps,ytit_ps,ztit_ps,gtit_ps !axis-title and global title
     &  ,nfile_ps !file counter
     &  ,lunall_ps !LUNs of plotfiles
     &  ,fileall_ps !plotfiles
     &  ,kbbxl_ps,kbbxr_ps,kbbyb_ps,kbbyt_ps !bounding box (kbbxl_ps,kbbxr_ps,...)
     &  ,inewpage_ps !new page
     &  ,kBox_ps,kDate_ps ! to plot surrounding box and date on top of plots
     &  ,ygti_ps !y-offset of global title for hplot
     &  ,gsiz_ps !character size of global title for hplot
     &  ,kFrameColor_ps,kFrameRed_ps,kFrameBlue_ps,kFrameGreen_ps
     &  ,kLineColor_ps,kLineRed_ps,kLineBlue_ps,kLineGreen_ps
     &  ,kFillColor_ps,kFillRed_ps,kFillBlue_ps,kFillGreen_ps
     &  ,kTextColor_ps,kTextRed_ps,kTextBlue_ps,kTextGreen_ps
     &  ,kMarkerColor_ps,kMarkerRed_ps,kMarkerBlue_ps,kMarkerGreen_ps
     &  ,xcen_ps,ycen_ps !center on canvas
     &  ,scaletxt_ps ! Scaling for text entities
     &  ,ihigzmode_ps !HIGZ-compatibility mode
     &  ,idrawgtit_ps !draw global title on top of page
     &  ,tsiz_ps ! axis-title-size in HIGZ-mode
     &  ,lunbase_ps
     &  ,coln_ps,cnLine_ps,cnFill_ps,cnMarker_ps,cnText_ps,cnFrame_ps !color norms
     &  ,log10x_ps,log10y_ps,log10z_ps !flags for logarithmic axis
     &  ,aLabHeight_ps !character height of axis label
     &  ,clips_ps,nclips_ps
     &  ,radtodeg_ps
     &  ,mode3d_ps
     &  ,isbox_ps,isdate_ps,isoffdate_ps
     &  ,xoffexp_ps !xoffset of power term
     &  ,yoffexp_ps !yoffset of power term
     &  ,xpcorn_ps,ypcorn_ps
     &  ,ifbox_ps,xfb_ps,yfb_ps
*KEND.

      real x1,y1,x,y
      integer ic,ir,ig,ib
      character(*) text
      character(2048) cline

      itouched_ps=1

      if (icolor_ps.ne.kTextColor_ps) then
        call mshplt_get_text_color(ic,ir,ig,ib)
        call mshplt_set_text_color(ic,ir,ig,ib)
      endif

      x1=xleft_ps+scalex_ps*(x-wxmin_ps)
      y1=ybottom_ps+scaley_ps*(y-wymin_ps)

      write(cline,*)
     &  x1,y1,' moveto ',
     &  tang_ps,' rotate ',
     &  '(',text(1:len_trim(text)),') show',
     &  -tang_ps,' rotate'
      call mshplt_fill_buff(cline)

      return
      end
*CMZ :  0.01/03 23/09/2014  08.55.40  by  Michael Scheer
*CMZ :  0.01/02 05/09/2014  15.41.43  by  Michael Scheer
*CMZ :  0.00/04 14/08/2014  12.58.41  by  Michael Scheer
*CMZ :  0.00/03 06/08/2014  09.32.19  by  Michael Scheer
*CMZ :  0.00/02 09/07/2014  14.48.35  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_get_clip(n,x,y)

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      real x(*),y(*)
      integer n,i

      if (n.lt.0) then
        n=nclips_ps
        return
      endif

      n=nclips_ps
      do i=1,nclips_ps
        x(n)=clips_ps(1,n)
        y(n)=clips_ps(2,n)
      enddo

      return
      end
*CMZ :  1.00/01 23/09/2014  19.03.09  by  Michael Scheer
*CMZ :  0.01/03 23/09/2014  09.38.20  by  Michael Scheer
*CMZ :  0.01/02 05/09/2014  15.41.43  by  Michael Scheer
*CMZ :  0.01/00 24/08/2014  14.25.20  by  Michael Scheer
*CMZ :  0.00/04 14/08/2014  12.58.41  by  Michael Scheer
*CMZ :  0.00/02 10/07/2014  15.51.37  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_marker_raw(n,x,y)

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      real x(*),y(*),x1,y1,yoffset,r
      integer n,i,mtyp
      character(2048) cline
      character(4) chmark

      itouched_ps=1

      if (mtyp_ps.eq.9) then
        do i=1,n
          r=rmsiz_ps*scaletxt_ps
          call mshplt_circle_raw(x(i),y(i),r)
        enddo
        return
      endif

      mtyp=nint(rmtyp_ps(1,mtyp_ps))
      chmark=chch_ps(mtyp)
      yoffset=rmtyp_ps(2,mtyp_ps)*rmsiz_ps*scaletxt_ps

      write(cline,'(a,f12.5,a)')'/Symbol findfont ',rmsiz_ps*scaletxt_ps,
     &  ' scalefont setfont'
      call mshplt_fill_buff(cline)

      do i=1,n

        x1=xleft_ps+scalex_ps*(x(i)-wxmin_ps)
        y1=ybottom_ps+scaley_ps*(y(i)-wymin_ps)

        write(cline,'(2f12.5,a,f12.5,a)')x1,y1,
     &    ' moveto ('//chmark(1:len_trim(chmark))
     &    //') dup stringwidth pop 2 div neg 0 rmoveto 0 ',
     &    yoffset,' rmoveto show'
        call mshplt_fill_buff(cline)

      enddo

      write(cline,'(a,f12.5,a)')'/Helvetica findfont ',chhe_ps,
     &  ' scalefont setfont'
      call mshplt_fill_buff(cline)

      return
      end
*CMZ :  0.01/03 23/09/2014  09.13.55  by  Michael Scheer
*CMZ :  0.01/02 15/09/2014  15.09.38  by  Michael Scheer
*CMZ :  0.00/04 11/08/2014  18.23.32  by  Michael Scheer
*CMZ :  0.00/02 07/07/2014  12.22.09  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_circle_raw(x,y,r)

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      real x,y,r,x1,y1

      x1=xleft_ps+scalex_ps*(x-wxmin_ps)
      y1=ybottom_ps+scaley_ps*(y-wymin_ps)

      write(cline_ps,*)'newpath ',x1,y1,r*2.*chhe_ps,' 0 360 arc stroke'
      call mshplt_fill_buff(cline_ps)

      return
      end
*CMZ :  1.02/00 30/09/2014  11.27.05  by  Michael Scheer
*CMZ :  1.00/01 24/09/2014  09.11.35  by  Michael Scheer
*CMZ :  0.01/03 23/09/2014  13.49.18  by  Michael Scheer
*CMZ :  0.01/02 05/09/2014  15.41.43  by  Michael Scheer
*CMZ :  0.00/03 05/08/2014  15.53.35  by  Michael Scheer
*CMZ :  0.00/02 08/07/2014  11.10.50  by  Michael Scheer
*CMZ :  0.00/01 07/07/2014  12.00.23  by  Michael Scheer
*CMZ : 00.00/02 04/07/2014  21.11.57  by  Michael Scheer
*-- Author :    Michael Scheer   11/06/2014
      subroutine mshplt_view_axis(vxmin,vxmax,vymin,vymax,smin,smax,title,
     &  iticside,anglab,titang,offtit,offlab)

c Draw an axis from (vxmin,vxmax) to (vymin,vymax) in view system
c with labeling from smin to smax, and title title

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      real
     &  vxmin,vxmax,vymin,vymax,anglab,
     &  smin,smax,titang,offtit,offlab

      real sizlab(100),xlabrel(100),ylaboffset(100),anglabrel(100),ticlen(100),
     &  ticposrel(100),titoff,titangrel,ticangrel(100),titsiz,titposrel

      integer iticside
      integer nlab,ntic

      character(*) title
      character(12) chlab(100)

      nlab=-100
      sizlab=aLabHeight_ps*(1-inolabs_ps)
      ylaboffset=offlab
      anglabrel=anglab
      ticlen=ticsiz_ps
      if (iticside.lt.0) ticlen=-ticlen
      ticangrel=90.
      titsiz=tsiz_ps*(1-inolabs_ps)
      titoff=offtit
      titangrel=titang
      titposrel=0.5

      call mshplt_axis_taylor(
     &  vxmin,vxmax,vymin,vymax,smin,smax,
     &  nlab,chlab,sizlab,xlabrel,ylaboffset,
     &  xoffexp_ps,yoffexp_ps,
     &  anglabrel,
     &  ntic,ticlen,ticposrel,ticangrel,
     &  titsiz,titposrel,titoff,titangrel,title
     &  )

      return
      end
*CMZ :  1.01/00 24/09/2014  15.22.29  by  Michael Scheer
*CMZ :  1.00/01 23/09/2014  17.15.48  by  Michael Scheer
*CMZ :  0.01/02 05/09/2014  15.41.43  by  Michael Scheer
*CMZ :  0.00/04 11/08/2014  17.53.18  by  Michael Scheer
*CMZ :  0.00/02 07/07/2014  12.22.09  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_marker_single(x,y)

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      real x(1),y(1)

      call mshplt_marker(1,x,y)

      return
      end
*CMZ :  1.02/00 02/10/2014  13.30.19  by  Michael Scheer
*CMZ :  1.00/01 23/09/2014  17.16.19  by  Michael Scheer
*CMZ :  0.01/02 05/09/2014  15.41.43  by  Michael Scheer
*CMZ :  0.00/04 11/08/2014  17.53.18  by  Michael Scheer
*CMZ :  0.00/02 07/07/2014  12.22.09  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_marker_single_raw(x,y)

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      real x(1),y(1)

      call mshplt_marker_raw(1,x(1),y(1))

      return
      end
*CMZ :  1.02/00 01/10/2014  14.10.39  by  Michael Scheer
*CMZ :  1.01/02 26/09/2014  13.40.07  by  Michael Scheer
*CMZ :  1.00/01 24/09/2014  13.22.59  by  Michael Scheer
*CMZ :  0.01/03 22/09/2014  21.17.48  by  Michael Scheer
*CMZ :  0.01/02 22/09/2014  19.54.46  by  Michael Scheer
*-- Author :    Michael Scheer   08/09/2014
      subroutine mshplt_view_log_axis(
     &  xorig, yorig, xylength, smin, smax, ang,
     &  chhe, ticheight, title,iticside,ilabside,reltitang,offtit,
     &  rellabang,offlab)

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      real xorig, yorig, xylength, smin, smax, ang,chhe,ticheight
      real titang,offtit,offlab,sminl,smaxl,reltitang,rellabang,
     &  x(2),y(2),aLabSide,xdec,ydec,nexp,tango,chheo,alabang

      integer iticside,ilabside,j,l,idecmin,idecmax
      character(*) title

      integer ntic,nlab,ndec,i
      real xmin,xmax,ymin,ymax,sinang,cosang,dxy,tside,costit,sintit

      character(8) c8

      call mshplt_flush_buff

      write(lun_ps,'(a)')'% begin of mshplt_view_log_axis'

      tango=ang_ps
      chheo=chhe_ps

      if (iticside.lt.0) then
        tside=1.
      else if (iticside.gt.0) then
        tside=-1.
      else
        tside=0.
      endif

      if (ilabside.lt.0) then
        aLabSide=-1.
      else if (ilabside.gt.0) then
        aLabSide=0.5
      else
        aLabSide=0.
      endif

      if (smin.le.0.0) then
        sminl=-30.
      else
        sminl=alog10(smin)
      endif

      if (smax.le.0.0) then
        smaxl=30.
      else
        smaxl=alog10(smax)
      endif

      cosang=cos(ang/radtodeg_ps)
      sinang=sin(ang/radtodeg_ps)

c      wtov=1./(scalex_ps*cosang+scaley_ps*sinang)

      xmin=xorig
      xmax=xorig+cosang*xylength
      ymin=yorig
      ymax=yorig+sinang*xylength

      call mshplt_line_raw(xmin,ymin,xmax,ymax)

      costit=cos((ang+rellabang)/radtodeg_ps)
      sintit=sin((ang+rellabang)/radtodeg_ps)

      call mshplt_view_to_world(0.,0.,x(1),y(1))
      call mshplt_view_to_world(costit,sintit,x(2),y(2))
      alabang=atan2(y(2)-y(1),x(2)-x(1))*radtodeg_ps

      costit=cos((ang+reltitang)/radtodeg_ps)
      sintit=sin((ang+reltitang)/radtodeg_ps)

      call mshplt_view_to_world(0.,0.,x(1),y(1))
      call mshplt_view_to_world(costit,sintit,x(2),y(2))
      titang=atan2(y(2)-y(1),x(2)-x(1))*radtodeg_ps

      idecmin=int(sminl)
      idecmax=int(smaxl)

      if (idecmax.lt.30.and.10**idecmax.lt.smax) then
        idecmax=idecmax+1
      endif

      ndec=idecmax-idecmin

      nlab=ndec*9
      ntic=nlab

      dxy=xylength/ndec

      write(c8,'(i8)')max(abs(idecmin),abs(idecmax))
      do l=1,8
        if(c8(l:l).ne.' ') then
          nexp=8-l+1
          goto 1
        endif
      enddo

 1    continue

      call mshplt_set_text_angle(alabang)
      call mshplt_set_character_height(aLabHeight_ps)

      do i=1,ndec

        xdec=xorig+i*dxy*cosang
        ydec=yorig+i*dxy*sinang

        x(1)=xdec
        x(2)=x(1)+ticheight*tside*sinang/scalex_ps
        y(1)=ydec
        y(2)=y(1)-ticheight*tside*cosang/scaley_ps

        call mshplt_pline_raw(2,x,y)

        if (ilabside.ne.0.and.aLabHeight_ps.gt.0.0) then

          x(2)=x(2)
     &      -iLabSide*sign(1.,90.0-theta_ps)*sinang/scalex_ps*(offlab+(0.8*nexp*aLabHeight_ps))

          y(2)=ydec
     &      -cosang/scaley_ps*(offlab+(0.8*nexp*aLabHeight_ps))

          call mshplt_text_raw(x(2),y(2),"10")

          write(c8,'(i8)')idecmin+i
          do l=1,8
            if(c8(l:l).ne.' ') goto 11
          enddo
 11       continue

          call mshplt_text_upper(c8(l:8))

        endif

        do j=2,9
          x(1)=xorig+((i-1)*dxy+alog10(float(j))*dxy)*cosang
          x(2)=x(1)+ticheight/scalex_ps*sinang*tside*0.8
          y(1)=yorig+((i-1)*dxy+alog10(float(j))*dxy)*sinang
          y(2)=y(1)-ticheight/scaley_ps*cosang*tside*0.8
          call mshplt_pline_raw(2,x,y)
        enddo

      enddo !ndec

      call mshplt_set_character_height(chhe)

      call mshplt_set_text_angle(titang)

      if (iLabSide.ne.0.and.aLabHeight_ps.ne.0.0) then

        x(1)=xorig
     &    +0.5*xylength*cosang
     &    -len_trim(title)*chhe/4.*costit/scalex_ps
     &    +iLabSide*(offtit+0.8*nexp*aLabHeight_ps)*abs(sintit)/scalex_ps

        y(1)=yorig
     &    +0.5*xylength*sinang
     &    -len_trim(title)*chhe/4.*sintit/scaley_ps
     &    -(offtit+0.8*nexp*aLabHeight_ps)*costit/scaley_ps

      else
        goto 9999
      endif

      call mshplt_set_character_height(chhe)
      call mshplt_text_raw(x(1),y(1),title(1:len_trim(title)))
      call mshplt_set_character_height(chheo)
      call mshplt_set_text_angle(tango)

 9999 call mshplt_flush_buff
      write(lun_ps,'(a)')'% end of mshplt_log_axis'

      return
      end
*CMZ :  1.03/01 09/10/2014  15.43.42  by  Michael Scheer
*CMZ :  1.01/00 24/09/2014  14.06.28  by  Michael Scheer
*CMZ :  0.01/02 05/09/2014  15.41.43  by  Michael Scheer
*CMZ :  0.00/04 11/08/2014  17.53.18  by  Michael Scheer
*CMZ :  0.00/02 07/07/2014  12.22.09  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_set_date_position(offx,offy)

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      real offx,offy

      offdatex_ps=offx
      offdatey_ps=offy

      isoffdate_ps=1

      return
      end
*CMZ :  1.01/02 26/09/2014  14.00.09  by  Michael Scheer
*CMZ :  1.01/00 24/09/2014  14.46.38  by  Michael Scheer
*CMZ :  1.00/01 23/09/2014  14.36.32  by  Michael Scheer
*CMZ :  0.01/02 05/09/2014  15.41.43  by  Michael Scheer
*CMZ :  0.00/03 06/08/2014  11.12.41  by  Michael Scheer
*CMZ :  0.00/02 09/07/2014  16.03.13  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_pline_3d(n,x,y,z)

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      real x(*),y(*),z(*)
      real, dimension(:), allocatable :: xn,yn,zn,xx,yy
      integer n,i

      allocate(xn(n))
      allocate(yn(n))
      allocate(zn(n))
      allocate(xx(n))
      allocate(yy(n))

      mode3d_ps=1

      do i=1,n
        if (log10x_ps.eq.0) then
          xn(i)=-0.5+(x(i)-xmin3d_ps)/(xmax3d_ps-xmin3d_ps)
        else
          xn(i)=-0.5+(alog10(x(i))-xmin3d_ps)/(xmax3d_ps-xmin3d_ps)
        endif
        if (log10y_ps.eq.0) then
          yn(i)=-0.5+(y(i)-ymin3d_ps)/(ymax3d_ps-ymin3d_ps)
        else
          yn(i)=-0.5+(alog10(y(i))-ymin3d_ps)/(ymax3d_ps-ymin3d_ps)
        endif
        if (log10z_ps.eq.0) then
          zn(i)=-0.5+(z(i)-zmin3d_ps)/(zmax3d_ps-zmin3d_ps)
        else
          zn(i)=-0.5+(alog10(z(i))-zmin3d_ps)/(zmax3d_ps-zmin3d_ps)
        endif
      enddo

      call mshplt_3dto2d(n,xn,yn,zn,xx,yy)
      call mshplt_pline(n,xx,yy)

      mode3d_ps=0

      deallocate(xn)
      deallocate(yn)
      deallocate(zn)
      deallocate(xx)
      deallocate(yy)

      return
      end
*CMZ :  1.02/00 30/09/2014  13.52.50  by  Michael Scheer
*CMZ :  1.01/01 25/09/2014  09.16.46  by  Michael Scheer
*CMZ :  1.01/00 24/09/2014  14.46.38  by  Michael Scheer
*CMZ :  1.00/01 23/09/2014  14.36.32  by  Michael Scheer
*CMZ :  0.01/02 05/09/2014  15.41.43  by  Michael Scheer
*CMZ :  0.00/03 06/08/2014  11.12.41  by  Michael Scheer
*CMZ :  0.00/02 09/07/2014  16.03.13  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_marker_3d(n,x,y,z)

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      real x(*),y(*),z(*)
      real, dimension(:), allocatable :: xn,yn,zn,xx,yy
      integer n,i

      allocate(xn(n))
      allocate(yn(n))
      allocate(zn(n))
      allocate(xx(n))
      allocate(yy(n))

      do i=1,n
        if (log10x_ps.eq.0) then
          xn(i)=-0.5+(x(i)-xmin3d_ps)/(xmax3d_ps-xmin3d_ps)
        else
          xn(i)=-0.5+(alog10(x(i))-xmin3d_ps)/(xmax3d_ps-xmin3d_ps)
        endif
        if (log10y_ps.eq.0) then
          yn(i)=-0.5+(y(i)-ymin3d_ps)/(ymax3d_ps-ymin3d_ps)
        else
          yn(i)=-0.5+(alog10(y(i))-ymin3d_ps)/(ymax3d_ps-ymin3d_ps)
        endif
        if (log10z_ps.eq.0) then
          zn(i)=-0.5+(z(i)-zmin3d_ps)/(zmax3d_ps-zmin3d_ps)
        else
          zn(i)=-0.5+(alog10(z(i))-zmin3d_ps)/(zmax3d_ps-zmin3d_ps)
        endif
      enddo

      call mshplt_3dto2d(n,xn,yn,zn,xx,yy)
      call mshplt_marker_raw(n,xx,yy)

      deallocate(xn)
      deallocate(yn)
      deallocate(zn)
      deallocate(xx)
      deallocate(yy)

      return
      end
*CMZ :  1.01/01 25/09/2014  09.02.51  by  Michael Scheer
*CMZ :  1.01/00 24/09/2014  15.22.29  by  Michael Scheer
*CMZ :  1.00/01 23/09/2014  17.15.48  by  Michael Scheer
*CMZ :  0.01/02 05/09/2014  15.41.43  by  Michael Scheer
*CMZ :  0.00/04 11/08/2014  17.53.18  by  Michael Scheer
*CMZ :  0.00/02 07/07/2014  12.22.09  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_marker_single_3d(x,y,z)

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      real x(1),y(1),z(1)

      call mshplt_marker_3d(1,x,y,z)

      return
      end
*CMZ :  1.02/00 02/10/2014  13.33.43  by  Michael Scheer
*CMZ :  1.01/02 25/09/2014  09.19.29  by  Michael Scheer
*CMZ :  1.01/01 25/09/2014  09.15.18  by  Michael Scheer
*CMZ :  1.01/00 24/09/2014  14.46.38  by  Michael Scheer
*CMZ :  1.00/01 23/09/2014  14.36.32  by  Michael Scheer
*CMZ :  0.01/02 05/09/2014  15.41.43  by  Michael Scheer
*CMZ :  0.00/03 06/08/2014  11.12.41  by  Michael Scheer
*CMZ :  0.00/02 09/07/2014  16.03.13  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_text_3d(x,y,z,text)

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      real x,y,z,xn(1),yn(1),zn(1),xx(1),yy(1)
      character(*) text

      xn=-0.5+(x-xmin3d_ps)/(xmax3d_ps-xmin3d_ps)
      yn=-0.5+(y-ymin3d_ps)/(ymax3d_ps-ymin3d_ps)

      if (log10z_ps.eq.0) then
        zn=-0.5+(z-zmin3d_ps)/(zmax3d_ps-zmin3d_ps)
      else
        zn=-0.5+(alog10(z)-zmin3d_ps)/(zmax3d_ps-zmin3d_ps)
      endif

      call mshplt_3dto2d(1,xn,yn,zn,xx,yy)
      call mshplt_text_raw(xx(1),yy(1),text)

      return
      end
*CMZ :  1.02/00 29/09/2014  13.44.41  by  Michael Scheer
*CMZ :  1.01/02 27/09/2014  16.19.04  by  Michael Scheer
*CMZ :  0.01/02 05/09/2014  15.41.43  by  Michael Scheer
*CMZ :  0.00/04 11/08/2014  17.53.18  by  Michael Scheer
*CMZ :  0.00/02 07/07/2014  12.22.09  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_set_box(ibox)

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      integer ibox

      kbox_ps=ibox

      isbox_ps=1

      return
      end
*CMZ :  1.03/00 06/10/2014  18.03.21  by  Michael Scheer
*CMZ :  1.02/00 29/09/2014  10.18.31  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_filled_area(n,x,y,ioutlined)

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      real x(*),y(*),dx,dy
      integer ioutlined,i,n,icolo,iro,igo,ibo

      if (n.lt.3) return

      call mshplt_get_color(icolo,iro,igo,ibo)

      call mshplt_set_color(
     &  klinecolor_ps,klinered_ps,klinegreen_ps,klineblue_ps)

      dx=x(1)
      dy=y(1)
      call mshplt_view_to_world(dx,dy,dx,dy)
      write(cline_ps,*)'newpath ',dx,dy,' moveto'
      call mshplt_fill_buff(cline_ps)

      do i=1,n-1
        dx=(x(i+1)-x(i))*scalex_ps
        dy=(y(i+1)-y(i))*scaley_ps
        write(cline_ps,*)dx,dy,' rlineto'
        call mshplt_fill_buff(cline_ps)
      enddo

      call mshplt_fill_buff('closepath')

      if (ioutlined.lt.0) then
        call mshplt_fill_buff('stroke')
        goto 9999
      endif

      if (ioutlined.eq.0) then
        call mshplt_set_color(
     &    kFillColor_ps,kFillRed_ps,kFillGreen_ps,kFillBlue_ps)
        call mshplt_fill_buff('fill')
      else
        call mshplt_fill_buff('gsave')
        call mshplt_set_color(
     &    kFillColor_ps,kFillRed_ps,kFillGreen_ps,kFillBlue_ps)
        call mshplt_fill_buff('fill grestore')
        call mshplt_set_color(
     &    kLineColor_ps,kLineRed_ps,kLineGreen_ps,kLineBlue_ps)
        call mshplt_fill_buff('stroke')
      endif

9999  call mshplt_set_color(icolo,iro,igo,ibo)

      return
      end
*CMZ :  1.02/00 29/09/2014  10.07.08  by  Michael Scheer
*CMZ :  0.01/02 05/09/2014  15.41.43  by  Michael Scheer
*CMZ :  0.00/06 19/08/2014  15.09.19  by  Michael Scheer
*CMZ :  0.00/04 07/08/2014  15.52.16  by  Michael Scheer
*CMZ :  0.00/03 04/08/2014  15.42.32  by  Michael Scheer
*CMZ :  0.00/02 07/07/2014  12.22.09  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_set_fill_color(icolor,ired,igreen,iblue)

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      integer icolor,ired,igreen,iblue

      kFillColor_ps=icolor

      if (icolor.gt.0) then
        if (icolor.eq.1) then
          kFillRed_ps=0
          kFillGreen_ps=0
          kFillBlue_ps=0
        else if (icolor.eq.2) then
          kFillRed_ps=1
          kFillGreen_ps=0
          kFillBlue_ps=0
        else if (icolor.eq.3) then
          kFillRed_ps=0
          kFillGreen_ps=1
          kFillBlue_ps=0
        else if (icolor.eq.4) then
          kFillRed_ps=0
          kFillGreen_ps=0
          kFillBlue_ps=1
        else if (icolor.eq.5) then
          kFillRed_ps=1
          kFillGreen_ps=1
          kFillBlue_ps=0
        else if (icolor.eq.6) then
          kFillRed_ps=1
          kFillGreen_ps=0
          kFillBlue_ps=1
        else if (icolor.eq.7) then
          kFillRed_ps=0
          kFillGreen_ps=1
          kFillBlue_ps=1
        else if (icolor.eq.8) then
          kFillRed_ps=35
          kFillGreen_ps=85
          kFillBlue_ps=33
        endif
      else
        kFillColor_ps=-1
        kFillRed_ps=ired
        kFillGreen_ps=igreen
        kFillBlue_ps=iblue
      endif !icolor.gt.0

      cnFill_ps=sqrt(float(kFillRed_ps**2+kFillGreen_ps**2+kFillBlue_ps**2))
      if (cnFill_ps.le.0.) cnFill_ps=1.

      write(cline_ps,*)
     &  kFillRed_ps/cnFill_ps,
     &  kFillGreen_ps/cnFill_ps,
     &  kFillBlue_ps/cnFill_ps,
     &  ' setrgbcolor'
      call mshplt_fill_buff(cline_ps)

      return
      end
*CMZ :  1.02/00 29/09/2014  10.09.09  by  Michael Scheer
*CMZ :  0.01/02 05/09/2014  15.41.43  by  Michael Scheer
*CMZ :  0.00/06 19/08/2014  15.01.54  by  Michael Scheer
*CMZ :  0.00/04 08/08/2014  16.34.41  by  Michael Scheer
*CMZ :  0.00/03 04/08/2014  15.42.32  by  Michael Scheer
*CMZ :  0.00/02 07/07/2014  12.22.09  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_get_fill_color(icolor,ired,igreen,iblue)

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      integer icolor,ired,igreen,iblue

      icolor=kFillColor_ps
      ired=kFillRed_ps
      igreen=kFillGreen_ps
      iblue=kFillBlue_ps

      return
      end
*CMZ :  1.02/00 29/09/2014  11.01.39  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_filled_mesh_3d(n,x,y,z,ioutlined)

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      real x(*),y(*),z(*),xp(4),yp(4),xn(4),yn(4),zn(4)
      integer ioutlined,n,i

      if (n.gt.4.or.n.lt.3) return

      do i=1,n
        if (log10z_ps.eq.0) then
          zn(i)=-0.5+(z(i)-zmin3d_ps)/(zmax3d_ps-zmin3d_ps)
        else
          zn(i)=-0.5+(alog10(z(i))-zmin3d_ps)/(zmax3d_ps-zmin3d_ps)
        endif
        if (log10x_ps.eq.0) then
          xn(i)=-0.5+(x(i)-xmin3d_ps)/(xmax3d_ps-xmin3d_ps)
        else
          xn(i)=-0.5+(alog10(x(i))-xmin3d_ps)/(xmax3d_ps-xmin3d_ps)
        endif
        if (log10y_ps.eq.0) then
          yn(i)=-0.5+(y(i)-ymin3d_ps)/(ymax3d_ps-ymin3d_ps)
        else
          yn(i)=-0.5+(alog10(y(i))-ymin3d_ps)/(ymax3d_ps-ymin3d_ps)
        endif
      enddo

      call mshplt_3dto2d(n,xn,yn,zn,xp,yp)
      call mshplt_filled_area(n,xp,yp,ioutlined)

      return
      end
*CMZ :  1.02/00 29/09/2014  13.44.30  by  Michael Scheer
*CMZ :  0.01/02 05/09/2014  15.41.43  by  Michael Scheer
*CMZ :  0.00/04 11/08/2014  17.53.18  by  Michael Scheer
*CMZ :  0.00/02 07/07/2014  12.22.09  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_set_date(idate)

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      integer idate

      if (idate.ne.0) then
        kDate_ps=1
      else
        kDate_ps=0
      endif

      isdate_ps=1

      return
      end
*CMZ :  1.02/00 30/09/2014  12.39.33  by  Michael Scheer
*CMZ :  0.01/02 05/09/2014  15.41.43  by  Michael Scheer
*CMZ :  0.00/06 21/08/2014  16.58.38  by  Michael Scheer
*CMZ :  0.00/04 07/08/2014  11.18.06  by  Michael Scheer
*CMZ :  0.00/02 07/07/2014  12.22.09  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_set_power_offset(x,y)

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      real x,y

      xoffexp_ps=x
      yoffexp_ps=y

      return
      end
*CMZ :  1.02/00 30/09/2014  12.40.15  by  Michael Scheer
*CMZ :  0.01/02 05/09/2014  15.41.43  by  Michael Scheer
*CMZ :  0.00/06 21/08/2014  16.58.38  by  Michael Scheer
*CMZ :  0.00/04 07/08/2014  11.18.06  by  Michael Scheer
*CMZ :  0.00/02 07/07/2014  12.22.09  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_get_power_offset(x,y)

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      real x,y

      x=xoffexp_ps
      y=yoffexp_ps

      return
      end
*CMZ :          07/08/2018  09.08.25  by  Michael Scheer
*CMZ :  1.03/00 07/10/2014  10.55.15  by  Michael Scheer
*CMZ :  1.02/01 05/10/2014  14.25.34  by  Michael Scheer
*CMZ :  1.02/00 03/10/2014  12.39.05  by  Michael Scheer
*-- Author :    Michael Scheer   30/09/2014
      subroutine mshplt_surf(nx,xmin,xmax,ny,ymin,ymax,z,chopt)

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      real xmin,xmax,ymin,ymax,xp(4),yp(4),zp(4),dx,dy,z(*),
     &  dz10,zmin,zmax

      integer nx,ny,ix,iy,ioutlined,lenc,i,iz,icoltile,
     &  kco,kro,kbo,kgo,kc,imcol,mco,mro,mgo,mbo,isurf

      character(*) chopt
      character c1
      character(4) c4
      character(400) copt
      integer ib
      equivalence (ib,c1)

      call mshplt_fill_buff('% begin of mshplt_surf')

      ioutlined=0
      icoltile=0
      imcol=0
      isurf=0

      lenc=min(400,len_trim(chopt))

      if (lenc.lt.4) then
        copt='linesurf'
        lenc=8
      else
        copt=chopt(1:min(lenc,400))
      endif

      ib=0
      do i=1,lenc
        c1=copt(i:i)
        if (ib.lt.97) then
          ib=ib+32
          copt(i:i)=c1
        endif
      enddo

      do i=1,lenc
        c4=copt(i:i+4)
        if (c4.eq.'line') ioutlined=1
        if (c4.eq.'tile') then
          icoltile=1
          isurf=1
        endif
        if (c4.eq.'mark') imcol=1
        if (c4.eq.'surf') isurf=1
      enddo

      if (imcol+icoltile+isurf.eq.0) then
        isurf=1
        ioutlined=1
      endif

      call mshplt_get_fill_color(kco,kro,kgo,kbo)
      call mshplt_get_marker_color(mco,mro,mgo,mbo)

c      if (imcol.ne.0) then
c        call mshplt_set_fill_color(-1,1,1,1)
c      endif

      zmin=1.0e30
      zmax=1.0e-30
      iz=0
      do iy=1,ny
        do ix=1,nx
          iz=iz+1
          if (z(iz).lt.zmin) zmin=z(iz)
          if (z(iz).gt.zmax) zmax=z(iz)
        enddo
      enddo

      if (log10z_ps.eq.0) then
        dz10=(zmax-zmin)*1.001/7.
      else
        dz10=(alog10(zmax)-alog10(zmin))*1.001/7.
      endif
      if (dz10.eq.0.0) dz10=1.

      if (nx.lt.1.or.ny.lt.1) then
        return
      else if (nx.eq.1.and.ny.eq.1) then

        xp(1)=xmin
        yp(1)=ymin
        zp(1)=z(1)
        call mshplt_marker_3d(1,xp,yp,zp)

      else if (nx.gt.1.and.ny.eq.1) then
        dx=(xmax-xmin)/(nx-1)
        do ix=1,nx-1
          xp(1)=xmin+(ix-1)*dx
          xp(2)=xp(1)+dx
          yp(1)=ymin
          yp(2)=ymin
          zp(1)=z(ix)
          zp(2)=z(ix+1)
          call mshplt_pline_3d(2,xp,yp,zp)
        enddo

      else if (ny.gt.1.and.nx.eq.1) then
        dy=(ymax-ymin)/(ny-1)
        do iy=1,ny-1
          yp(1)=ymin+(iy-1)*dy
          yp(2)=yp(1)+dy
          xp(1)=xmin
          xp(2)=xmin
          zp(1)=z(iy)
          zp(2)=z(iy+1)
          call mshplt_pline_3d(2,xp,yp,zp)
        enddo

      else

        dx=(xmax-xmin)/(nx-1)
        dy=(ymax-ymin)/(ny-1)

        if (theta_ps.le.90.) then

          if (phi_ps.le.90.) then

            do iy=ny-1,1,-1
              do ix=nx-1,1,-1

                xp(1)=xmin+(ix-1)*dx
                xp(2)=xp(1)+dx
                xp(3)=xp(2)
                xp(4)=xp(1)

                yp(1)=ymin+(iy-1)*dy
                yp(2)=yp(1)
                yp(3)=yp(2)+dy
                yp(4)=yp(3)

                zp(1)=z(nx*(iy-1)+ix)
                zp(2)=z(nx*(iy-1)+ix+1)
                zp(3)=z(nx*iy+ix+1)
                zp(4)=z(nx*iy+ix)

                if (icoltile.ne.0) then
                  if(log10z_ps.eq.0) then
                    kc=mod(int(((zp(1)+zp(2)+zp(3)+zp(4))/4.-zmin)/dz10)+1,8)
                  else
                    kc=mod(
     &                int(
     &                (alog10((zp(1)+zp(2)+zp(3)+zp(4))/4.)-
     &                alog10(zmin))/dz10)+1,8)
                  endif
                  call mshplt_set_fill_color(kc,0,0,0)
                endif

                if (isurf.ne.0) call mshplt_filled_mesh_3d(4,xp,yp,zp,ioutlined)

                if (imcol.ne.0) then
                  do i=1,4
                    if(log10z_ps.eq.0) then
                      kc=mod(int((zp(i)-zmin)/dz10)+1,8)
                    else
                      kc=mod(int((alog10(zp(i))-alog10(zmin))/dz10)+1,8)
                    endif
                    call mshplt_set_marker_color(kc,0,0,0)
                    call mshplt_marker_3d(1,xp(i),yp(i),zp(i))
                  enddo
                endif
              enddo
            enddo

          else if (phi_ps.le.180.) then

            do ix=nx-1,1,-1
              do iy=1,ny-1

                xp(1)=xmin+(ix-1)*dx
                xp(2)=xp(1)+dx
                xp(3)=xp(2)
                xp(4)=xp(1)

                yp(1)=ymin+(iy-1)*dy
                yp(2)=yp(1)
                yp(3)=yp(2)+dy
                yp(4)=yp(3)

                zp(1)=z(nx*(iy-1)+ix)
                zp(2)=z(nx*(iy-1)+ix+1)
                zp(3)=z(nx*iy+ix+1)
                zp(4)=z(nx*iy+ix)

                if (icoltile.ne.0) then
                  kc=mod(int(((zp(1)+zp(2)+zp(3)+zp(4))/4.-zmin)/dz10)+1,8)
                  call mshplt_set_fill_color(kc,0,0,0)
                endif

                if (isurf.ne.0) call mshplt_filled_mesh_3d(4,xp,yp,zp,ioutlined)

                if (imcol.ne.0) then
                  do i=1,4
                    kc=mod(int((zp(i)-zmin)/dz10)+1,8)
                    call mshplt_set_marker_color(kc,0,0,0)
                    call mshplt_marker_3d(1,xp(i),yp(i),zp(i))
                  enddo
                endif
              enddo
            enddo

          else if (phi_ps.le.270.) then

            do iy=1,ny-1
              do ix=1,nx-1

                xp(1)=xmin+(ix-1)*dx
                xp(2)=xp(1)+dx
                xp(3)=xp(2)
                xp(4)=xp(1)

                yp(1)=ymin+(iy-1)*dy
                yp(2)=yp(1)
                yp(3)=yp(2)+dy
                yp(4)=yp(3)

                zp(1)=z(nx*(iy-1)+ix)
                zp(2)=z(nx*(iy-1)+ix+1)
                zp(3)=z(nx*iy+ix+1)
                zp(4)=z(nx*iy+ix)

                if (icoltile.ne.0) then
                  kc=mod(int(((zp(1)+zp(2)+zp(3)+zp(4))/4.-zmin)/dz10)+1,8)
                  call mshplt_set_fill_color(kc,0,0,0)
                endif

                if (isurf.ne.0) call mshplt_filled_mesh_3d(4,xp,yp,zp,ioutlined)

                if (imcol.ne.0) then
                  do i=1,4
                    kc=mod(int((zp(i)-zmin)/dz10)+1,8)
                    call mshplt_set_marker_color(kc,0,0,0)
                    call mshplt_marker_3d(1,xp(i),yp(i),zp(i))
                  enddo
                endif
              enddo
            enddo

          else if (phi_ps.le.360.) then
            do ix=1,nx-1
              do iy=ny-1,1,-1

                xp(1)=xmin+(ix-1)*dx
                xp(2)=xp(1)+dx
                xp(3)=xp(2)
                xp(4)=xp(1)

                yp(1)=ymin+(iy-1)*dy
                yp(2)=yp(1)
                yp(3)=yp(2)+dy
                yp(4)=yp(3)

                zp(1)=z(nx*(iy-1)+ix)
                zp(2)=z(nx*(iy-1)+ix+1)
                zp(3)=z(nx*iy+ix+1)
                zp(4)=z(nx*iy+ix)

                if (icoltile.ne.0) then
                  kc=mod(int(((zp(1)+zp(2)+zp(3)+zp(4))/4.-zmin)/dz10)+1,8)
                  call mshplt_set_fill_color(kc,0,0,0)
                endif

                if (isurf.ne.0) call mshplt_filled_mesh_3d(4,xp,yp,zp,ioutlined)

                if (imcol.ne.0) then
                  do i=1,4
                    kc=mod(int((zp(i)-zmin)/dz10)+1,8)
                    call mshplt_set_marker_color(kc,0,0,0)
                    call mshplt_marker_3d(1,xp(i),yp(i),zp(i))
                  enddo
                endif
              enddo
            enddo

          endif !phi.le.90.

        else if (theta_ps.le.180.) then

          if (phi_ps.le.90.) then

            do iy=1,ny-1
              do ix=1,nx-1

                xp(1)=xmin+(ix-1)*dx
                xp(2)=xp(1)+dx
                xp(3)=xp(2)
                xp(4)=xp(1)

                yp(1)=ymin+(iy-1)*dy
                yp(2)=yp(1)
                yp(3)=yp(2)+dy
                yp(4)=yp(3)

                zp(1)=z(nx*(iy-1)+ix)
                zp(2)=z(nx*(iy-1)+ix+1)
                zp(3)=z(nx*iy+ix+1)
                zp(4)=z(nx*iy+ix)

                if (icoltile.ne.0) then
                  kc=mod(int(((zp(1)+zp(2)+zp(3)+zp(4))/4.-zmin)/dz10)+1,8)
                  call mshplt_set_fill_color(kc,0,0,0)
                endif

                if (isurf.ne.0) call mshplt_filled_mesh_3d(4,xp,yp,zp,ioutlined)

                if (imcol.ne.0) then
                  do i=1,4
                    kc=mod(int((zp(i)-zmin)/dz10)+1,8)
                    call mshplt_set_marker_color(kc,0,0,0)
                    call mshplt_marker_3d(1,xp(i),yp(i),zp(i))
                  enddo
                endif
              enddo
            enddo

          else if (phi_ps.le.180.) then

            do ix=1,nx-1
              do iy=ny-1,1,-1

                xp(1)=xmin+(ix-1)*dx
                xp(2)=xp(1)+dx
                xp(3)=xp(2)
                xp(4)=xp(1)

                yp(1)=ymin+(iy-1)*dy
                yp(2)=yp(1)
                yp(3)=yp(2)+dy
                yp(4)=yp(3)

                zp(1)=z(nx*(iy-1)+ix)
                zp(2)=z(nx*(iy-1)+ix+1)
                zp(3)=z(nx*iy+ix+1)
                zp(4)=z(nx*iy+ix)

                if (icoltile.ne.0) then
                  kc=mod(int(((zp(1)+zp(2)+zp(3)+zp(4))/4.-zmin)/dz10)+1,8)
                  call mshplt_set_fill_color(kc,0,0,0)
                endif

                if (isurf.ne.0) call mshplt_filled_mesh_3d(4,xp,yp,zp,ioutlined)

                if (imcol.ne.0) then
                  do i=1,4
                    kc=mod(int((zp(i)-zmin)/dz10)+1,8)
                    call mshplt_set_marker_color(kc,0,0,0)
                    call mshplt_marker_3d(1,xp(i),yp(i),zp(i))
                  enddo
                endif
              enddo
            enddo

          else if (phi_ps.le.270.) then

            do iy=ny-1,1,-1
              do ix=nx-1,1,-1

                xp(1)=xmin+(ix-1)*dx
                xp(2)=xp(1)+dx
                xp(3)=xp(2)
                xp(4)=xp(1)

                yp(1)=ymin+(iy-1)*dy
                yp(2)=yp(1)
                yp(3)=yp(2)+dy
                yp(4)=yp(3)

                zp(1)=z(nx*(iy-1)+ix)
                zp(2)=z(nx*(iy-1)+ix+1)
                zp(3)=z(nx*iy+ix+1)
                zp(4)=z(nx*iy+ix)

                if (icoltile.ne.0) then
                  kc=mod(int(((zp(1)+zp(2)+zp(3)+zp(4))/4.-zmin)/dz10)+1,8)
                  call mshplt_set_fill_color(kc,0,0,0)
                endif

                if (isurf.ne.0) call mshplt_filled_mesh_3d(4,xp,yp,zp,ioutlined)

                if (imcol.ne.0) then
                  do i=1,4
                    kc=mod(int((zp(i)-zmin)/dz10)+1,8)
                    call mshplt_set_marker_color(kc,0,0,0)
                    call mshplt_marker_3d(1,xp(i),yp(i),zp(i))
                  enddo
                endif
              enddo
            enddo

          else if (phi_ps.le.360.) then

            do ix=nx-1,1,-1
              do iy=1,ny-1

                xp(1)=xmin+(ix-1)*dx
                xp(2)=xp(1)+dx
                xp(3)=xp(2)
                xp(4)=xp(1)

                yp(1)=ymin+(iy-1)*dy
                yp(2)=yp(1)
                yp(3)=yp(2)+dy
                yp(4)=yp(3)

                zp(1)=z(nx*(iy-1)+ix)
                zp(2)=z(nx*(iy-1)+ix+1)
                zp(3)=z(nx*iy+ix+1)
                zp(4)=z(nx*iy+ix)

                if (icoltile.ne.0) then
                  kc=mod(int(((zp(1)+zp(2)+zp(3)+zp(4))/4.-zmin)/dz10)+1,8)
                  call mshplt_set_fill_color(kc,0,0,0)
                endif

                if (isurf.ne.0) call mshplt_filled_mesh_3d(4,xp,yp,zp,ioutlined)

                if (imcol.ne.0) then
                  do i=1,4
                    kc=mod(int((zp(i)-zmin)/dz10)+1,8)
                    call mshplt_set_marker_color(kc,0,0,0)
                    call mshplt_marker_3d(1,xp(i),yp(i),zp(i))
                  enddo
                endif
              enddo
            enddo

          endif !phi.le.90.

        endif !theta

      endif

      if (ifbox_ps.ne.0) call mshplt_draw_forground_box

      call mshplt_set_marker_color(mco,mro,mgo,mbo)
      call mshplt_set_fill_color(kco,kro,kgo,kbo)

      call mshplt_fill_buff('% end of mshplt_surf')

      return
      end
*CMZ :  1.02/00 02/10/2014  15.01.51  by  Michael Scheer
*CMZ :  1.01/02 27/09/2014  16.19.04  by  Michael Scheer
*CMZ :  0.01/02 05/09/2014  15.41.43  by  Michael Scheer
*CMZ :  0.00/04 11/08/2014  17.53.18  by  Michael Scheer
*CMZ :  0.00/02 07/07/2014  12.22.09  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_get_box(ibox)

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      integer ibox

      ibox=kbox_ps

      return
      end
*CMZ :  1.03/00 06/10/2014  18.03.21  by  Michael Scheer
*CMZ :  1.02/01 05/10/2014  15.49.28  by  Michael Scheer
*CMZ :  0.01/02 05/09/2014  15.41.43  by  Michael Scheer
*CMZ :  0.00/04 11/08/2014  17.53.18  by  Michael Scheer
*CMZ :  0.00/02 07/07/2014  12.22.09  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_box_3d(xc,yc,zc,xl,yl,zl,
     &  kcol,kred,kgreen,kblue,ioutlined)

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      real x(8),y(8),z(8),xp(4),yp(4),zp(4),
     &  xmin,ymin,zmin,xmax,ymax,zmax,
     &  xc,yc,zc,xl,yl,zl

      integer kcol(3),kred(3),kgreen(3),kblue(3),ioutlined,kc,kr,kg,kb

      ymin=1.0e30
      ymax=-1.0e30

      xmin=xc-xl/2.
      xmax=xc+xl/2.
      ymin=yc-yl/2.
      ymax=yc+yl/2.
      zmin=zc-zl/2.
      zmax=zc+zl/2.

      x(1)=xmin
      x(2)=xmax
      x(3)=xmax
      x(4)=xmin
      x(5:8)=x(1:4)

      y(1)=ymin
      y(2)=ymin
      y(3)=ymax
      y(4)=ymax
      y(5:8)=y(1:4)

      z(1:4)=zmin
      z(5:8)=zmax

      call mshplt_get_fill_color(kc,kr,kg,kb)

      if(theta_ps.le.90.0)then
        call mshplt_set_fill_color(kcol(1),kred(1),kgreen(1),kblue(1))
        call mshplt_filled_mesh_3d(4,x(5),y(5),z(5),ioutlined)
        if (phi_ps.le.90.) then
          xp(1)=x(8)
          xp(2)=x(4)
          xp(3)=x(1)
          xp(4)=x(5)
          yp(1)=y(8)
          yp(2)=y(4)
          yp(3)=y(1)
          yp(4)=y(5)
          zp(1)=z(8)
          zp(2)=z(4)
          zp(3)=z(1)
          zp(4)=z(5)
          call mshplt_set_fill_color(kcol(2),kred(2),kgreen(2),kblue(2))
          call mshplt_filled_mesh_3d(4,xp,yp,zp,ioutlined)
          xp(1)=x(1)
          xp(2)=x(2)
          xp(3)=x(6)
          xp(4)=x(5)
          yp(1)=y(1)
          yp(2)=y(2)
          yp(3)=y(6)
          yp(4)=y(5)
          zp(1)=z(1)
          zp(2)=z(2)
          zp(3)=z(6)
          zp(4)=z(5)
          call mshplt_set_fill_color(kcol(3),kred(3),kgreen(3),kblue(3))
          call mshplt_filled_mesh_3d(4,xp,yp,zp,ioutlined)
        else if (phi_ps.le.180.) then
          xp(1)=x(8)
          xp(2)=x(4)
          xp(3)=x(1)
          xp(4)=x(5)
          yp(1)=y(8)
          yp(2)=y(4)
          yp(3)=y(1)
          yp(4)=y(5)
          zp(1)=z(8)
          zp(2)=z(4)
          zp(3)=z(1)
          zp(4)=z(5)
          call mshplt_set_fill_color(kcol(2),kred(2),kgreen(2),kblue(2))
          call mshplt_filled_mesh_3d(4,xp,yp,zp,ioutlined)
          xp(1)=x(7)
          xp(2)=x(3)
          xp(3)=x(4)
          xp(4)=x(8)
          yp(1)=y(7)
          yp(2)=y(3)
          yp(3)=y(4)
          yp(4)=y(8)
          zp(1)=z(7)
          zp(2)=z(3)
          zp(3)=z(4)
          zp(4)=z(8)
          call mshplt_set_fill_color(kcol(3),kred(3),kgreen(3),kblue(3))
          call mshplt_filled_mesh_3d(4,xp,yp,zp,ioutlined)
        else if (phi_ps.le.270.) then
          xp(1)=x(6)
          xp(2)=x(2)
          xp(3)=x(3)
          xp(4)=x(7)
          yp(1)=y(6)
          yp(2)=y(2)
          yp(3)=y(3)
          yp(4)=y(7)
          zp(1)=z(6)
          zp(2)=z(2)
          zp(3)=z(3)
          zp(4)=z(7)
          call mshplt_set_fill_color(kcol(2),kred(2),kgreen(2),kblue(2))
          call mshplt_filled_mesh_3d(4,xp,yp,zp,ioutlined)
          xp(1)=x(7)
          xp(2)=x(3)
          xp(3)=x(4)
          xp(4)=x(8)
          yp(1)=y(7)
          yp(2)=y(3)
          yp(3)=y(4)
          yp(4)=y(8)
          zp(1)=z(7)
          zp(2)=z(3)
          zp(3)=z(4)
          zp(4)=z(8)
          call mshplt_set_fill_color(kcol(3),kred(3),kgreen(3),kblue(3))
          call mshplt_filled_mesh_3d(4,xp,yp,zp,ioutlined)
        else if (phi_ps.le.360.) then
          xp(1)=x(6)
          xp(2)=x(2)
          xp(3)=x(3)
          xp(4)=x(7)
          yp(1)=y(6)
          yp(2)=y(2)
          yp(3)=y(3)
          yp(4)=y(7)
          zp(1)=z(6)
          zp(2)=z(2)
          zp(3)=z(3)
          zp(4)=z(7)
          call mshplt_set_fill_color(kcol(2),kred(2),kgreen(2),kblue(2))
          call mshplt_filled_mesh_3d(4,xp,yp,zp,ioutlined)
          xp(1)=x(5)
          xp(2)=x(1)
          xp(3)=x(2)
          xp(4)=x(6)
          yp(1)=y(5)
          yp(2)=y(1)
          yp(3)=y(2)
          yp(4)=y(6)
          zp(1)=z(5)
          zp(2)=z(1)
          zp(3)=z(2)
          zp(4)=z(6)
          call mshplt_set_fill_color(kcol(3),kred(3),kgreen(3),kblue(3))
          call mshplt_filled_mesh_3d(4,xp,yp,zp,ioutlined)
        endif
      else !theta
        call mshplt_set_fill_color(kcol(1),kred(1),kgreen(1),kblue(1))
        call mshplt_filled_mesh_3d(4,x,y,z,ioutlined)
        if (phi_ps.le.90.) then
          xp(1)=x(7)
          xp(2)=x(3)
          xp(3)=x(4)
          xp(4)=x(8)
          yp(1)=y(7)
          yp(2)=y(3)
          yp(3)=y(4)
          yp(4)=y(8)
          zp(1)=z(7)
          zp(2)=z(3)
          zp(3)=z(4)
          zp(4)=z(8)
          call mshplt_set_fill_color(kcol(2),kred(2),kgreen(2),kblue(2))
          call mshplt_filled_mesh_3d(4,xp,yp,zp,ioutlined)
          xp(1)=x(6)
          xp(2)=x(2)
          xp(3)=x(3)
          xp(4)=x(7)
          yp(1)=y(6)
          yp(2)=y(2)
          yp(3)=y(3)
          yp(4)=y(7)
          zp(1)=z(6)
          zp(2)=z(2)
          zp(3)=z(3)
          zp(4)=z(7)
          call mshplt_set_fill_color(kcol(3),kred(3),kgreen(3),kblue(3))
          call mshplt_filled_mesh_3d(4,xp,yp,zp,ioutlined)
        else if (phi_ps.le.180.) then
          xp(1)=x(1)
          xp(2)=x(5)
          xp(3)=x(6)
          xp(4)=x(2)
          yp(1)=y(1)
          yp(2)=y(5)
          yp(3)=y(6)
          yp(4)=y(2)
          zp(1)=z(1)
          zp(2)=z(5)
          zp(3)=z(6)
          zp(4)=z(2)
          call mshplt_set_fill_color(kcol(2),kred(2),kgreen(2),kblue(2))
          call mshplt_filled_mesh_3d(4,xp,yp,zp,ioutlined)
          xp(1)=x(6)
          xp(2)=x(2)
          xp(3)=x(3)
          xp(4)=x(7)
          yp(1)=y(6)
          yp(2)=y(2)
          yp(3)=y(3)
          yp(4)=y(7)
          zp(1)=z(6)
          zp(2)=z(2)
          zp(3)=z(3)
          zp(4)=z(7)
          call mshplt_set_fill_color(kcol(3),kred(3),kgreen(3),kblue(3))
          call mshplt_filled_mesh_3d(4,xp,yp,zp,ioutlined)
        else if (phi_ps.le.270.) then
          xp(1)=x(1)
          xp(2)=x(4)
          xp(3)=x(8)
          xp(4)=x(5)
          yp(1)=y(1)
          yp(2)=y(4)
          yp(3)=y(8)
          yp(4)=y(5)
          zp(1)=z(1)
          zp(2)=z(4)
          zp(3)=z(8)
          zp(4)=z(5)
          call mshplt_set_fill_color(kcol(2),kred(2),kgreen(2),kblue(2))
          call mshplt_filled_mesh_3d(4,xp,yp,zp,ioutlined)
          xp(1)=x(1)
          xp(2)=x(2)
          xp(3)=x(6)
          xp(4)=x(5)
          yp(1)=y(1)
          yp(2)=y(2)
          yp(3)=y(6)
          yp(4)=y(5)
          zp(1)=z(1)
          zp(2)=z(2)
          zp(3)=z(6)
          zp(4)=z(5)
          call mshplt_set_fill_color(kcol(3),kred(3),kgreen(3),kblue(3))
          call mshplt_filled_mesh_3d(4,xp,yp,zp,ioutlined)
        else if (phi_ps.le.360.) then
          xp(1)=x(1)
          xp(2)=x(5)
          xp(3)=x(8)
          xp(4)=x(4)
          yp(1)=y(1)
          yp(2)=y(5)
          yp(3)=y(8)
          yp(4)=y(4)
          zp(1)=z(1)
          zp(2)=z(5)
          zp(3)=z(8)
          zp(4)=z(4)
          call mshplt_set_fill_color(kcol(2),kred(2),kgreen(2),kblue(2))
          call mshplt_filled_mesh_3d(4,xp,yp,zp,ioutlined)
          xp(1)=x(4)
          xp(2)=x(8)
          xp(3)=x(7)
          xp(4)=x(3)
          yp(1)=y(4)
          yp(2)=y(8)
          yp(3)=y(7)
          yp(4)=y(3)
          zp(1)=z(4)
          zp(2)=z(8)
          zp(3)=z(7)
          zp(4)=z(3)
          call mshplt_set_fill_color(kcol(3),kred(3),kgreen(3),kblue(3))
          call mshplt_filled_mesh_3d(4,xp,yp,zp,ioutlined)
        endif !phi
      endif !theta

      call mshplt_set_fill_color(kc,kr,kg,kb)

      return
      end
*CMZ :  1.02/01 05/10/2014  15.10.20  by  Michael Scheer
*CMZ :  1.02/00 03/10/2014  12.39.05  by  Michael Scheer
*-- Author :    Michael Scheer   30/09/2014
      subroutine mshplt_draw_forground_box

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      integer ix,iy,kco,kro,kgo,kbo

      call mshplt_get_line_color(kco,kro,kgo,kbo)
      call mshplt_set_line_color(
     &  kframecolor_ps,kframered_ps,kframegreen_ps,kframeblue_ps)

      if (ifbox_ps.eq.0) return

      if (theta_ps.le.90.) then

        if (phi_ps.le.90.) then

          call mshplt_line_raw(
     &      xpcorn_ps(1),ypcorn_ps(1),xpcorn_ps(5),ypcorn_ps(5))
          call mshplt_line_raw(
     &      xpcorn_ps(5),ypcorn_ps(5),xpcorn_ps(8),ypcorn_ps(8))
          call mshplt_line_raw(
     &      xpcorn_ps(5),ypcorn_ps(5),xpcorn_ps(6),ypcorn_ps(6))

        else if (phi_ps.le.180.) then

          call mshplt_line_raw(
     &      xpcorn_ps(7),ypcorn_ps(7),xpcorn_ps(8),ypcorn_ps(8))
          call mshplt_line_raw(
     &      xpcorn_ps(5),ypcorn_ps(5),xpcorn_ps(8),ypcorn_ps(8))
          call mshplt_line_raw(
     &      xpcorn_ps(4),ypcorn_ps(4),xpcorn_ps(8),ypcorn_ps(8))

        else if (phi_ps.le.270.) then

          call mshplt_line_raw(
     &      xpcorn_ps(7),ypcorn_ps(7),xpcorn_ps(8),ypcorn_ps(8))
          call mshplt_line_raw(
     &      xpcorn_ps(3),ypcorn_ps(3),xpcorn_ps(7),ypcorn_ps(7))
          call mshplt_line_raw(
     &      xpcorn_ps(6),ypcorn_ps(6),xpcorn_ps(7),ypcorn_ps(7))

        else if (phi_ps.le.360.) then

          call mshplt_line_raw(
     &      xpcorn_ps(5),ypcorn_ps(5),xpcorn_ps(6),ypcorn_ps(6))
          call mshplt_line_raw(
     &      xpcorn_ps(6),ypcorn_ps(6),xpcorn_ps(7),ypcorn_ps(7))
          call mshplt_line_raw(
     &      xpcorn_ps(2),ypcorn_ps(2),xpcorn_ps(6),ypcorn_ps(6))

        endif !phi.le.90.

      else if (theta_ps.le.180.) then

        if (phi_ps.le.90.) then

          call mshplt_pline_raw(3,xpcorn_ps(2),ypcorn_ps(2))
          call mshplt_line_raw(
     &      xpcorn_ps(3),ypcorn_ps(3),xpcorn_ps(7),ypcorn_ps(7))

        else if (phi_ps.le.180.) then

          ix=1
          iy=3
          call mshplt_pline_raw(3,xpcorn_ps(ix),ypcorn_ps(ix))
          call mshplt_line_raw(
     &      xpcorn_ps(ix+1),ypcorn_ps(ix+1),xpcorn_ps(ix+5),ypcorn_ps(ix+5))

        else if (phi_ps.le.270.) then

          ix=0
          iy=ix+2
          call mshplt_pline_raw_closed(4,xpcorn_ps(1),ypcorn_ps(1))
          call mshplt_line_raw(
     &      xpcorn_ps(ix+1),ypcorn_ps(ix+1),xpcorn_ps(ix+5),ypcorn_ps(ix+5))

        else if (phi_ps.le.360.) then

          ix=3
          iy=3
          call mshplt_pline_raw_closed(4,xpcorn_ps(1),ypcorn_ps(1))
          call mshplt_line_raw(
     &      xpcorn_ps(ix+1),ypcorn_ps(ix+1),xpcorn_ps(ix+5),ypcorn_ps(ix+5))

        endif !phi.le.90.

      endif !theta

      call mshplt_set_line_color(kco,kro,kgo,kbo)

      return
      end
*CMZ :  1.03/01 07/10/2014  14.12.42  by  Michael Scheer
*CMZ :  1.03/00 07/10/2014  10.36.01  by  Michael Scheer
*CMZ :  1.02/01 05/10/2014  15.49.53  by  Michael Scheer
*CMZ :  1.02/00 03/10/2014  12.39.05  by  Michael Scheer
*-- Author :    Michael Scheer   30/09/2014
      subroutine mshplt_lego(nx,xmin,xmax,ny,ymin,ymax,z,ilego)

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      real dx,dy,xmin,xmax,ymin,ymax,z(*),dz10,zmin,zmax,zc,zl,zbase

      integer nx,ny,ix,iy,ioutlined,i,iz,
     &  kco,kro,kbo,kgo,ilego,levels,nseg

      integer kcol(3),kred(3),kgreen(3),kblue(3)

      data ioutlined/1/

      call mshplt_fill_buff('% begin of mshplt_lego')

      call mshplt_get_fill_color(kco,kro,kgo,kbo)

      kcol=kco
      kred=kro
      kgreen=kgreen
      kblue=kbo

      levels=1
      if (ilego.gt.1) levels=ilego

      zmin=1.0e30
      zmax=-1.0e30
      iz=0
      do iy=1,ny
        do ix=1,nx
          iz=iz+1
          if (z(iz).lt.zmin) zmin=z(iz)
          if (z(iz).gt.zmax) zmax=z(iz)
        enddo
      enddo

c obsolete{
c      if(zmin.eq.zmax) then
c        zbase=0.0
c      else if (zmin.ge.0.0.or.zmax.le.0.0) then
c        zbase=zmin
c      else
c        zbase=0.0
c      endif
c
c      dz10=(zmax-zmin)*1.001/max(1,levels)
c obsolete}

      if (log10z_ps.eq.0) then
        zbase=0.0
        dz10=max(abs(zmax),abs(zmin),zmax-zmin)*1.001/max(1,levels)
      else
        zbase=10**int(zmin3d_ps)
        dz10=(zmax/zbase)**(1./max(1,levels))
      endif

      dx=(xmax-xmin)/max(1,(nx-1))
      dy=(ymax-ymin)/max(1,(ny-1))

      if (theta_ps.le.90.) then

        if (phi_ps.le.90.) then

          do iy=ny,1,-1
            do ix=nx,1,-1

              iz=nx*(iy-1)+ix
              if (levels.gt.1) then
                if (log10z_ps.eq.0) then
                  nseg=int(abs(z(iz))/dz10)+1
                else
                  nseg=0
                  do i=1,levels
                    if (zbase*dz10**i.ge.z(iz)) goto 122
                  enddo
 122              nseg=min(levels,i)
                endif
              else
                nseg=1
              endif
              do i=1,nseg
                if (levels.gt.1) then
                  kred=0
                  kgreen=0
                  kblue=0
                  zl=dz10
                  if (log10z_ps.eq.0) then
                    if (z(iz).ge.zbase) then
                      if (zbase+i*dz10.gt.z(iz)) then
                        zl=z(iz)-(zbase+dz10*(i-1))
                      endif
                      zc=zbase+dz10*(i-1)+zl/2.
                    else
                      if (z(iz)+i*dz10.gt.zbase) zl=zbase-(z(iz)+dz10*(i-1))
                      zc=z(iz)+dz10*(i-1)+zl/2.
                    endif
                    kcol=mod(int((zc-zmin)/dz10)+1,8)
                  else
                    zl=zbase*(levels**i-levels**(i-1))
                    zc=zbase*levels**(i-1)+zl/2.
                    zl=zbase*(dz10**i-dz10**(i-1))
                    zc=zbase*dz10**(i-1)+zl/2.
                    if (zc+zl/2.0.gt.z(iz)) then
                      zl=z(iz)-zbase*dz10**(i-1)
                      zc=z(iz)-zl/2.
                    endif
c obsolete                    kcol=mod(int(alog10(zc+zl/2.)-0.001)+1,8)
                    kcol=mod(i+1,8)
                  endif

                  call mshplt_box_3d(
     &              xmin+(ix-1)*dx,ymin+(iy-1)*dy,zc,
     &              dx,dy,zl,kcol,kred,kgreen,kblue,0)
                  if (i.eq.nseg.and.ioutlined.ne.0) then
                    zc=(z(iz)+zbase)/2.
                    zl=abs(z(iz)-zbase)
                    call mshplt_box_3d(
     &                xmin+(ix-1)*dx,ymin+(iy-1)*dy,zc,
     &                dx,dy,zl,kcol,kred,kgreen,kblue,-1)
                  endif
                else
                  zc=(z(iz)+zbase)/2.
                  zl=abs(z(iz)-zbase)
                  call mshplt_box_3d(
     &              xmin+(ix-1)*dx,ymin+(iy-1)*dy,zc,
     &              dx,dy,zl,kcol,kred,kgreen,kblue,1)
                endif
              enddo
c              call mshplt_stop
            enddo
          enddo

        else if (phi_ps.le.180.) then

            do ix=nx,1,-1
              do iy=1,ny

              iz=nx*(iy-1)+ix
              if (levels.gt.1) then
                if (log10z_ps.eq.0) then
                  nseg=int(abs(z(iz))/dz10)+1
                else
                  nseg=0
                  do i=1,levels
                    if (zbase*dz10**i.ge.z(iz)) goto 124
                  enddo
 124              nseg=min(levels,i)
                endif
              else
                nseg=1
              endif
              do i=1,nseg
                if (levels.gt.1) then
                  kred=0
                  kgreen=0
                  kblue=0
                  zl=dz10
                  if (log10z_ps.eq.0) then
                    if (z(iz).ge.zbase) then
                      if (zbase+i*dz10.gt.z(iz)) then
                        zl=z(iz)-(zbase+dz10*(i-1))
                      endif
                      zc=zbase+dz10*(i-1)+zl/2.
                    else
                      if (z(iz)+i*dz10.gt.zbase) zl=zbase-(z(iz)+dz10*(i-1))
                      zc=z(iz)+dz10*(i-1)+zl/2.
                    endif
                    kcol=mod(int((zc-zmin)/dz10)+1,8)
                  else
                    zl=zbase*(levels**i-levels**(i-1))
                    zc=zbase*levels**(i-1)+zl/2.
                    zl=zbase*(dz10**i-dz10**(i-1))
                    zc=zbase*dz10**(i-1)+zl/2.
                    if (zc+zl/2.0.gt.z(iz)) then
                      zl=z(iz)-zbase*dz10**(i-1)
                      zc=z(iz)-zl/2.
                    endif
c obsolete                    kcol=mod(int(alog10(zc+zl/2.)-0.001)+1,8)
                    kcol=mod(i+1,8)
                  endif

                  call mshplt_box_3d(
     &              xmin+(ix-1)*dx,ymin+(iy-1)*dy,zc,
     &              dx,dy,zl,kcol,kred,kgreen,kblue,0)
                  if (i.eq.nseg.and.ioutlined.ne.0) then
                    zc=(z(iz)+zbase)/2.
                    zl=abs(z(iz)-zbase)
                    call mshplt_box_3d(
     &                xmin+(ix-1)*dx,ymin+(iy-1)*dy,zc,
     &                dx,dy,zl,kcol,kred,kgreen,kblue,-1)
                  endif
                else
                  zc=(z(iz)+zbase)/2.
                  zl=abs(z(iz)-zbase)
                  call mshplt_box_3d(
     &              xmin+(ix-1)*dx,ymin+(iy-1)*dy,zc,
     &              dx,dy,zl,kcol,kred,kgreen,kblue,1)
                endif
              enddo
c              call mshplt_stop
              enddo
            enddo

          else if (phi_ps.le.270.) then

            do iy=1,ny
              do ix=1,nx

              iz=nx*(iy-1)+ix
              if (levels.gt.1) then
                if (log10z_ps.eq.0) then
                  nseg=int(abs(z(iz))/dz10)+1
                else
                  nseg=0
                  do i=1,levels
                    if (zbase*dz10**i.ge.z(iz)) goto 125
                  enddo
 125              nseg=min(levels,i)
                endif
              else
                nseg=1
              endif
              do i=1,nseg
                if (levels.gt.1) then
                  kred=0
                  kgreen=0
                  kblue=0
                  zl=dz10
                  if (log10z_ps.eq.0) then
                    if (z(iz).ge.zbase) then
                      if (zbase+i*dz10.gt.z(iz)) then
                        zl=z(iz)-(zbase+dz10*(i-1))
                      endif
                      zc=zbase+dz10*(i-1)+zl/2.
                    else
                      if (z(iz)+i*dz10.gt.zbase) zl=zbase-(z(iz)+dz10*(i-1))
                      zc=z(iz)+dz10*(i-1)+zl/2.
                    endif
                    kcol=mod(int((zc-zmin)/dz10)+1,8)
                  else
                    zl=zbase*(levels**i-levels**(i-1))
                    zc=zbase*levels**(i-1)+zl/2.
                    zl=zbase*(dz10**i-dz10**(i-1))
                    zc=zbase*dz10**(i-1)+zl/2.
                    if (zc+zl/2.0.gt.z(iz)) then
                      zl=z(iz)-zbase*dz10**(i-1)
                      zc=z(iz)-zl/2.
                    endif
c obsolete                    kcol=mod(int(alog10(zc+zl/2.)-0.001)+1,8)
                    kcol=mod(i+1,8)
                  endif

                  call mshplt_box_3d(
     &              xmin+(ix-1)*dx,ymin+(iy-1)*dy,zc,
     &              dx,dy,zl,kcol,kred,kgreen,kblue,0)
                  if (i.eq.nseg.and.ioutlined.ne.0) then
                    zc=(z(iz)+zbase)/2.
                    zl=abs(z(iz)-zbase)
                    call mshplt_box_3d(
     &                xmin+(ix-1)*dx,ymin+(iy-1)*dy,zc,
     &                dx,dy,zl,kcol,kred,kgreen,kblue,-1)
                  endif
                else
                  zc=(z(iz)+zbase)/2.
                  zl=abs(z(iz)-zbase)
                  call mshplt_box_3d(
     &              xmin+(ix-1)*dx,ymin+(iy-1)*dy,zc,
     &              dx,dy,zl,kcol,kred,kgreen,kblue,1)
                endif
              enddo
c              call mshplt_stop
              enddo
            enddo

          else if (phi_ps.le.360.) then
            do ix=1,nx
              do iy=ny,1,-1

              iz=nx*(iy-1)+ix
              if (levels.gt.1) then
                if (log10z_ps.eq.0) then
                  nseg=int(abs(z(iz))/dz10)+1
                else
                  nseg=0
                  do i=1,levels
                    if (zbase*dz10**i.ge.z(iz)) goto 127
                  enddo
 127              nseg=min(levels,i)
                endif
              else
                nseg=1
              endif
              do i=1,nseg
                if (levels.gt.1) then
                  kred=0
                  kgreen=0
                  kblue=0
                  zl=dz10
                  if (log10z_ps.eq.0) then
                    if (z(iz).ge.zbase) then
                      if (zbase+i*dz10.gt.z(iz)) then
                        zl=z(iz)-(zbase+dz10*(i-1))
                      endif
                      zc=zbase+dz10*(i-1)+zl/2.
                    else
                      if (z(iz)+i*dz10.gt.zbase) zl=zbase-(z(iz)+dz10*(i-1))
                      zc=z(iz)+dz10*(i-1)+zl/2.
                    endif
                    kcol=mod(int((zc-zmin)/dz10)+1,8)
                  else
                    zl=zbase*(levels**i-levels**(i-1))
                    zc=zbase*levels**(i-1)+zl/2.
                    zl=zbase*(dz10**i-dz10**(i-1))
                    zc=zbase*dz10**(i-1)+zl/2.
                    if (zc+zl/2.0.gt.z(iz)) then
                      zl=z(iz)-zbase*dz10**(i-1)
                      zc=z(iz)-zl/2.
                    endif
c obsolete                    kcol=mod(int(alog10(zc+zl/2.)-0.001)+1,8)
                    kcol=mod(i+1,8)
                  endif

                  call mshplt_box_3d(
     &              xmin+(ix-1)*dx,ymin+(iy-1)*dy,zc,
     &              dx,dy,zl,kcol,kred,kgreen,kblue,0)
                  if (i.eq.nseg.and.ioutlined.ne.0) then
                    zc=(z(iz)+zbase)/2.
                    zl=abs(z(iz)-zbase)
                    call mshplt_box_3d(
     &                xmin+(ix-1)*dx,ymin+(iy-1)*dy,zc,
     &                dx,dy,zl,kcol,kred,kgreen,kblue,-1)
                  endif
                else
                  zc=(z(iz)+zbase)/2.
                  zl=abs(z(iz)-zbase)
                  call mshplt_box_3d(
     &              xmin+(ix-1)*dx,ymin+(iy-1)*dy,zc,
     &              dx,dy,zl,kcol,kred,kgreen,kblue,1)
                endif
              enddo
c              call mshplt_stop
              enddo
            enddo

          endif !phi.le.90.

        else if (theta_ps.le.180.) then

          if (phi_ps.le.90.) then

            do iy=1,ny
              do ix=1,nx

              iz=nx*(iy-1)+ix
              if (levels.gt.1) then
                if (log10z_ps.eq.0) then
                  nseg=int(abs(z(iz))/dz10)+1
                else
                  nseg=0
                  do i=1,levels
                    if (zbase*dz10**i.ge.z(iz)) goto 150
                  enddo
 150              nseg=min(levels,i)
                endif
              else
                nseg=1
              endif
              do i=nseg,1,-1
                if (levels.gt.1) then
                  kred=0
                  kgreen=0
                  kblue=0
                  zl=dz10
                  if (log10z_ps.eq.0) then
                    if (z(iz).ge.zbase) then
                      if (zbase+i*dz10.gt.z(iz)) then
                        zl=z(iz)-(zbase+dz10*(i-1))
                      endif
                      zc=zbase+dz10*(i-1)+zl/2.
                    else
                      if (z(iz)+i*dz10.gt.zbase) zl=zbase-(z(iz)+dz10*(i-1))
                      zc=z(iz)+dz10*(i-1)+zl/2.
                    endif
                    kcol=mod(int((zc-zmin)/dz10)+1,8)
                  else
                    zl=zbase*(levels**i-levels**(i-1))
                    zc=zbase*levels**(i-1)+zl/2.
                    zl=zbase*(dz10**i-dz10**(i-1))
                    zc=zbase*dz10**(i-1)+zl/2.
                    if (zc+zl/2.0.gt.z(iz)) then
                      zl=z(iz)-zbase*dz10**(i-1)
                      zc=z(iz)-zl/2.
                    endif
c obsolete                    kcol=mod(int(alog10(zc+zl/2.)-0.001)+1,8)
                    kcol=mod(i+1,8)
                  endif

                  call mshplt_box_3d(
     &              xmin+(ix-1)*dx,ymin+(iy-1)*dy,zc,
     &              dx,dy,zl,kcol,kred,kgreen,kblue,0)
                  if (i.eq.1.and.ioutlined.ne.0) then
                    zc=(z(iz)+zbase)/2.
                    zl=abs(z(iz)-zbase)
                    call mshplt_box_3d(
     &                xmin+(ix-1)*dx,ymin+(iy-1)*dy,zc,
     &                dx,dy,zl,kcol,kred,kgreen,kblue,-1)
                  endif
                else
                  zc=(z(iz)+zbase)/2.
                  zl=abs(z(iz)-zbase)
                  call mshplt_box_3d(
     &              xmin+(ix-1)*dx,ymin+(iy-1)*dy,zc,
     &              dx,dy,zl,kcol,kred,kgreen,kblue,1)
                endif
              enddo
c              call mshplt_stop
              enddo
            enddo

          else if (phi_ps.le.180.) then

            do ix=1,nx
              do iy=ny,1,-1

              iz=nx*(iy-1)+ix
              if (levels.gt.1) then
                if (log10z_ps.eq.0) then
                  nseg=int(abs(z(iz))/dz10)+1
                else
                  nseg=0
                  do i=1,levels
                    if (zbase*dz10**i.ge.z(iz)) goto 128
                  enddo
 128              nseg=min(levels,i)
                endif
              else
                nseg=1
              endif
              do i=nseg,1,-1
                if (levels.gt.1) then
                  kred=0
                  kgreen=0
                  kblue=0
                  zl=dz10
                  if (log10z_ps.eq.0) then
                    if (z(iz).ge.zbase) then
                      if (zbase+i*dz10.gt.z(iz)) then
                        zl=z(iz)-(zbase+dz10*(i-1))
                      endif
                      zc=zbase+dz10*(i-1)+zl/2.
                    else
                      if (z(iz)+i*dz10.gt.zbase) zl=zbase-(z(iz)+dz10*(i-1))
                      zc=z(iz)+dz10*(i-1)+zl/2.
                    endif
                    kcol=mod(int((zc-zmin)/dz10)+1,8)
                  else
                    zl=zbase*(levels**i-levels**(i-1))
                    zc=zbase*levels**(i-1)+zl/2.
                    zl=zbase*(dz10**i-dz10**(i-1))
                    zc=zbase*dz10**(i-1)+zl/2.
                    if (zc+zl/2.0.gt.z(iz)) then
                      zl=z(iz)-zbase*dz10**(i-1)
                      zc=z(iz)-zl/2.
                    endif
c obsolete                    kcol=mod(int(alog10(zc+zl/2.)-0.001)+1,8)
                    kcol=mod(i+1,8)
                  endif

                  call mshplt_box_3d(
     &              xmin+(ix-1)*dx,ymin+(iy-1)*dy,zc,
     &              dx,dy,zl,kcol,kred,kgreen,kblue,0)
                  if (i.eq.1.and.ioutlined.ne.0) then
                    zc=(z(iz)+zbase)/2.
                    zl=abs(z(iz)-zbase)
                    call mshplt_box_3d(
     &                xmin+(ix-1)*dx,ymin+(iy-1)*dy,zc,
     &                dx,dy,zl,kcol,kred,kgreen,kblue,-1)
                  endif
                else
                  zc=(z(iz)+zbase)/2.
                  zl=abs(z(iz)-zbase)
                  call mshplt_box_3d(
     &              xmin+(ix-1)*dx,ymin+(iy-1)*dy,zc,
     &              dx,dy,zl,kcol,kred,kgreen,kblue,1)
                endif
              enddo
            enddo
          enddo

        else if (phi_ps.le.270.) then

          do iy=ny,1,-1
            do ix=nx,1,-1

              iz=nx*(iy-1)+ix
              if (levels.gt.1) then
                if (log10z_ps.eq.0) then
                  nseg=int(abs(z(iz))/dz10)+1
                else
                  nseg=0
                  do i=1,levels
                    if (zbase*dz10**i.ge.z(iz)) goto 129
                  enddo
 129              nseg=min(levels,i)
                endif
              else
                nseg=1
              endif
              do i=nseg,1,-1
                if (levels.gt.1) then
                  kred=0
                  kgreen=0
                  kblue=0
                  zl=dz10
                  if (log10z_ps.eq.0) then
                    if (z(iz).ge.zbase) then
                      if (zbase+i*dz10.gt.z(iz)) then
                        zl=z(iz)-(zbase+dz10*(i-1))
                      endif
                      zc=zbase+dz10*(i-1)+zl/2.
                    else
                      if (z(iz)+i*dz10.gt.zbase) zl=zbase-(z(iz)+dz10*(i-1))
                      zc=z(iz)+dz10*(i-1)+zl/2.
                    endif
                    kcol=mod(int((zc-zmin)/dz10)+1,8)
                  else
                    zl=zbase*(levels**i-levels**(i-1))
                    zc=zbase*levels**(i-1)+zl/2.
                    zl=zbase*(dz10**i-dz10**(i-1))
                    zc=zbase*dz10**(i-1)+zl/2.
                    if (zc+zl/2.0.gt.z(iz)) then
                      zl=z(iz)-zbase*dz10**(i-1)
                      zc=z(iz)-zl/2.
                    endif
c obsolete                    kcol=mod(int(alog10(zc+zl/2.)-0.001)+1,8)
                    kcol=mod(i+1,8)
                  endif

                  call mshplt_box_3d(
     &              xmin+(ix-1)*dx,ymin+(iy-1)*dy,zc,
     &              dx,dy,zl,kcol,kred,kgreen,kblue,0)
                  if (i.eq.1.and.ioutlined.ne.0) then
                    zc=(z(iz)+zbase)/2.
                    zl=abs(z(iz)-zbase)
                    call mshplt_box_3d(
     &                xmin+(ix-1)*dx,ymin+(iy-1)*dy,zc,
     &                dx,dy,zl,kcol,kred,kgreen,kblue,-1)
                  endif
                else
                  zc=(z(iz)+zbase)/2.
                  zl=abs(z(iz)-zbase)
                  call mshplt_box_3d(
     &              xmin+(ix-1)*dx,ymin+(iy-1)*dy,zc,
     &              dx,dy,zl,kcol,kred,kgreen,kblue,1)
                endif
              enddo
            enddo
          enddo

        else if (phi_ps.le.360.) then

          do ix=nx,1,-1
            do iy=1,ny

              iz=nx*(iy-1)+ix
              if (levels.gt.1) then
                if (log10z_ps.eq.0) then
                  nseg=int(abs(z(iz))/dz10)+1
                else
                  nseg=0
                  do i=1,levels
                    if (zbase*dz10**i.ge.z(iz)) goto 130
                  enddo
 130              nseg=min(levels,i)
                endif
              else
                nseg=1
              endif
              do i=nseg,1,-1
                if (levels.gt.1) then
                  kred=0
                  kgreen=0
                  kblue=0
                  zl=dz10
                  if (log10z_ps.eq.0) then
                    if (z(iz).ge.zbase) then
                      if (zbase+i*dz10.gt.z(iz)) then
                        zl=z(iz)-(zbase+dz10*(i-1))
                      endif
                      zc=zbase+dz10*(i-1)+zl/2.
                    else
                      if (z(iz)+i*dz10.gt.zbase) zl=zbase-(z(iz)+dz10*(i-1))
                      zc=z(iz)+dz10*(i-1)+zl/2.
                    endif
                    kcol=mod(int((zc-zmin)/dz10)+1,8)
                  else
                    zl=zbase*(levels**i-levels**(i-1))
                    zc=zbase*levels**(i-1)+zl/2.
                    zl=zbase*(dz10**i-dz10**(i-1))
                    zc=zbase*dz10**(i-1)+zl/2.
                    if (zc+zl/2.0.gt.z(iz)) then
                      zl=z(iz)-zbase*dz10**(i-1)
                      zc=z(iz)-zl/2.
                    endif
c obsolete                    kcol=mod(int(alog10(zc+zl/2.)-0.001)+1,8)
                    kcol=mod(i+1,8)
                  endif

                  call mshplt_box_3d(
     &              xmin+(ix-1)*dx,ymin+(iy-1)*dy,zc,
     &              dx,dy,zl,kcol,kred,kgreen,kblue,0)
                  if (i.eq.1.and.ioutlined.ne.0) then
                    zc=(z(iz)+zbase)/2.
                    zl=abs(z(iz)-zbase)
                    call mshplt_box_3d(
     &                xmin+(ix-1)*dx,ymin+(iy-1)*dy,zc,
     &                dx,dy,zl,kcol,kred,kgreen,kblue,-1)
                  endif
                else
                  zc=(z(iz)+zbase)/2.
                  zl=abs(z(iz)-zbase)
                  call mshplt_box_3d(
     &              xmin+(ix-1)*dx,ymin+(iy-1)*dy,zc,
     &              dx,dy,zl,kcol,kred,kgreen,kblue,1)
                endif
              enddo
            enddo
          enddo

        endif !phi.le.90.

      endif !theta

      if (ifbox_ps.ne.0) call mshplt_draw_forground_box

      call mshplt_set_fill_color(kco,kro,kgo,kbo)

      call mshplt_fill_buff('% end of mshplt_lego')

      return
      end
*CMZ :  1.03/01 10/10/2014  13.08.28  by  Michael Scheer
*CMZ :  0.01/03 23/09/2014  12.00.19  by  Michael Scheer
*CMZ :  0.01/02 05/09/2014  15.41.43  by  Michael Scheer
*CMZ :  0.00/04 14/08/2014  12.58.41  by  Michael Scheer
*CMZ :  0.00/03 06/08/2014  09.32.19  by  Michael Scheer
*CMZ :  0.00/02 09/07/2014  14.48.35  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_reset_clipping

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      nclips_ps=0

      call mshplt_fill_buff('initclip') !reset clipping

      return
      end
*CMZ :          25/09/2016  11.32.32  by  Michael Scheer
*CMZ :  1.03/01 08/10/2014  14.16.19  by  Michael Scheer
*CMZ :  1.03/00 06/10/2014  16.06.01  by  Michael Scheer
*CMZ :  1.02/01 05/10/2014  10.00.08  by  Michael Scheer
*CMZ :  1.02/00 03/10/2014  11.12.46  by  Michael Scheer
*CMZ :  1.01/02 26/09/2014  13.31.57  by  Michael Scheer
*CMZ :  1.01/00 25/09/2014  08.49.12  by  Michael Scheer
*CMZ :  1.00/01 24/09/2014  13.15.20  by  Michael Scheer
*CMZ :  0.01/03 23/09/2014  13.56.57  by  Michael Scheer
*CMZ :  0.01/02 05/09/2014  15.41.43  by  Michael Scheer
*CMZ :  0.00/06 21/08/2014  16.57.22  by  Michael Scheer
*CMZ :  0.00/04 14/08/2014  12.58.41  by  Michael Scheer
*CMZ :  0.00/03 06/08/2014  11.22.28  by  Michael Scheer
*CMZ :  0.00/02 09/07/2014  16.03.13  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_frame3d_xzy(xminin,xmaxin,zminin,zmaxin,yminin,ymaxin,
     &  xtit,ztit,ytit,chopt)

c ONYL FOR SPECIAL PURPOSES. TO GET REAL RESULT,
C YOU MUST PLOT YMAX-(Y-YMIN) INSTEAD OF Y DATA INTO THIS FRAME

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.


      real xmin,xmax,ymin,ymax,zmin,zmax,
     &  xminin,xmaxin,yminin,ymaxin,zminin,zmaxin,
     &  xcorn(8),ycorn(8),zcorn(8),ang,xyplen,
     &  ticheight,chhe,titang,rellabang

      integer ilinestyleo,ilinecoloro,iaxis,ibox,ic,iticside,ilabside
      integer kred,kgreen,kblue,kcolor,ifirst
      integer kredo,kgreeno,kblueo

      character(*) xtit,ytit,ztit,chopt
      character(2048) xtitd,ytitd,ztitd,choptd

      data xcorn/-0.5,0.5,0.5,-0.5,-0.5,0.5,0.5,-0.5/
      data ycorn/-0.5,-0.5,0.5,0.5,-0.5,-0.5,0.5,0.5/
      data zcorn/-0.5,-0.5,-0.5,-0.5,0.5,0.5,0.5,0.5/

      xmax=xmaxin
      xmin=xminin
      ymax=ymaxin
      ymin=yminin
      zmax=zmaxin
      zmin=zminin

      if (log10x_ps.ne.0) then

        if (xminin.le.0.0) then
          xmin=-30.
        else
          xmin=float(int(alog10(xminin)))
          if (xmin.le.0.0.and.
     &      (10.**xmin-xminin)/10.**xmin.gt.1.0e-5) xmin=xmin-1.
        endif

        if (xmaxin.le.xminin) then
          xmax=xmin+1
        else
          xmax=float(int(alog10(xmax)))
          if (xmax.ge.0.0.and.
     &      (xmaxin-10.**xmax)/10.**xmax.gt.1.0e-5)
     &      xmax=xmax+1.
        endif
      endif

      if (log10y_ps.ne.0) then

        if (yminin.le.0.0) then
          ymin=-30.
        else
          ymin=float(int(alog10(yminin)))
          if (ymin.le.0.0.and.
     &      (10.**ymin-yminin)/10.**ymin.gt.1.0e-5) ymin=ymin-1.
        endif

        if (ymaxin.le.yminin) then
          ymax=ymin+1
        else
          ymax=float(int(alog10(ymax)))
          if (ymax.ge.0.0.and.
     &      (ymaxin-10.**ymax)/10.**ymax.gt.1.0e-5)
     &      ymax=ymax+1.
        endif
      endif

      if (log10z_ps.ne.0) then

        if (zminin.le.0.0) then
          zmin=-30.
        else
          zmin=float(int(alog10(zminin)))
          if (zmin.le.0.0.and.
     &      (10.**zmin-zminin)/10.**zmin.gt.1.0e-5) zmin=zmin-1.
        endif

        if (zmaxin.le.zminin) then
          zmax=zmin+1
        else
          zmax=float(int(alog10(zmax)))
          if (zmax.ge.0.0.and.
     &      (zmaxin-10.**zmax)/10.**zmax.gt.1.0e-5)
     &      zmax=zmax+1.
        endif
      endif

      ticheight=ticsiz_ps
      chhe=chhe_ps

      write(lun_ps,'(a)')'% begin of mshplt_frame3d'

      if (ihigzmode_ps.ne.0.and.isameframe_ps.eq.0.and.kzone_ps.gt.0) then
        if (kzone_ps.lt.nzone_ps) then
          ifirst=kzone_ps+1
          call mshplt_zone(nxzone_ps,nyzone_ps,ifirst,'s')
          else
            if (itouched_ps.ne.0) call mshplt_newpage
          call mshplt_zone(nxzone_ps,nyzone_ps,1,' ')
        endif
      endif

      kcolor=kFramecolor_ps
      kgreen=kFramegreen_ps
      kred=kFramered_ps
      kblue=kFrameblue_ps

      call mshplt_set_frame_color(kcolor,kgreen,kred,kblue)

      if (nzone_ps.le.0) then
        xcorn(1)=-0.5
        xcorn(2)=0.5
        xcorn(3)=0.5
        xcorn(4)=-0.5
        xcorn(5)=-0.5
        xcorn(6)=0.5
        xcorn(7)=0.5
        xcorn(8)=-0.5
        ycorn(1)=-0.5
        ycorn(2)=-0.5
        ycorn(3)=0.5
        ycorn(4)=0.5
        ycorn(5)=-0.5
        ycorn(6)=-0.5
        ycorn(7)=0.5
        ycorn(8)=0.5
        zcorn(1)=-0.5
        zcorn(2)=-0.5
        zcorn(3)=-0.5
        zcorn(4)=-0.5
        zcorn(5)=0.5
        zcorn(6)=0.5
        zcorn(7)=0.5
        zcorn(8)=0.5
      else
        xcorn(1)=-0.625
        xcorn(2)=0.625
        xcorn(3)=0.625
        xcorn(4)=-0.625
        xcorn(5)=-0.625
        xcorn(6)=0.625
        xcorn(7)=0.625
        xcorn(8)=-0.625
        ycorn(1)=-0.625
        ycorn(2)=-0.625
        ycorn(3)=0.625
        ycorn(4)=0.625
        ycorn(5)=-0.625
        ycorn(6)=-0.625
        ycorn(7)=0.625
        ycorn(8)=0.625
        zcorn(1)=-0.625
        zcorn(2)=-0.625
        zcorn(3)=-0.625
        zcorn(4)=-0.625
        zcorn(5)=0.625
        zcorn(6)=0.625
        zcorn(7)=0.625
        zcorn(8)=0.625
      endif

      xtitd=xtit
      ytitd=ytit
      ztitd=ztit
      choptd=chopt

      if(chopt.eq.'') choptd='AB'

      iaxis=0
      ibox=0
      ifbox_ps=0
      do ic=1,len_trim(choptd)
        if (choptd(ic:ic).eq.'A'.or.choptd(ic:ic).eq.'a') iaxis=1
        if (choptd(ic:ic).eq.'B'.or.choptd(ic:ic).eq.'b') then
          ibox=1
          ifbox_ps=1
        endif
      enddo

      ilinestyleo=ilinestyle_ps
      call mshplt_get_line_color(ilinecoloro,kredo,kblueo,kgreeno)

      call mshplt_3dto2d(8,xcorn,ycorn,zcorn,xpcorn_ps,ypcorn_ps)

      wxmin_ps=min(xpcorn_ps(1),xpcorn_ps(2),xpcorn_ps(3),xpcorn_ps(4),
     &  xpcorn_ps(5),xpcorn_ps(6),xpcorn_ps(7),xpcorn_ps(8))
      wxmax_ps=max(xpcorn_ps(1),xpcorn_ps(2),xpcorn_ps(3),xpcorn_ps(4),
     &  xpcorn_ps(5),xpcorn_ps(6),xpcorn_ps(7),xpcorn_ps(8))

      wymin_ps=min(ypcorn_ps(1),ypcorn_ps(2),ypcorn_ps(3),ypcorn_ps(4),
     &  ypcorn_ps(5),ypcorn_ps(6),ypcorn_ps(7),ypcorn_ps(8))
      wymax_ps=max(ypcorn_ps(1),ypcorn_ps(2),ypcorn_ps(3),ypcorn_ps(4),
     &  ypcorn_ps(5),ypcorn_ps(6),ypcorn_ps(7),ypcorn_ps(8))

      scalex_ps=xsiz_ps/(wxmax_ps-wxmin_ps)
      scaley_ps=ysiz_ps/(wymax_ps-wymin_ps)

      if (log10x_ps.eq.0) then
        xmin3d_ps=xmin
        xmax3d_ps=xmax
      else
        if (xmin.le.0.0) xmin=1.0e-30
        if (xmax.le.xmin) xmax=1.0e30
        xmin3d_ps=alog10(xmin)
        xmax3d_ps=alog10(xmax)
      endif

      if (log10y_ps.eq.0) then
        ymin3d_ps=ymin
        ymax3d_ps=ymax
      else
        if (ymin.le.0.0) ymin=1.0e-30
        if (ymax.le.ymin) ymax=1.0e30
        ymin3d_ps=alog10(ymin)
        ymax3d_ps=alog10(ymax)
      endif

      if (log10z_ps.eq.0) then
        zmin3d_ps=zmin
        zmax3d_ps=zmax
      else
        if (zmin.le.0.0) zmin=1.0e-30
        if (zmax.le.zmin) zmax=1.0e30
        zmin3d_ps=alog10(zmin)
        zmax3d_ps=alog10(zmax)
      endif

      call mshplt_set_line_style(1)

      if (iaxis.ne.1.or.theta_ps.gt.90.or.phi_ps.gt.90.) then
        print*,"*** Warning in mshplt_frame3d_xzy: This routine is not ready"
        print*,"for iaxis=0 or theta>90 or phi>90, sorry!"
      endif

      if (iaxis.eq.1) then

        if (theta_ps.le.90.0) then
          if (phi_ps.le.90.0) then

            if (log10x_ps.eq.0) then
              call mshplt_view_axis(
     &          xpcorn_ps(1),xpcorn_ps(2),ypcorn_ps(1),ypcorn_ps(2),
     &          xmin,xmax,xtitd(1:len_trim(xtitd)),-1,0.,0.,
     &          xtitoff_ps,
     &          xlaboff_ps) !x-axis
            else
              xyplen=sqrt(
     &          (xpcorn_ps(2)-xpcorn_ps(1))**2+
     &          (ypcorn_ps(2)-ypcorn_ps(1))**2
     &          )
              ang=atan2(
     &          (ypcorn_ps(2)-ypcorn_ps(1)),
     &          (xpcorn_ps(2)-xpcorn_ps(1))
     &          )*radtodeg_ps
              iticside=-1
              ilabside=1
              titang=0.
              rellabang=-ang
              call mshplt_view_log_axis(
     &          xpcorn_ps(1),ypcorn_ps(1),xyplen,
     &          xminin,xmaxin,ang,
     &          chhe,ticheight,xtit,
     &          iticside,ilabside,
     &          titang,ytitoff_ps,
     &          rellabang,xlaboff_ps) !x-axis
            endif

            if (log10y_ps.eq.0) then
              call mshplt_view_axis(
     &          xpcorn_ps(4),xpcorn_ps(1),ypcorn_ps(4),ypcorn_ps(1),
     &          ymin,ymax,ytitd(1:len_trim(ytitd)),-1,0.,0.,
     &          ytitoff_ps,
     &          ylaboff_ps) !y-axis
            else
              xyplen=sqrt(
     &          (xpcorn_ps(4)-xpcorn_ps(1))**2+
     &          (ypcorn_ps(4)-ypcorn_ps(1))**2
     &          )
              ang=atan2(
     &          (ypcorn_ps(1)-ypcorn_ps(4)),
     &          (xpcorn_ps(1)-xpcorn_ps(4))
     &          )*radtodeg_ps
              iticside=+1
              ilabside=1
              titang=0.
              rellabang=-ang
              call mshplt_view_log_axis(
     &          xpcorn_ps(4),ypcorn_ps(4),xyplen,
     &          yminin,ymaxin,ang,
     &          chhe,ticheight,ytit,
     &          iticside,ilabside,
     &          titang,ztitoff_ps,
     &          rellabang,2.*zlaboff_ps) !y-axis
            endif

            if (log10z_ps.eq.0) then
              call mshplt_view_axis(
     &          xpcorn_ps(4),xpcorn_ps(8),ypcorn_ps(4),ypcorn_ps(8),
     &          zmin,zmax,ztitd(1:len_trim(ztitd)),1,-90.,0.,
     &          -ztitoff_ps,-zlaboff_ps) !z-axis
            else
              xyplen=sqrt(
     &          (xpcorn_ps(8)-xpcorn_ps(4))**2+
     &          (ypcorn_ps(8)-ypcorn_ps(4))**2
     &          )
              ang=atan2(
     &          (ypcorn_ps(8)-ypcorn_ps(4)),
     &          (xpcorn_ps(8)-xpcorn_ps(4))
     &          )*radtodeg_ps
              iticside=1
              ilabside=-1
              titang=0.
              rellabang=-90.
              call mshplt_view_log_axis(
     &          xpcorn_ps(4),ypcorn_ps(4),xyplen,
     &          zminin,zmaxin,ang,
     &          chhe,ticheight,ztit,
     &          iticside,ilabside,
     &          titang,ztitoff_ps,
     &          rellabang,-zlaboff_ps*1.75) !z-axis
            endif

          else if (phi_ps.le.180.0) then

            if (log10x_ps.eq.0) then
              call mshplt_view_axis(
     &          xpcorn_ps(4),xpcorn_ps(3),ypcorn_ps(4),ypcorn_ps(3),
     &          xmin,xmax,xtitd(1:len_trim(xtitd)),1,180.,180.,
     &          xtitoff_ps,
     &          -xlaboff_ps) !x-axis
            else
              xyplen=sqrt(
     &          (xpcorn_ps(3)-xpcorn_ps(4))**2+
     &          (ypcorn_ps(3)-ypcorn_ps(4))**2
     &          )
              ang=atan2(
     &          (ypcorn_ps(3)-ypcorn_ps(4)),
     &          (xpcorn_ps(3)-xpcorn_ps(4))
     &          )*radtodeg_ps
              iticside=1
              ilabside=-1
              titang=180.
              rellabang=-ang
              call mshplt_view_log_axis(
     &          xpcorn_ps(4),ypcorn_ps(4),xyplen,
     &          xminin,xmaxin,ang,
     &          chhe,ticheight,xtit,
     &          iticside,ilabside,
     &          titang,ztitoff_ps,
     &          rellabang,-2.*xlaboff_ps) !x-axis
            endif

            if (log10y_ps.eq.0) then
              call mshplt_view_axis(
     &          xpcorn_ps(1),xpcorn_ps(4),ypcorn_ps(1),ypcorn_ps(4),
     &          ymin,ymax,ytitd(1:len_trim(ytitd)),1,180.,180.,
     &          ytitoff_ps,
     &          -ylaboff_ps) !y-axis
            else
              xyplen=sqrt(
     &          (xpcorn_ps(4)-xpcorn_ps(1))**2+
     &          (ypcorn_ps(4)-ypcorn_ps(1))**2
     &          )
              ang=atan2(
     &          (ypcorn_ps(4)-ypcorn_ps(1)),
     &          (xpcorn_ps(4)-xpcorn_ps(1))
     &          )*radtodeg_ps
              iticside=1
              ilabside=1
              titang=180.
              rellabang=-ang
              call mshplt_view_log_axis(
     &          xpcorn_ps(1),ypcorn_ps(1),xyplen,
     &          yminin,ymaxin,ang,
     &          chhe,ticheight,ytit,
     &          iticside,ilabside,
     &          titang,ztitoff_ps,
     &          rellabang,-1.5*zlaboff_ps) !y-axis
            endif

            if (log10z_ps.eq.0) then
              call mshplt_view_axis(
     &          xpcorn_ps(3),xpcorn_ps(7),ypcorn_ps(3),ypcorn_ps(7),
     &          zmin,zmax,ztitd(1:len_trim(ztitd)),1,-90.,0.,
     &          -ztitoff_ps,-zlaboff_ps) !z-axis
            else
              xyplen=sqrt(
     &          (xpcorn_ps(7)-xpcorn_ps(3))**2+
     &          (ypcorn_ps(7)-ypcorn_ps(3))**2
     &          )
              ang=atan2(
     &          (ypcorn_ps(7)-ypcorn_ps(3)),
     &          (xpcorn_ps(7)-xpcorn_ps(3))
     &          )*radtodeg_ps
              iticside=1
              ilabside=-1
              titang=0.
              rellabang=-90.
              call mshplt_view_log_axis(
     &          xpcorn_ps(3),ypcorn_ps(3),xyplen,
     &          zminin,zmaxin,ang,
     &          chhe,ticheight,ztit,
     &          iticside,ilabside,
     &          titang,ztitoff_ps,
     &          rellabang,zlaboff_ps*0.6) !z-axis
            endif

          else if (phi_ps.le.270.0) then

            if (log10x_ps.eq.0) then
              call mshplt_view_axis(
     &          xpcorn_ps(4),xpcorn_ps(3),ypcorn_ps(4),ypcorn_ps(3),
     &          xmin,xmax,xtitd(1:len_trim(xtitd)),1,180.,180.,
     &          xtitoff_ps,
     &          -xlaboff_ps) !x-axis
            else
              xyplen=sqrt(
     &          (xpcorn_ps(3)-xpcorn_ps(4))**2+
     &          (ypcorn_ps(3)-ypcorn_ps(4))**2
     &          )
              ang=atan2(
     &          (ypcorn_ps(3)-ypcorn_ps(4)),
     &          (xpcorn_ps(3)-xpcorn_ps(4))
     &          )*radtodeg_ps
              iticside=1
              ilabside=-1
              titang=180.
              rellabang=-ang
              call mshplt_view_log_axis(
     &          xpcorn_ps(4),ypcorn_ps(4),xyplen,
     &          xminin,xmaxin,ang,
     &          chhe,ticheight,xtit,
     &          iticside,ilabside,
     &          titang,ztitoff_ps*0.7,
     &          rellabang,-2.*xlaboff_ps) !x-axis
            endif

            if (log10y_ps.eq.0) then
              call mshplt_view_axis(
     &          xpcorn_ps(2),xpcorn_ps(3),ypcorn_ps(2),ypcorn_ps(3),
     &          ymin,ymax,ytitd(1:len_trim(ytitd)),-1,0.,0.,
     &          1.1*ytitoff_ps,
     &          1.25*ylaboff_ps) !y-axis
            else
              xyplen=sqrt(
     &          (xpcorn_ps(3)-xpcorn_ps(2))**2+
     &          (ypcorn_ps(3)-ypcorn_ps(2))**2
     &          )
              ang=atan2(
     &          (ypcorn_ps(3)-ypcorn_ps(2)),
     &          (xpcorn_ps(3)-xpcorn_ps(2))
     &          )*radtodeg_ps
              iticside=-1
              ilabside=-1
              titang=0.
              rellabang=-ang
              call mshplt_view_log_axis(
     &          xpcorn_ps(2),ypcorn_ps(2),xyplen,
     &          yminin,ymaxin,ang,
     &          chhe,ticheight,ytit,
     &          iticside,ilabside,
     &          titang,ztitoff_ps,
     &          rellabang,1.25*zlaboff_ps) !y-axis
            endif

            if (log10z_ps.eq.0) then
              call mshplt_view_axis(
     &          xpcorn_ps(2),xpcorn_ps(6),ypcorn_ps(2),ypcorn_ps(6),
     &          zmin,zmax,ztitd(1:len_trim(ztitd)),1,-90.,0.,
     &          -ztitoff_ps,-zlaboff_ps) !z-axis
            else
              xyplen=sqrt(
     &          (xpcorn_ps(6)-xpcorn_ps(2))**2+
     &          (ypcorn_ps(6)-ypcorn_ps(2))**2
     &          )
              ang=atan2(
     &          (ypcorn_ps(6)-ypcorn_ps(2)),
     &          (xpcorn_ps(6)-xpcorn_ps(2))
     &          )*radtodeg_ps
              iticside=1
              ilabside=-1
              titang=0.
              rellabang=-90.
              call mshplt_view_log_axis(
     &          xpcorn_ps(2),ypcorn_ps(2),xyplen,
     &          zminin,zmaxin,ang,
     &          chhe,ticheight,ztit,
     &          iticside,ilabside,
     &          titang,ztitoff_ps,
     &          rellabang,-3.*zlaboff_ps*0.6) !z-axis
            endif

          else if (phi_ps.le.360.0) then

            if (log10x_ps.eq.0) then
              call mshplt_view_axis(
     &          xpcorn_ps(1),xpcorn_ps(2),ypcorn_ps(1),ypcorn_ps(2),
     &          xmin,xmax,xtitd(1:len_trim(xtitd)),-1,0.,0.,
     &          1.1*xtitoff_ps,
     &          1.2*xlaboff_ps) !x-axis
            else
              xyplen=sqrt(
     &          (xpcorn_ps(1)-xpcorn_ps(2))**2+
     &          (ypcorn_ps(1)-ypcorn_ps(2))**2
     &          )
              ang=atan2(
     &          (ypcorn_ps(2)-ypcorn_ps(1)),
     &          (xpcorn_ps(2)-xpcorn_ps(1))
     &          )*radtodeg_ps
              iticside=-1
              ilabside=-1
              titang=0.
              rellabang=-ang
              call mshplt_view_log_axis(
     &          xpcorn_ps(1),ypcorn_ps(1),xyplen,
     &          xminin,xmaxin,ang,
     &          chhe,ticheight,xtit,
     &          iticside,ilabside,
     &          titang,ztitoff_ps*0.7,
     &          rellabang,xlaboff_ps) !x-axis
            endif

            if (log10y_ps.eq.0) then
              call mshplt_view_axis(
     &          xpcorn_ps(2),xpcorn_ps(3),ypcorn_ps(2),ypcorn_ps(3),
     &          ymin,ymax,ytitd(1:len_trim(ytitd)),-1,0.,0.,
     &          1.2*ytitoff_ps,
     &          1.1*ylaboff_ps) !y-axis
            else
              xyplen=sqrt(
     &          (xpcorn_ps(3)-xpcorn_ps(2))**2+
     &          (ypcorn_ps(3)-ypcorn_ps(2))**2
     &          )
              ang=atan2(
     &          (ypcorn_ps(3)-ypcorn_ps(2)),
     &          (xpcorn_ps(3)-xpcorn_ps(2))
     &          )*radtodeg_ps
              iticside=-1
              ilabside=1
              titang=0.
              rellabang=-ang
              call mshplt_view_log_axis(
     &          xpcorn_ps(2),ypcorn_ps(2),xyplen,
     &          yminin,ymaxin,ang,
     &          chhe,ticheight,ytit,
     &          iticside,ilabside,
     &          titang,ztitoff_ps*0.8,
     &          rellabang,zlaboff_ps) !y-axis
            endif

            if (log10z_ps.eq.0) then
              call mshplt_view_axis(
     &          xpcorn_ps(1),xpcorn_ps(5),ypcorn_ps(1),ypcorn_ps(5),
     &          zmin,zmax,ztitd(1:len_trim(ztitd)),1,-90.,0.,
     &          -ztitoff_ps,-zlaboff_ps) !z-axis
            else
              xyplen=sqrt(
     &          (xpcorn_ps(5)-xpcorn_ps(1))**2+
     &          (ypcorn_ps(5)-ypcorn_ps(1))**2
     &          )
              ang=atan2(
     &          (ypcorn_ps(5)-ypcorn_ps(1)),
     &          (xpcorn_ps(5)-xpcorn_ps(1))
     &          )*radtodeg_ps
              iticside=1
              ilabside=-1
              titang=0.
              rellabang=-90.
              call mshplt_view_log_axis(
     &          xpcorn_ps(1),ypcorn_ps(1),xyplen,
     &          zminin,zmaxin,ang,
     &          chhe,ticheight,ztit,
     &          iticside,ilabside,
     &          titang,ztitoff_ps,
     &          rellabang,-3.*zlaboff_ps*0.6) !z-axis
            endif

          endif !if (phi_ps.le.90.0) then

        else if (theta_ps.le.180.0) then

          if (phi_ps.le.90.0) then

            if (log10x_ps.eq.0) then
              call mshplt_view_axis(
     &          xpcorn_ps(8),xpcorn_ps(7),ypcorn_ps(8),ypcorn_ps(7),
     &          xmin,xmax,xtitd(1:len_trim(xtitd)),-1,0.,0.,
     &          xtitoff_ps,
     &          xlaboff_ps) !x-axis
            else
              xyplen=sqrt(
     &          (xpcorn_ps(7)-xpcorn_ps(8))**2+
     &          (ypcorn_ps(7)-ypcorn_ps(8))**2
     &          )
              ang=atan2(
     &          (ypcorn_ps(7)-ypcorn_ps(8)),
     &          (xpcorn_ps(7)-xpcorn_ps(8))
     &          )*radtodeg_ps
              iticside=-1
              ilabside=1
              titang=0.
              rellabang=-ang
              call mshplt_view_log_axis(
     &          xpcorn_ps(8),ypcorn_ps(8),xyplen,
     &          xminin,xmaxin,ang,
     &          chhe,ticheight,xtit,
     &          iticside,ilabside,
     &          titang,ytitoff_ps,
     &          rellabang,xlaboff_ps) !x-axis
            endif

            if (log10y_ps.eq.0) then
              call mshplt_view_axis(
     &          xpcorn_ps(6),xpcorn_ps(7),ypcorn_ps(6),ypcorn_ps(7),
     &          ymin,ymax,ytitd(1:len_trim(ytitd)),1,180.,180.,
     &          ytitoff_ps,
     &          -ylaboff_ps) !y-axis
            else
              xyplen=sqrt(
     &          (xpcorn_ps(7)-xpcorn_ps(6))**2+
     &          (ypcorn_ps(7)-ypcorn_ps(6))**2
     &          )
              ang=atan2(
     &          (ypcorn_ps(7)-ypcorn_ps(6)),
     &          (xpcorn_ps(7)-xpcorn_ps(6))
     &          )*radtodeg_ps
              iticside=1
              ilabside=1
              titang=180.
              rellabang=-ang
              call mshplt_view_log_axis(
     &          xpcorn_ps(6),ypcorn_ps(6),xyplen,
     &          yminin,ymaxin,ang,
     &          chhe,ticheight,ytit,
     &          iticside,ilabside,
     &          titang,ztitoff_ps,
     &          rellabang,-1.5*zlaboff_ps) !y-axis
            endif

            if (log10z_ps.eq.0) then
              call mshplt_view_axis(
     &          xpcorn_ps(4),xpcorn_ps(8),ypcorn_ps(4),ypcorn_ps(8),
     &          zmin,zmax,ztitd(1:len_trim(ztitd)),1,+90.,0.,
     &          ztitoff_ps,zlaboff_ps) !z-axis
            else
              xyplen=sqrt(
     &          (xpcorn_ps(8)-xpcorn_ps(4))**2+
     &          (ypcorn_ps(8)-ypcorn_ps(4))**2
     &          )
              ang=atan2(
     &          (ypcorn_ps(8)-ypcorn_ps(4)),
     &          (xpcorn_ps(8)-xpcorn_ps(4))
     &          )*radtodeg_ps
              iticside=-1
              ilabside=-1
              titang=0.
              rellabang=-ang
              call mshplt_view_log_axis(
     &          xpcorn_ps(4),ypcorn_ps(4),xyplen,
     &          zminin,zmaxin,ang,
     &          chhe,ticheight,ztit,
     &          iticside,ilabside,
     &          titang,ztitoff_ps,
     &          rellabang,-1.75*zlaboff_ps) !z-axis
            endif

          else if (phi_ps.le.180.0) then

            if (log10x_ps.eq.0) then
              call mshplt_view_axis(
     &          xpcorn_ps(5),xpcorn_ps(6),ypcorn_ps(5),ypcorn_ps(6),
     &          xmin,xmax,xtitd(1:len_trim(xtitd)),1,180.,180.,
     &          1.1*xtitoff_ps,
     &          -xlaboff_ps) !x-axis
            else
              xyplen=sqrt(
     &          (xpcorn_ps(6)-xpcorn_ps(5))**2+
     &          (ypcorn_ps(6)-ypcorn_ps(5))**2
     &          )
              ang=atan2(
     &          (ypcorn_ps(6)-ypcorn_ps(5)),
     &          (xpcorn_ps(6)-xpcorn_ps(5))
     &          )*radtodeg_ps
              iticside=1
              ilabside=1
              titang=180.
              rellabang=-ang
              call mshplt_view_log_axis(
     &          xpcorn_ps(5),ypcorn_ps(5),xyplen,
     &          xminin,xmaxin,ang,
     &          chhe,ticheight,xtit,
     &          iticside,ilabside,
     &          titang,ztitoff_ps,
     &          rellabang,-1.5*xlaboff_ps) !x-axis
            endif

            if (log10y_ps.eq.0) then
              call mshplt_view_axis(
     &          xpcorn_ps(6),xpcorn_ps(7),ypcorn_ps(6),ypcorn_ps(7),
     &          ymin,ymax,ytitd(1:len_trim(ytitd)),1,180.,180.,
     &          ytitoff_ps,
     &          -ylaboff_ps) !y-axis
            else
              xyplen=sqrt(
     &          (xpcorn_ps(7)-xpcorn_ps(6))**2+
     &          (ypcorn_ps(7)-ypcorn_ps(6))**2
     &          )
              ang=atan2(
     &          (ypcorn_ps(7)-ypcorn_ps(6)),
     &          (xpcorn_ps(7)-xpcorn_ps(6))
     &          )*radtodeg_ps
              iticside=1
              ilabside=1
              titang=180.
              rellabang=-ang
              call mshplt_view_log_axis(
     &          xpcorn_ps(6),ypcorn_ps(6),xyplen,
     &          yminin,ymaxin,ang,
     &          chhe,ticheight,ytit,
     &          iticside,ilabside,
     &          titang,ztitoff_ps,
     &          rellabang,-2.*zlaboff_ps) !y-axis
            endif

            if (log10z_ps.eq.0) then
              call mshplt_view_axis(
     &          xpcorn_ps(3),xpcorn_ps(7),ypcorn_ps(3),ypcorn_ps(7),
     &          zmin,zmax,ztitd(1:len_trim(ztitd)),-1,+90.,180.,
     &          -ztitoff_ps,zlaboff_ps) !z-axis
            else
              xyplen=sqrt(
     &          (xpcorn_ps(7)-xpcorn_ps(3))**2+
     &          (ypcorn_ps(7)-ypcorn_ps(3))**2
     &          )
              ang=atan2(
     &          (ypcorn_ps(7)-ypcorn_ps(3)),
     &          (xpcorn_ps(7)-xpcorn_ps(3))
     &          )*radtodeg_ps
              iticside=-1
              ilabside=-1
              titang=0.
              rellabang=+90.
              call mshplt_view_log_axis(
     &          xpcorn_ps(3),ypcorn_ps(3),xyplen,
     &          zminin,zmaxin,ang,
     &          chhe,ticheight,ztit,
     &          iticside,ilabside,
     &          titang,ztitoff_ps,
     &          rellabang,-1.75*zlaboff_ps) !z-axis
            endif

          else if (phi_ps.le.270.0) then

            if (log10x_ps.eq.0) then
              call mshplt_view_axis(
     &          xpcorn_ps(5),xpcorn_ps(6),ypcorn_ps(5),ypcorn_ps(6),
     &          xmin,xmax,xtitd(1:len_trim(xtitd)),1,180.,180.,
     &          1.2*xtitoff_ps,
     &          -1.25*xlaboff_ps) !x-axis
            else
              xyplen=sqrt(
     &          (xpcorn_ps(6)-xpcorn_ps(5))**2+
     &          (ypcorn_ps(6)-ypcorn_ps(5))**2
     &          )
              ang=atan2(
     &          (ypcorn_ps(6)-ypcorn_ps(5)),
     &          (xpcorn_ps(6)-xpcorn_ps(5))
     &          )*radtodeg_ps
              iticside=1
              ilabside=1
              titang=180.
              rellabang=-ang
              call mshplt_view_log_axis(
     &          xpcorn_ps(5),ypcorn_ps(5),xyplen,
     &          xminin,xmaxin,ang,
     &          chhe,ticheight,xtit,
     &          iticside,ilabside,
     &          titang,ztitoff_ps*0.7,
     &          rellabang,-2.*xlaboff_ps) !x-axis
            endif

            if (log10y_ps.eq.0) then
              call mshplt_view_axis(
     &          xpcorn_ps(5),xpcorn_ps(8),ypcorn_ps(5),ypcorn_ps(8),
     &          ymin,ymax,ytitd(1:len_trim(ytitd)),-1,0.,0.,
     &          1.1*ytitoff_ps,
     &          1.25*ylaboff_ps) !y-axis
            else
              xyplen=sqrt(
     &          (xpcorn_ps(8)-xpcorn_ps(5))**2+
     &          (ypcorn_ps(8)-ypcorn_ps(5))**2
     &          )
              ang=atan2(
     &          (ypcorn_ps(8)-ypcorn_ps(5)),
     &          (xpcorn_ps(8)-xpcorn_ps(5))
     &          )*radtodeg_ps
              iticside=-1
              ilabside=1
              titang=0.
              rellabang=-ang
              call mshplt_view_log_axis(
     &          xpcorn_ps(5),ypcorn_ps(5),xyplen,
     &          yminin,ymaxin,ang,
     &          chhe,ticheight,ytit,
     &          iticside,ilabside,
     &          titang,ztitoff_ps,
     &          rellabang,zlaboff_ps*0.5) !y-axis
            endif

            if (log10z_ps.eq.0) then
              call mshplt_view_axis(
     &          xpcorn_ps(2),xpcorn_ps(6),ypcorn_ps(2),ypcorn_ps(6),
     &          zmin,zmax,ztitd(1:len_trim(ztitd)),-1,+90.,0.,
     &          ztitoff_ps,zlaboff_ps) !z-axis
            else
              xyplen=sqrt(
     &          (xpcorn_ps(6)-xpcorn_ps(2))**2+
     &          (ypcorn_ps(6)-ypcorn_ps(2))**2
     &          )
              ang=atan2(
     &          (ypcorn_ps(6)-ypcorn_ps(2)),
     &          (xpcorn_ps(6)-xpcorn_ps(2))
     &          )*radtodeg_ps
              iticside=-1
              ilabside=-1
              titang=0.
              rellabang=90.
              call mshplt_view_log_axis(
     &          xpcorn_ps(2),ypcorn_ps(2),xyplen,
     &          zminin,zmaxin,ang,
     &          chhe,ticheight,ztit,
     &          iticside,ilabside,
     &          titang,ztitoff_ps,
     &          rellabang,-1.75*zlaboff_ps) !z-axis
            endif

          else if (phi_ps.le.360.0) then

            if (log10x_ps.eq.0) then
              call mshplt_view_axis(
     &          xpcorn_ps(8),xpcorn_ps(7),ypcorn_ps(8),ypcorn_ps(7),
     &          xmin,xmax,xtitd(1:len_trim(xtitd)),-1,0.,0.,
     &          1.1*xtitoff_ps,
     &          1.2*xlaboff_ps) !x-axis
            else
              xyplen=sqrt(
     &          (xpcorn_ps(7)-xpcorn_ps(8))**2+
     &          (ypcorn_ps(7)-ypcorn_ps(8))**2
     &          )
              ang=atan2(
     &          (ypcorn_ps(7)-ypcorn_ps(8)),
     &          (xpcorn_ps(7)-xpcorn_ps(8))
     &          )*radtodeg_ps
              iticside=-1
              ilabside=1
              titang=0.
              rellabang=-ang
              call mshplt_view_log_axis(
     &          xpcorn_ps(8),ypcorn_ps(8),xyplen,
     &          xminin,xmaxin,ang,
     &          chhe,ticheight,xtit,
     &          iticside,ilabside,
     &          titang,ztitoff_ps*0.7,
     &          rellabang,xlaboff_ps*0.5) !x-axis
            endif

            if (log10y_ps.eq.0) then
              call mshplt_view_axis(
     &          xpcorn_ps(5),xpcorn_ps(8),ypcorn_ps(5),ypcorn_ps(8),
     &          ymin,ymax,ytitd(1:len_trim(ytitd)),-1,0.,0.,
     &          1.2*ytitoff_ps,
     &          1.1*ylaboff_ps) !y-axis
            else
              xyplen=sqrt(
     &          (xpcorn_ps(8)-xpcorn_ps(5))**2+
     &          (ypcorn_ps(8)-ypcorn_ps(5))**2
     &          )
              ang=atan2(
     &          (ypcorn_ps(8)-ypcorn_ps(5)),
     &          (xpcorn_ps(8)-xpcorn_ps(5))
     &          )*radtodeg_ps
              iticside=-1
              ilabside=1
              titang=0.
              rellabang=-ang
              call mshplt_view_log_axis(
     &          xpcorn_ps(5),ypcorn_ps(5),xyplen,
     &          yminin,ymaxin,ang,
     &          chhe,ticheight,ytit,
     &          iticside,ilabside,
     &          titang,ztitoff_ps,
     &          rellabang,zlaboff_ps) !y-axis
            endif

            if (log10z_ps.eq.0) then
              call mshplt_view_axis(
     &          xpcorn_ps(1),xpcorn_ps(5),ypcorn_ps(1),ypcorn_ps(5),
     &          zmin,zmax,ztitd(1:len_trim(ztitd)),-1,90.,0.,
     &          ztitoff_ps,zlaboff_ps) !z-axis
            else
              xyplen=sqrt(
     &          (xpcorn_ps(5)-xpcorn_ps(1))**2+
     &          (ypcorn_ps(5)-ypcorn_ps(1))**2
     &          )
              ang=atan2(
     &          (ypcorn_ps(5)-ypcorn_ps(1)),
     &          (xpcorn_ps(5)-xpcorn_ps(1))
     &          )*radtodeg_ps
              iticside=-1
              ilabside=-1
              titang=0.
              rellabang=90.
              call mshplt_view_log_axis(
     &          xpcorn_ps(1),ypcorn_ps(1),xyplen,
     &          zminin,zmaxin,ang,
     &          chhe,ticheight,ztit,
     &          iticside,ilabside,
     &          titang,ztitoff_ps,
     &          rellabang,-1.75*zlaboff_ps) !z-axis
            endif

          endif !if (phi_ps.le.90.0) then
        endif !if (theta_ps.le.90.0) then

      endif !iaxis.eq.1

      if (ibox.ne.0) then
        call mshplt_line_raw(xpcorn_ps(1),ypcorn_ps(1),xpcorn_ps(2),ypcorn_ps(2))
        call mshplt_line_raw(xpcorn_ps(1),ypcorn_ps(1),xpcorn_ps(4),ypcorn_ps(4))
        call mshplt_line_raw(xpcorn_ps(4),ypcorn_ps(4),xpcorn_ps(8),ypcorn_ps(8))
        call mshplt_line_raw(xpcorn_ps(4),ypcorn_ps(4),xpcorn_ps(3),ypcorn_ps(3))
        call mshplt_line_raw(xpcorn_ps(2),ypcorn_ps(2),xpcorn_ps(3),ypcorn_ps(3))
        call mshplt_line_raw(xpcorn_ps(1),ypcorn_ps(1),xpcorn_ps(5),ypcorn_ps(5))
        call mshplt_line_raw(xpcorn_ps(2),ypcorn_ps(2),xpcorn_ps(6),ypcorn_ps(6))
        call mshplt_line_raw(xpcorn_ps(3),ypcorn_ps(3),xpcorn_ps(7),ypcorn_ps(7))
        call mshplt_pline_raw_closed(4,xpcorn_ps(5),ypcorn_ps(5))
      endif

      call mshplt_get_line_color(ilinecoloro,kredo,kblueo,kgreeno)
      call mshplt_set_line_style(ilinestyleo)

      if (kzone_ps.lt.0) then
        kzone_ps=ifirst_ps
      endif

      write(lun_ps,'(a)')'% end of mshplt_frame3d'

      return
      end
*CMZ :          26/05/2017  10.39.29  by  Michael Scheer
*CMZ :  0.01/02 05/09/2014  15.41.43  by  Michael Scheer
*CMZ :  0.00/06 19/08/2014  15.09.19  by  Michael Scheer
*CMZ :  0.00/04 07/08/2014  15.52.16  by  Michael Scheer
*CMZ :  0.00/03 04/08/2014  15.42.32  by  Michael Scheer
*CMZ :  0.00/02 07/07/2014  12.22.09  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_set_line_width(width)

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      real width
      rlinewidth_ps=width !cm
      write(cline_ps,*)rlinewidth_ps,' setlinewidth'
      call mshplt_fill_buff(cline_ps)

      return
      end
*CMZ :          26/05/2017  10.03.21  by  Michael Scheer
*CMZ :  0.01/02 05/09/2014  15.41.43  by  Michael Scheer
*CMZ :  0.00/06 19/08/2014  15.09.19  by  Michael Scheer
*CMZ :  0.00/04 07/08/2014  15.52.16  by  Michael Scheer
*CMZ :  0.00/03 04/08/2014  15.42.32  by  Michael Scheer
*CMZ :  0.00/02 07/07/2014  12.22.09  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_get_line_width(width)

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      real width
      width=rlinewidth_ps

      return
      end
*CMZ :          14/12/2017  16.53.50  by  Michael Scheer
*-- Author :    Michael Scheer   14/12/2017

! This software is under the GNU General Public License:

!******************************************************************************
!
!      Copyright 2013 Helmholtz-Zentrum Berlin (HZB)
!      Hahn-Meitner-Platz 1
!      D-14109 Berlin
!      Germany
!
!      Author Michael Scheer, Michael.Scheer@Helmholtz-Berlin.de
!
! -----------------------------------------------------------------------
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy (wave_gpl.txt) of the GNU General Public
!    License along with this program.
!    If not, see <http://www.gnu.org/licenses/>.
!
!    Dieses Programm ist Freie Software: Sie koennen es unter den Bedingungen
!    der GNU General Public License, wie von der Free Software Foundation,
!    Version 3 der Lizenz oder (nach Ihrer Option) jeder spaeteren
!    veroeffentlichten Version, weiterverbreiten und/oder modifizieren.
!
!    Dieses Programm wird in der Hoffnung, dass es nuetzlich sein wird, aber
!    OHNE JEDE GEWAEHRLEISTUNG, bereitgestellt; sogar ohne die implizite
!    Gewaehrleistung der MARKTFAEHIGKEIT oder EIGNUNG FueR EINEN BESTIMMTEN ZWECK.
!    Siehe die GNU General Public License fuer weitere Details.
!
!    Sie sollten eine Kopie (wave_gpl.txt) der GNU General Public License
!    zusammen mit diesem Programm erhalten haben. Wenn nicht,
!    siehe <http://www.gnu.org/licenses/>.
!
!******************************************************************************
*CMZ :          31/10/2018  12.03.15  by  Michael Scheer
*-- Author :    Michael Scheer   31/07/2018
      subroutine mshplt_frame_auto(n,x,y,xtit,ytit,chopt)

      implicit none

c call mshplt_frame after auto-scaling

      character(*) xtit,ytit,chopt

      real x(n),y(n)
      real :: xmin, xmax, ymin, ymax, dx, dy

      integer i,n

      xmin=1.0e30
      xmax=-1.0e30
      ymin=1.0e30
      ymax=-1.0e30

      do i=1,n
        if (x(i).lt.xmin) xmin=x(i)
        if (x(i).gt.xmax) xmax=x(i)
        if (y(i).lt.ymin) ymin=y(i)
        if (y(i).gt.ymax) ymax=y(i)
      enddo

      dx=(xmax-xmin)*0.05
      dy=(ymax-ymin)*0.05

      if (dx.eq.0.0) then
        dx=1.
      else if (dx/(abs(xmax)+abs(xmin)).lt.0.001) then
        dx=(abs(xmax)+abs(xmin))*0.0001
      endif

      if (dy.eq.0.0) then
        dy=1.
      else if (dy/(abs(ymax)+abs(ymin)).lt.0.001) then
        dy=(abs(ymax)+abs(ymin))*0.0001
      endif


      call mshplt_frame(xmin-dx,xmax+dx,ymin-dy,ymax+dy,xtit,ytit,chopt)

      return
      end
*CMZ :          06/08/2018  15.37.07  by  Michael Scheer
*-- Author :    Michael Scheer   01/08/2018
      subroutine mshplt_frame3d_xzy_auto(n,x,z,y,
     &  xtit,ztit,ytit,chopt,istatus)

      implicit none

      integer istatus,n,i
      real x(n),y(n),z(n),xmin,ymin,zmin,xmax,ymax,zmax,dx,dy,dz
      character(*) xtit,ztit,ytit,chopt

      xmin=1.0e30
      xmax=-1.0e30
      ymin=1.0e30
      ymax=-1.0e30

      do i=1,n
        if (x(i).lt.xmin) xmin=x(i)
        if (x(i).gt.xmax) xmax=x(i)
        if (y(i).lt.ymin) ymin=y(i)
        if (y(i).gt.ymax) ymax=y(i)
        if (z(i).lt.zmin) zmin=z(i)
        if (z(i).gt.zmax) zmax=z(i)
      enddo

      dx=(xmax-xmin)*0.05
      dy=(ymax-ymin)*0.05
      dz=(zmax-zmin)*0.05

      if(dx.eq.0.0) dx=1.
      if(dy.eq.0.0) dy=1.
      if(dz.eq.0.0) dz=1.

      call mshplt_frame3d_xzy(xmin,xmax,zmin,zmax,ymin,ymax,
     &  xtit,ztit,ytit,chopt)

      istatus=0

      end
*CMZ :          06/08/2018  15.36.01  by  Michael Scheer
*-- Author :    Michael Scheer   01/08/2018
      subroutine mshplt_frame3d_auto(n,x,y,z,
     &  xtit,ytit,ztit,chopt,istatus)

      implicit none

      integer istatus,n,i
      real x(n),y(n),z(n),xmin,ymin,zmin,xmax,ymax,zmax,dx,dy,dz
      character(*) xtit,ztit,ytit,chopt

      xmin=1.0e30
      xmax=-1.0e30
      ymin=1.0e30
      ymax=-1.0e30

      do i=1,n
        if (x(i).lt.xmin) xmin=x(i)
        if (x(i).gt.xmax) xmax=x(i)
        if (y(i).lt.ymin) ymin=y(i)
        if (y(i).gt.ymax) ymax=y(i)
        if (z(i).lt.zmin) zmin=z(i)
        if (z(i).gt.zmax) zmax=z(i)
      enddo

      dx=(xmax-xmin)*0.05
      dy=(ymax-ymin)*0.05
      dz=(zmax-zmin)*0.05

      if (dx.eq.0.0) dx=1.
      if (dy.eq.0.0) dy=1.
      if (dz.eq.0.0) dz=1.

      call mshplt_frame3d(xmin,xmax,ymin,ymax,zmin,zmax,
     &  xtit,ytit,ztit,chopt)

      istatus=0

      end
*CMZ :          01/08/2018  12.27.21  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_set_zlaboff(offset)

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      real offset
      zlaboff_ps=offset

      return
      end
*CMZ :          01/08/2018  12.28.34  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_get_zlaboff(offset)

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      real offset
      offset=zlaboff_ps

      return
      end
*CMZ :          01/08/2018  12.30.47  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_set_ztitoff(offset)

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      real offset
      ztitoff_ps=offset

      return
      end
*CMZ :          01/08/2018  12.30.47  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_get_ztitoff(offset)

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      real offset
      offset=ztitoff_ps

      return
      end
*CMZ :          01/08/2018  15.24.42  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_set_ylaboff(offset)

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      real offset
      ylaboff_ps=offset

      return
      end
*CMZ :          01/08/2018  15.24.42  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_get_ylaboff(offset)

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      real offset
      offset=ylaboff_ps

      return
      end
*CMZ :          01/08/2018  15.24.42  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_set_ytitoff(offset)

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      real offset
      ytitoff_ps=offset

      return
      end
*CMZ :          01/08/2018  15.24.42  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_get_ytitoff(offset)

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      real offset
      offset=ytitoff_ps

      return
      end
*CMZ :          06/08/2018  08.40.38  by  Michael Scheer
*-- Author :    Michael Scheer   06/08/2018
      subroutine mshplt_graph(n,x,y,xtit,ytit,chopt)

      implicit none

      integer n
      real x(n),y(n)

      character(*) xtit,ytit,chopt

      call mshplt_frame_auto(n,x,y,xtit,ytit,chopt)
      call mshplt_pline(n,x,y)

      return
      end
*CMZ :          23/01/2025  15.57.45  by  Michael Scheer
*CMZ :  0.01/03 23/09/2014  09.22.45  by  Michael Scheer
*CMZ :  0.01/02 15/09/2014  15.09.38  by  Michael Scheer
*CMZ :  0.00/04 11/08/2014  18.23.32  by  Michael Scheer
*CMZ :  0.00/02 07/07/2014  12.22.09  by  Michael Scheer
*-- Author :    Michael Scheer   07/07/2014
      subroutine mshplt_arc(x,y,rin,rout,phi1,phi2)

      implicit none

*KEEP,mshpltincl.
      include 'mshplt.cmn'
*KEND.

      real, parameter :: pi=3.14159,grarad=pi/180.
      real :: x,y,rin,rout,x1,y1,phi1,phi2,xpl(2),ypl(2)

      character(32) cphi1,cphi2,carc

      write(cphi1,*) mod(nint(phi1),360)
      cphi1=trim(adjustl(cphi1))

      write(cphi2,*) mod(nint(phi2),360)
      cphi2=trim(adjustl(cphi2))

      if (log10x_ps.eq.0) then
        x1=xleft_ps+scalex_ps*(x-wxmin_ps)
      else
        x1=xleft_ps+scalex_ps*(alog10(x)-wxmin_ps)
      endif
      if (log10y_ps.eq.0) then
        y1=ybottom_ps+scaley_ps*(y-wymin_ps)
      else
        y1=ybottom_ps+scaley_ps*(alog10(y)-wymin_ps)
      endif

      carc=' ' // trim(cphi1) // ' ' // trim(cphi2) // ' arc stroke'
      write(cline_ps,*)'newpath ',x1,y1,rin*2.*chhe_ps,trim(carc)
      write(cline_ps,*)'newpath ',x1,y1,rout*2.*chhe_ps,trim(carc)

      xpl(1)=rin*cos(phi1*grarad)
      ypl(1)=rin*sin(phi1*grarad)
      xpl(2)=rout*cos(phi1*grarad)
      ypl(2)=rout*sin(phi1*grarad)

      call mshplt_pline(2,xpl,ypl)

      xpl(1)=rin*cos(phi2*grarad)
      ypl(1)=rin*sin(phi2*grarad)
      xpl(2)=rout*cos(phi2*grarad)
      ypl(2)=rout*sin(phi2*grarad)

      call mshplt_pline(2,xpl,ypl)

      call mshplt_fill_buff(cline_ps)

      return
      end
*CMZ :          31/07/2018  11.59.03  by  Michael Scheer
*CMZ :  1.01/00 24/09/2014  14.19.07  by  Michael Scheer
*CMZ :  0.01/02 05/09/2014  15.42.16  by  Michael Scheer
*CMZ :  0.01/01 27/08/2014  15.12.25  by  Michael Scheer
*CMZ :  0.01/00 25/08/2014  14.37.06  by  Michael Scheer
*CMZ :  0.00/06 22/08/2014  11.44.36  by  Michael Scheer
*CMZ :  0.00/04 18/08/2014  08.59.57  by  Michael Scheer
*CMZ :  1.17/00 06/05/2014  19.24.52  by  Michael Scheer
*CMZ :  1.16/04 16/04/2014  16.30.12  by  Michael Scheer
*-- Author :    Michael Scheer   16/04/2014
      subroutine mplint(jdev)

      implicit none

      integer idev,jdev

*KEEP,MPLOTINCL.
      include 'mplot.cmn'
*KEND.
      integer lunbase
      logical lexist
      real xsiz,ysiz
      integer ibbxl,ibbyb,ibbxr,ibbyt,lenfile,lenview,lenkview,i

      idev=jdev



      lunbase=100000
      viewer_mshplt=''
      fileeps_mshplt='mshplt.eps'
      xsiz=-20.
      ysiz=-20.
      rescale_mshplt=1.
      inquire(file='.mshplt.cnf',exist=lexist)
      if (lexist) then
c        call util_get_free_lun(lunbase)
        open(newunit=lunbase,file='.mshplt.cnf',status='old')
        read(lunbase,*)idev
        read(lunbase,*)xsiz,ysiz
        read(lunbase,*)rescale_mshplt
        read(lunbase,*)ibbxl, ibbyb,ibbxr,ibbyt
        read(lunbase,'(a)')fileeps_mshplt
        if(idev.ne.0) then
          read(lunbase,'(a)')viewer_mshplt
          read(lunbase,'(a)')viewer_kill_mshplt
        endif
        close(lunbase)
      endif

      lenfile=0
      do i=1,len_trim(fileeps_mshplt)
        if (fileeps_mshplt(i:i).eq.'!') goto 99
        lenfile=i
      enddo

99    continue

      lenview=0
      do i=1,len_trim(viewer_mshplt)
        if (viewer_mshplt(i:i).eq.'!') goto 99
        lenview=i
      enddo

      lenkview=0
      do i=1,len_trim(viewer_kill_mshplt)
        if (viewer_kill_mshplt(i:i).eq.'!') goto 99
        lenkview=i
      enddo

      if (idev.ne.0) idev=1

      ! negative size values indicate HIGZ mode for mshplt
      call mshplt_init(idev,
     &  xsiz,ysiz,
     &  ibbxl,ibbyb,ibbxr,ibbyt,
     &  fileeps_mshplt(1:lenfile)
     &  ,viewer_mshplt(1:lenview)
     &  ,viewer_kill_mshplt(1:lenkview),rescale_mshplt
     &  )
      call mplzon(1,1,1,'')
      return
      end
*CMZ :  0.00/06 19/08/2014  18.59.39  by  Michael Scheer
*CMZ :  1.17/00 17/04/2014  18.23.45  by  Michael Scheer
*CMZ :  1.16/04 16/04/2014  14.56.04  by  Michael Scheer
*-- Author :    Michael Scheer   16/04/2014
      subroutine mpl3(n,x,y,z)

      implicit none

      integer n
     &  ,icolor,ired,igreen,iblue
      real x(n),y(n),z(n)



      call mshplt_get_line_color(
     &  icolor,ired,igreen,iblue
     &  )
      call mshplt_set_line_color(
     &  icolor,ired,igreen,iblue
     &  )

      call mshplt_pline3d(n,x,y,z)

      call mshplt_get_color(
     &  icolor,ired,igreen,iblue
     &  )
      call mshplt_set_color(
     &  icolor,ired,igreen,iblue
     &  )

      return
      end
*CMZ :  0.00/06 20/08/2014  09.34.38  by  Michael Scheer
*CMZ :  1.17/00 01/05/2014  09.49.35  by  Michael Scheer
*CMZ :  1.16/04 16/04/2014  14.57.09  by  Michael Scheer
*-- Author :    Michael Scheer   16/04/2014
      subroutine mpl(n,x,y)

      implicit none

      integer n
      real x(n),y(n)


      call mshplt_pline(n,x,y)
      return
      end
*CMZ :  0.01/00 23/08/2014  11.19.02  by  Michael Scheer
*CMZ :  0.00/06 20/08/2014  09.19.28  by  Michael Scheer
*CMZ :  1.17/00 06/05/2014  16.11.26  by  Michael Scheer
*CMZ :  1.16/04 16/04/2014  15.09.15  by  Michael Scheer
*-- Author :    Michael Scheer   16/04/2014
      subroutine mplax(xtit,ytit)

      implicit none

      character*(*) xtit,ytit
      character xtitd,ytitd

      xtitd=xtit(1:1)
      ytitd=ytit(1:1)


      call mshplt_draw_axis_titles(xtit,ytit)

      return
      end
*CMZ :  0.00/06 20/08/2014  09.18.19  by  Michael Scheer
*CMZ :  1.17/00 17/04/2014  18.25.11  by  Michael Scheer
*CMZ :  1.16/04 16/04/2014  15.11.09  by  Michael Scheer
*-- Author :    Michael Scheer   16/04/2014
      subroutine muwk(iwkid,irflg)

      implicit none

      integer iwkid,irflg
      integer iwkidd,irflgd

      iwkidd=iwkid
      irflgd=irflg


      call mshplt_flush_buff
      return
      end
*CMZ :  1.03/02 08/12/2015  13.47.56  by  Michael Scheer
*CMZ :  1.02/00 03/10/2014  10.05.12  by  Michael Scheer
*CMZ :  0.00/06 19/08/2014  14.25.00  by  Michael Scheer
*CMZ :  1.17/00 26/05/2014  11.00.23  by  Michael Scheer
*CMZ :  1.16/04 16/04/2014  15.13.52  by  Michael Scheer
*-- Author :    Michael Scheer   16/04/2014
      subroutine mgset(pname,rval)

      implicit none

      character(4) pname
      real rval



      if (pname.eq.'MSCF') then
        call mshplt_igset(pname,rval/10.)
      else
        call mshplt_igset(pname,rval)
      endif
      return
      end
*CMZ :  0.00/06 19/08/2014  19.48.10  by  Michael Scheer
*CMZ :  1.17/00 22/04/2014  15.41.26  by  Michael Scheer
*CMZ :  1.16/04 16/04/2014  15.16.07  by  Michael Scheer
*-- Author :    Michael Scheer   16/04/2014
      subroutine mtx(x,y,chars)

      implicit none

      real x,y
      character*(*) chars


      call mshplt_text(x,y,chars)

      return
      end
*CMZ :          23/02/2022  10.34.24  by  Michael Scheer
*CMZ :  0.00/06 20/08/2014  12.45.42  by  Michael Scheer
*CMZ :  0.00/04 14/08/2014  11.23.35  by  Michael Scheer
*CMZ :  1.17/00 22/04/2014  15.58.15  by  Michael Scheer
*CMZ :  1.16/04 16/04/2014  15.22.38  by  Michael Scheer
*-- Author :    Michael Scheer   16/04/2014
      subroutine mplfra(xmin,xmax,ymin,ymax,chopt)

      implicit none

      real xmin,xmax,ymin,ymax
      integer i
      character*(*) chopt
      character(6) copt



      copt='LrBtc'

      do i=1,len_trim(chopt)
        if (chopt(i:i).eq.'S'.or.chopt(i:i).eq.'s') then
          copt(6:6)='S'
        endif
      enddo

      do i=1,len_trim(chopt)
        if (chopt(i:i).eq.'A'.or.chopt(i:i).eq.'a') then
c          copt(1:1)='l'
c          copt(3:3)='b'
        endif
      enddo

      do i=1,len_trim(chopt)
        if (chopt(i:i).eq.'B'.or.chopt(i:i).eq.'b') then
          copt(1:5)='xxxxx'
        endif
      enddo

      call mshplt_frame(xmin,xmax,ymin,ymax,'','',copt)

      return
      end
*CMZ :  1.03/01 07/11/2014  19.48.39  by  Michael Scheer
*CMZ :  0.00/06 20/08/2014  09.15.17  by  Michael Scheer
*CMZ :  1.16/04 16/04/2014  15.25.49  by  Michael Scheer
*-- Author :    Michael Scheer   16/04/2014
      subroutine mzitoc(ival,cval)

      implicit none

      integer ival,k,i
      character*(*) cval
      character(12) c12

      write(c12,*)ival
      k=0
      cval=''
      do i=1,min(12,len(cval))
        if (c12(i:i).ne.' ') then
          k=k+1
          cval(k:k)=c12(i:i)
        endif
      enddo

      return
      end
*CMZ :  1.02/00 02/10/2014  15.04.37  by  Michael Scheer
*CMZ :  0.00/06 20/08/2014  11.28.27  by  Michael Scheer
*CMZ :  1.17/00 06/05/2014  09.27.26  by  Michael Scheer
*CMZ :  1.16/04 16/04/2014  15.27.58  by  Michael Scheer
*-- Author :    Michael Scheer   16/04/2014
      subroutine mplfr3(xmin,xmax,ymin,ymax,zmin,zmax,theta,phi,chopt)

      implicit none

      real xmin,xmax,ymin,ymax,zmin,zmax,theta,phi
      character*(*) chopt
      character(8) copt
      integer ko



      copt=chopt(1:len_trim(chopt))
      call mshplt_set_theta_phi(theta,phi)
      call mshplt_get_box(ko)
      call mshplt_set_box(0)
      call mshplt_frame3d(xmin,xmax,ymin,ymax,zmin,zmax,
     &  '','','',copt)
      call mshplt_set_box(ko)

      return
      end
*CMZ :  0.00/06 20/08/2014  09.17.11  by  Michael Scheer
*CMZ :  1.17/00 26/05/2014  09.48.06  by  Michael Scheer
*CMZ :  1.16/04 16/04/2014  15.31.14  by  Michael Scheer
*-- Author :    Michael Scheer   16/04/2014
      subroutine mpm(n,x,y)

      implicit none

      integer n
      real x(*),y(*)


      call mshplt_marker(n,x,y)

      return
      end
*CMZ :          07/08/2018  14.43.55  by  Michael Scheer
*CMZ :  0.00/04 13/08/2014  16.44.46  by  Michael Scheer
*CMZ :  1.17/00 26/05/2014  09.35.25  by  Michael Scheer
*CMZ :  1.16/04 16/04/2014  15.36.34  by  Michael Scheer
*-- Author :    Michael Scheer   16/04/2014
      subroutine mplzon(ix,iy,iwi,ckic)

      implicit none
      integer ix,iy,iwi
      character(*) ckic


      call mshplt_zone(ix,iy,iwi,ckic)
      return
      end
*CMZ :  0.00/06 19/08/2014  13.33.59  by  Michael Scheer
*CMZ :  1.17/00 17/04/2014  17.15.31  by  Michael Scheer
*CMZ :  1.16/04 16/04/2014  15.51.04  by  Michael Scheer
*-- Author :    Michael Scheer   16/04/2014
      subroutine mplset(case,arg)

      implicit none

      real arg
      character(4) case



      call mshplt_hplset(case,arg)

      return
      end
*CMZ :  0.01/00 22/08/2014  16.39.59  by  Michael Scheer
*CMZ :  0.00/06 20/08/2014  09.12.05  by  Michael Scheer
*CMZ :  1.17/00 06/05/2014  11.03.18  by  Michael Scheer
*CMZ :  1.16/04 16/04/2014  15.41.44  by  Michael Scheer
*-- Author :    Michael Scheer   16/04/2014
      subroutine mtitle(tit)

      implicit none

      character*(*) tit
      character tdum

      tdum=tit(1:1)



      call mshplt_title(tit)
      return
      end
*CMZ :  0.01/02 05/09/2014  15.42.16  by  Michael Scheer
*CMZ :  0.01/01 26/08/2014  08.34.13  by  Michael Scheer
*CMZ :  0.00/04 13/08/2014  15.55.54  by  Michael Scheer
*CMZ :  1.17/00 06/05/2014  19.07.28  by  Michael Scheer
*CMZ :  1.16/04 16/04/2014  15.43.29  by  Michael Scheer
*-- Author :    Michael Scheer   16/04/2014
      subroutine mplend

      implicit none

*KEEP,mplotincl.
      include 'mplot.cmn'
*KEND.


      call mshplt_end
      return
      end
*CMZ :  0.00/06 20/08/2014  09.16.07  by  Michael Scheer
*CMZ :  1.17/00 17/04/2014  17.36.30  by  Michael Scheer
*CMZ :  1.16/04 16/04/2014  15.45.07  by  Michael Scheer
*-- Author :    Michael Scheer   16/04/2014
      subroutine mgmeta(iunit,metafl)

      implicit none

      integer iunit,metafl
      integer iunitd,metafld

      iunitd=iunit
      metafld=metafl


      return
      end
*CMZ :  1.02/00 02/10/2014  16.21.50  by  Michael Scheer
*CMZ :  0.01/02 17/09/2014  13.06.35  by  Michael Scheer
*CMZ :  0.00/06 19/08/2014  13.52.14  by  Michael Scheer
*CMZ :  1.17/00 17/04/2014  17.28.58  by  Michael Scheer
*CMZ :  1.16/04 16/04/2014  15.50.46  by  Michael Scheer
*-- Author :    Michael Scheer   16/04/2014
      subroutine mplopt(choptn,n)

      implicit none

      integer n
c      character(4) chopt(*)
      character(*) choptn



      call mshplt_hplopt(choptn,n)

      return
      end
*CMZ :  0.00/06 20/08/2014  09.10.55  by  Michael Scheer
*CMZ :  1.16/04 16/04/2014  15.49.21  by  Michael Scheer
*-- Author :    Michael Scheer   16/04/2014
      subroutine mlimit(ndata)

      implicit none
      integer ndata,ndatadum

      ndatadum=ndata


      return
      end
*CMZ :  1.00/01 23/09/2014  14.28.02  by  Michael Scheer
*CMZ :  0.01/02 04/09/2014  14.47.39  by  Michael Scheer
*-- Author :    Michael Scheer   04/09/2014
      subroutine mplbox(xleft,xright,ybottom,ytop,chopt)
      real xleft,xright,ybottom,ytop
      character chopt,choptd

      choptd=chopt
      call mshplt_box(xleft,xright,ybottom,ytop)

      return
      end
*CMZ :  1.03/01 10/11/2014  14.37.11  by  Michael Scheer
*CMZ :  0.00/06 20/08/2014  09.15.17  by  Michael Scheer
*CMZ :  1.16/04 16/04/2014  15.25.49  by  Michael Scheer
*-- Author :    Michael Scheer   16/04/2014
      subroutine miztoc(ival,cval)

      implicit none

      integer ival,k,i
      character*(*) cval
      character(12) c12

      write(c12,*)ival
      k=0
      cval=''
      do i=1,min(12,len(cval))
        if (c12(i:i).ne.' ') then
          k=k+1
          cval(k:k)=c12(i:i)
        endif
      enddo

      return
      end
*CMZ :  1.03/01 10/10/2014  16.01.29  by  Michael Scheer
*CMZ : 00.00/07 21/07/2009  14.58.29  by  Michael Scheer
*CMZ : 00.00/06 12/07/2007  15.45.32  by  Michael Scheer
*-- Author :    Michael Scheer   12/07/2007
      subroutine util_get_free_lun(lun)

      integer lun,luni

      logical isopen

      luni=abs(lun)

1     continue
      inquire(unit=luni,opened=isopen)
      if (isopen) then
        luni=luni+1
        goto 1
      endif

      lun=luni

      return
      end
*CMZ :  1.00/01 24/09/2014  12.06.33  by  Michael Scheer
*CMZ :  0.01/03 23/09/2014  12.14.42  by  Michael Scheer
*CMZ : 00.00/06 07/03/2007  17.00.51  by  Michael Scheer
*CMZ : 00.00/05 07/03/2007  12.58.44  by  Michael Scheer
*-- Author :    Michael Scheer   07/03/2007
      subroutine util_string_split(cline,ndim,nwords,ipos,istat)

c Input:
c      cline
c      ndim: max. number of words

c Output:
c      nwords: number of words
c      ipos(1:2,iword): start and end position of word in cline

      implicit none

      integer ilen,istat,jx,nwords,ndim,ipos(2,*),iblank,in
      character(*) cline

      byte ic
      character c1
      equivalence(c1,ic)

      istat=-1

      ilen=len_trim(cline)

      in=0
      nwords=0

      do jx=1,ilen

        c1=cline(jx:jx)

        if (ic.le.32.or.ic.ge.127) then
c found invisible character
          iblank=1
        else
          iblank=0
        endif

        if (iblank.eq.1) then

          if (in.eq.1) then
c found end of word
            if (jx.lt.ilen) then
              ipos(2,nwords)=jx-1
            else
              ipos(2,nwords)=jx
            endif
            in=0
          endif

        else !iblank

          if (in.eq.0) then
c found beginning of word
            nwords=nwords+1
            if (nwords.gt.ndim) return
            ipos(1,nwords)=jx
            ipos(2,nwords)=ilen
            in=1
          endif

        endif !iblank

      enddo


      if (nwords.gt.0) then
        istat=0
      endif

      return
      end
*CMZ :  1.01/00 24/09/2014  13.42.26  by  Michael Scheer
*CMZ :  0.01/02 04/09/2014  12.59.12  by  Michael Scheer
*CMZ : 00.00/15 12/10/2013  12.19.05  by  Michael Scheer
*CMZ : 00.00/11 07/06/2011  14.38.25  by  Michael Scheer
*CMZ : 00.00/07 07/06/2011  13.46.51  by  Michael Scheer
*CMZ : 00.00/06 06/07/2007  16.55.15  by  Michael Scheer
*CMZ : 00.00/02 25/08/2006  15.16.22  by  Michael Scheer
*CMZ : 00.00/01 23/02/96  14.56.50  by  Michael Scheer
*CMZ : 00.00/00 10/01/95  15.27.54  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE UTIL_SPLINE_INTERPOLATION_f90(N,XA,YA,X,Y,Y2A,MODE)

C---  INTERPOLATES Y(X) VIA SPLINE

C--   INPUT:

C-       N: NUMBER OF X,Y-VALUES
C-       XA:   ARRAY OF X-VALUES
C-       YA:   ARRAY OF Y-VALUES
C-       X: Y(X) IS CALCULATED
C-       MODE: CONTROL FLAG:
C-             MODE.GE.0: USE VALUES OF LAST CALL TO START WITH
C-             MODE.LT.0: NEW INITIALIZATION

C--   OUTPUT:

C-       Y: Y(X) IS CALCULATED
C-     Y2A: SECOND DERIVATE OF YA(XA), I.E. SPLINE COEFFICIENTS

      IMPLICIT NONE

      double precision Y,X,XA1OLD,XANOLD,H,A,B
      double precision XA(*),YA(*),Y2A(*)

      INTEGER NOLD,N,KLO,KHI,KLOLD,K,MODE

      DATA KLOLD/1/,NOLD/-99/
      DATA XA1OLD/-9999.D0/,XANOLD/-9999./

      IF(N.lt.3) THEN
        STOP '*** ERROR IN UTIL_SPLINE_INTER_F90: LESS THAN 3 POINTS GIVEN***'
      ENDIF

      IF(XA(1).LT.XA(N).AND.(X.LT.XA(1).OR.X.GT.XA(N))
     &    .OR.
     &    XA(N).LT.XA(1).AND.(X.LT.XA(N).OR.X.GT.XA(1))) THEN
        WRITE(6,*)'XA(1), XA(N):',XA(1), XA(N)
        WRITE(6,*)'X:',X
        STOP '*** ERROR IN UTIL_SPLINE_INTER_F90: X OUT OF RANGE ***'
      ENDIF

      IF (MODE.LT.0.OR.KLOLD.GE.N) THEN
        CALL UTIL_SPLINE_COEF_F90(XA,YA,N,9999.0d0,9999.0d0,Y2A)
        KLO=1
      ELSE IF(NOLD.EQ.N
     &    .AND. XA(1).EQ.XA1OLD
     &    .AND. XA(N).EQ.XANOLD
     &    .AND. X.GT.XA(KLOLD)
     &    ) THEN
        KLO=KLOLD
      ELSE
        KLO=1
      ENDIF

      IF (X.LT.XA(KLO+1)) THEN
        KHI=KLO+1
        GOTO 2
      ENDIF

      KHI=N
1     IF (KHI-KLO.GT.1) THEN
        K=(KHI+KLO)/2
        IF(XA(K).GT.X)THEN
          KHI=K
        ELSE
          KLO=K
        ENDIF
        GOTO 1
      ENDIF

2     H=XA(KHI)-XA(KLO)

      IF (H.le.0.D0) THEN
        PRINT*,'*** ERROR IN UTIL_SPLINE_INTER_F90: BAD INPUT ***'
        STOP
      ENDIF

      A=(XA(KHI)-X)/H
      B=(X-XA(KLO))/H
      Y=A*YA(KLO)+B*YA(KHI)+
     &  (A*(A+1.D0)*(A-1.D0)*Y2A(KLO)+B*(B+1.D0)*
     &  (B-1.D0)*Y2A(KHI))*(H**2)/6.D0

      KLOLD=KLO
      NOLD=N
      XA1OLD=XA(1)
      XANOLD=XA(N)

      RETURN
      END
*CMZ :  1.01/00 24/09/2014  13.42.26  by  Michael Scheer
*CMZ : 00.00/07 07/05/2008  14.02.20  by  Michael Scheer
*CMZ :  2.47/07 14/04/2003  15.17.05  by  Michael Scheer
*CMZ :  2.20/09 03/04/2001  11.03.37  by  Michael Scheer
*CMZ :  2.15/00 28/04/2000  10.32.37  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.55.48  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.14.48  by  Michael Scheer
*-- Author :
      SUBROUTINE UTIL_SPLINE_COEF_F90(X,Y,N,YP1,YPN,Y2)

C--- CALCULATES SPLINE COEFFICIENTS

C--   INPUT:

C-       N: NUMBER OF X,Y-VALUES
C-       X: ARRAY OF X-VALUES
C-       Y: ARRAY OF Y-VALUES
C-       YP1:  SECOND DERIVATIVE AT FIRST X-VALUE
C-       YPN:  SECOND DERIVATIVE AT LAST X-VALUE

C--   OUPUT:

C-       Y2:   SPLINE-COEFFICIENTS

C--   WORKINGSPACE: AA(N),BB(N),CC(N),C(N)


      IMPLICIT NONE

      INTEGER N,J
      DOUBLE PRECISION  X(*),Y(*),Y2(*)
      DOUBLE PRECISION YP1,YPN

      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: AA,BB,CC,C

      double precision xx(3),yy(3),a(3),yp(3),xopt,yopt
      INTEGER ifail

      ALLOCATE(AA(N))
      ALLOCATE(BB(N))
      ALLOCATE(CC(N))
      ALLOCATE(C(N))

      IF (N.LT.3) then
        if (abs(yp1).eq.9999.0d0) then
          y2(1)=0.0d0
        else
          y2(1)=yp1
        endif
        if (abs(ypn).eq.9999.0d0) then
          y2(n)=0.0d0
        else
          y2(n)=ypn
        endif
        GOTO 999
      endif

      if (abs(yp1).eq.9999.0d0) then
        xx=x(1:3)
        yy=y(1:3)
        call UTIL_PARABEL(xx,yy,A,YP,XOPT,yopt,IFAIL)
        if (ifail.eq.0) then
          y2(1)=2.0d0*a(3)
        else
          y2(1)=0.0d0
        endif
      else
        Y2(1)=YP1
      endif

      if (abs(ypn).eq.9999.0d0) then
        xx=x(n-2:n)
        yy=y(n-2:n)
        call UTIL_PARABEL(xx,yy,A,YP,XOPT,yopt,IFAIL)
        if (ifail.eq.0) then
          y2(n)=2.0d0*a(3)
        else
          y2(N)=0.0d0
        endif
      else
        Y2(N)=YPN
      endif

      C(1)=Y2(1)
      C(N)=y2(n)

      BB(1)=1.D0
      CC(1)=0.D0
      CC(N)=1.D0

      DO J=2,N-1
          AA(J)=(X(J  )-X(J-1))/6.D0
          BB(J)=(X(J+1)-X(J-1))/3.D0
          CC(J)=(X(J+1)-X(J  ))/6.D0
          C(J)=(Y(J+1)-Y(J  ))/(X(J+1)-X(J  ))
     &          -(Y(J  )-Y(J-1))/(X(J  )-X(J-1))
      ENDDO !J

      DO J=2,N-1

          BB(J)=BB(J)-AA(J)*CC(J-1)
           C(J)= C(J)-AA(J)* C(J-1)
C030414          AA(J)=AA(J)-AA(J)*BB(J-1)

          CC(J)=CC(J)/BB(J)
           C(J)= C(J)/BB(J)
          BB(J)=1.D0

      ENDDO !J

      DO J=N-1,2,-1
         Y2(J)=C(J)-CC(J)*Y2(J+1)
      ENDDO

999   CONTINUE

      DEALLOCATE(AA)
      DEALLOCATE(BB)
      DEALLOCATE(CC)
      DEALLOCATE(C)

      RETURN
      END
*CMZ :  0.01/02 19/03/2014  12.30.26  by  Michael Scheer
*CMZ : 00.00/15 03/09/2012  09.27.13  by  Michael Scheer
*CMZ : 00.00/07 22/03/2010  15.28.00  by  Michael Scheer
*CMZ : 00.00/02 26/03/97  10.23.11  by  Michael Scheer
*CMZ : 00.00/00 10/01/95  15.27.40  by  Michael Scheer
*-- Author :
      SUBROUTINE UTIL_PARABEL(Xin,Yin,A,YP,XOPT,yopt,IFAIL)

C--- CALCULATES A(1),A(2),A(3), THE DERIVATIVES YP(X(1)),YP(X(2)),YP(X(3)),
C    AND THE EXTREMUM (XOPT,A(XOPT)) OF PARABOLA Y=A1+A2*X+A3*X**2
C    FROM COORDINATES OF THE THREE POINTS (X(1),Y(1)),(X(2),Y(2)),(X(3),Y(3))
C

      IMPLICIT NONE

      INTEGER IFAIL

      REAL*8 A(3),X(3),Y(3),DXM,DXP,x0,a1,a2,dxm2,dxp2,dxmax,dymax,
     &  xin(3),yin(3)
      REAL*8 DET,YP(3),XOPT,yopt,a22,fm,fp,f0

      IFAIL=0

c calculate f=a0+a1*(x-x0)+a2*(x-x0)**2
c  = a0 + a1*x - a0*x0 + a2*x**2 - 2*a2*x*x0 + a2*x0**2
c  = a0 + (a2*x0 -a0)*x0 + (a1 - 2*a2*x0 )*x+ a2*x**2

c change system: (x0,s0)->(0,0), i.e.
c calculate f=a1*dx+a2*dx**2
c  df/dx=a1+2*a2*dx_max =! 0, dx_max=-a1/2/a2

      x=xin
      y=yin

      call util_sort_func(3,x,y)

      x0=x(2)
      f0=y(2)

      fm=y(1)-f0
      fp=y(3)-f0

      dxm=x(1)-x0
      dxp=x(3)-x0

c fm=a1*dxm+a2*dxm**2
c fp=a1*dxp+a2*dxp**2

c (dxm dxm2) (a1) = (y(1))
c (dxp dxp2) (a2) = (y(3))

      dxm2=dxm*dxm
      dxp2=dxp*dxp

      det=dxm*dxp2-dxp*dxm2

      if (det.ne.0.0d0) then
        a1=(fm*dxp2-fp*dxm2)/det
        a2=(fp*dxm-fm*dxp)/det
      else
        ifail=1
        return
      endif

      if (a2.ne.0.0d0) then
        dxmax=-a1/(2.0d0*a2)
        dymax=(a1+a2*dxmax)*dxmax
        xopt=x0+dxmax
        yopt=f0+dymax
      endif

c calculate f=f0+a1*dx+a2*dx**2
c = a1*x - a1*x0 + a2*x**2 + a2*x0**2 - 2*a2*x*x0
c  f = f0 + (a2*x0 -a1)*x0 + (a1 - 2*a2*x0 )*x+ a2*x**2

      a22=2.0d0*a2

      a(1)=f0 + (a2*x0 -a1)*x0
      a(2)=a1 - a22*x0
      a(3)=a2

c calculate yp=a1+2*a2*dx

      yp(1)=a1+a22*dxm
      yp(2)=a1
      yp(3)=a1+a22*dxp

      RETURN
      END
*CMZ : 00.00/15 05/01/2012  13.52.39  by  Michael Scheer
*CMZ : 00.00/00 11/01/95  11.41.04  by  Michael Scheer
*-- Author :
      SUBROUTINE UTIL_SORT_FUNC(N,RA,YA)

C--- HEAPSORT ROUTINE; SEE NUMERICAL RECIPES 8.2 S 231
C--- ARRAY YA IS FUNCTION OF RA AND SORTED ACCORDINGLY

      IMPLICIT NONE

      INTEGER N,L,IR,I,J

      REAL*8 RA(N),RRA
      REAL*8 YA(N),YYA

      if (n.lt.2) return

      L=N/2+1
      IR=N
10    CONTINUE
        IF(L.GT.1)THEN
          L=L-1
          RRA=RA(L)
          YYA=YA(L)
        ELSE
          RRA=RA(IR)
          YYA=YA(IR)
          RA(IR)=RA(1)
          YA(IR)=YA(1)
          IR=IR-1
          IF(IR.EQ.1)THEN
            RA(1)=RRA
            YA(1)=YYA
            RETURN
          ENDIF
        ENDIF
        I=L
        J=L+L
20      IF(J.LE.IR)THEN
          IF(J.LT.IR)THEN
            IF(RA(J).LT.RA(J+1))J=J+1
          ENDIF
          IF(RRA.LT.RA(J))THEN
            RA(I)=RA(J)
            YA(I)=YA(J)
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        RA(I)=RRA
        YA(I)=YYA
      GO TO 10
      END
