*CMZ :  2.63/02 17/03/2008  16.12.18  by  Michael Scheer
*CMZ :  2.62/01 17/04/2007  13.58.17  by  Michael Scheer
*CMZ :  2.61/04 29/03/2007  15.30.35  by  Michael Scheer
*CMZ :  2.61/02 20/03/2007  16.55.13  by  Michael Scheer
*CMZ :  2.57/05 08/01/2007  18.02.51  by  Michael Scheer
*-- Author :    Michael Scheer   13/12/2006
      program wave_qpu_loop

      implicit none

      double precision b0hal,b0err,eta,b0halmn,b0halmx,b0errmn,b0errmx,
     &  etamn,etamx,db0err,db0hal,deta,energy,fmax,width

      integer neta,nb0err,nb0hal,icut,ib0err,ib0hal,ieta,idum,
     &  icutmn,icutmx,kcut,jcut,iavail(10000),icode,iread,ib0,ifibocut

      open(unit=99,file='wave_qpu_loop.in',status='old')

      read(99,*)neta,etamn,etamx
      read(99,*)icutmn,icutmx
      read(99,*)nb0hal,b0halmn,b0halmx
      read(99,*)nb0err,b0errmn,b0errmx

      close(99)

      call system('echo \*\*\*\*\*\*\*\*\*\*\ `date` >> wave_qpu_loop.dat')
      call system('echo \*\*\*\*\*\*\*\*\*\*\ `date` >> wave_qpu_loop_int.dat')

      open(unit=20,file='wave_qpu_loop.dat',status='unknown',access='append')

      db0hal=0.0d0
      db0err=0.0d0
      deta=0.0d0

      if (nb0hal.gt.1) db0hal=(b0halmx-b0halmn)/(nb0hal-1)
      if (nb0err.gt.1) db0err=(b0errmx-b0errmn)/(nb0err-1)
      if (neta.gt.1) deta=(etamx-etamn)/(neta-1)

      kcut=0
      b0hal=b0halmn-db0hal

      do ib0hal=1,nb0hal

        b0hal=b0hal+db0hal
        b0err=b0errmn-db0err
        ib0=0

        do ib0err=1,nb0err

          b0err=b0err+db0err
          eta=etamn-deta

          do ieta=1,neta

            eta=eta+deta

            icut=0
            kcut=1
            iavail=0
            jcut=icutmx

21          icut=icut+1

            if (icut.eq.1.and.ib0err.eq.1) then

              open(unit=99,file='wave_qpu_uname.dat',status='unknown')

              if (ib0.eq.0) then
                write(99,*)b0hal,' 0.0d0',' 1'
              else
                write(99,*)b0hal,' 0.0d0',' 0'
              endif

              write(99,*)eta,kcut

              close(99)

              call system('sh wave_qpu_loop.sh')

              iread=0
              open(unit=99,file='wave_qpu_uout.dat',status='old',err=99)

              read(99,*)ifibocut

              if (ifibocut.gt.0) then
                ib0=1
              endif

2             read(99,*,end=93)icode,energy,fmax,width
              iread=iread+1

              write(20,"(i8,g13.6,i4,'  ',5g12.5)")
     &          icode,eta,kcut,b0hal,' 0.0d0 ',energy,fmax,width

              goto 2
93            close(99)

              if (iread.gt.0) then
                call system('cat wave_qpu_uout_int.dat >> wave_qpu_loop_int.dat')
              endif

              open(unit=99,file='fibonacci_available_cuts.dat',status='old',err=92)

              read(99,*)
              read(99,*)

              jcut=0

11            read(99,*,end=12)iread

              jcut=jcut+1

              if (jcut.gt.10000) then
                stop '*** Error in wave_qpu_loop.exe: Dimension exceeded'
              endif

              iavail(jcut)=iread

              goto 11
12            close(99)

            endif !icut

            if (iavail(icut).ne.0) then
              kcut=iavail(icut)
            else
              if (icut.le.icutmx) write(6,*)'no cut for eta, icut:',eta,icut
              goto 888 ! next eta
            endif

            open(unit=99,file='wave_qpu_uname.dat',status='unknown')

            write(99,*)b0hal,b0err,' 1'
            write(99,*)eta,kcut

            close(99)

            call system('sh wave_qpu_loop.sh')

92          iread=0
            open(unit=99,file='wave_qpu_uout.dat',status='old',err=99)

            read(99,*)idum

1           read(99,*,end=91)icode,energy,fmax,width
            iread=iread+1
            write(20,"(i8,g13.6,i4,'  ',5g12.5)")
     &        icode,eta,kcut,b0hal,b0err,energy,fmax,width
            goto 1

91          close(99)

99          continue

            if (iread.gt.0) then
              call system('cat wave_qpu_uout_int.dat >> wave_qpu_loop_int.dat')
            endif

            if (icut.le.icutmx.and.icut.le.jcut) goto 21

888         continue
          enddo !ieta
        enddo !ib0err
      enddo !ib0hal

      close(20)

      stop
      end
