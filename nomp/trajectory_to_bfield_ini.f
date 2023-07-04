*CMZ :  4.00/11 21/04/2021  12.06.18  by  Michael Scheer
*-- Author :    Michael Scheer   21/04/2021
      subroutine trajectory_to_bfield_ini
      implicit none

*KEEP,contrl.
      include 'contrl.cmn'
*KEEP,myfiles.
      include 'myfiles.cmn'
*KEEP,phyconparam.
      include 'phyconparam.cmn'
*KEND.

      double precision, dimension (:,:), allocatable :: xyz,bxyz
      double precision x,y,z

      integer lunt,i,nstep,istat
      logical lexist

      character(128) chfilebx,chfileby,chfilebz

      inquire(file=trim(filetr),exist=lexist)
      if (lexist.eqv..false.) then
        print*,'*** Error in trajectory_to_bfield_ini: File'
        print*,trim(filetr)
        print*,"not found ***"
        write(lungfo,*)'*** Error in trajectory_to_bfield_ini: File'
        write(lungfo,*)trim(filetr)
        write(lungfo,*)"not found ***"
        return
      endif

      chfilebx='trajectory_to_bfield_bx.dat'
      chfileby='trajectory_to_bfield_by.dat'
      chfilebz='trajectory_to_bfield_bz.dat'

      filetbx=chfilebx
      filetb=chfileby
      filetbz=chfilebz

      irbtabxyz=1

      open(newunit=lunt,file=trim(filetr),status='old')

      nstep=0

      do while(.true.)
        call util_skip_comment(lunt)
        read(lunt,*,end=9)x,y,z
        nstep=nstep+1
      enddo
9     continue

      allocate(xyz(3,nstep),bxyz(3,nstep))

      rewind(lunt)
      do i=1,nstep
        call util_skip_comment(lunt)
        read(lunt,*)xyz(:,i)
      enddo
      close(lunt)

      call trajectory_to_bfield(nstep,emasskg1,echarge1,dmyenergy,
     &  xyz,bxyz,istat)

      open(newunit=lunt,file=chfilebx)
      write(lunt,*)icode,trim(code)
      write(lunt,*)'1. 1.'
      write(lunt,*)nstep
      do i=1,nstep
        write(lunt,*)xyz(1,i),bxyz(1,i)
      enddo
      close(lunt)

      open(newunit=lunt,file=chfileby)
      write(lunt,*)icode,trim(code)
      write(lunt,*)'1. 1.'
      write(lunt,*)nstep
      do i=1,nstep
        write(lunt,*)xyz(1,i),bxyz(2,i)
      enddo
      close(lunt)

      open(newunit=lunt,file=chfilebz)
      write(lunt,*)icode,trim(code)
      write(lunt,*)'1. 1.'
      write(lunt,*)nstep
      do i=1,nstep
        write(lunt,*)xyz(1,i),bxyz(3,i)
      enddo
      close(lunt)

      deallocate(xyz,bxyz)

      return
      end
