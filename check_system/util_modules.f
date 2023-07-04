*CMZ :          28/07/2018  09.43.27  by  Michael Scheer
*-- Author :    Michael Scheer   28/07/2018
*KEEP,COHERENT.

      module coherentmod

! Bunchcharge nicht mehr verwendet, da zu kompliziert...

      integer nein,nbunch,nstep,ncoef,mode,nsigs

      double precision xlam,xlen,xlenmn,xlenmx,dxlen,fmean,frms
      double precision c(0:1000),xlenfou,bunchcharge

      namelist/coherentn/
     &  mode,ncoef,c,nein,nbunch,nstep,xlenfou,bunchcharge
     &  ,xlam,xlen,xlenmn,xlenmx,dxlen,fmean,frms

      end module

