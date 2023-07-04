*CMZ :  3.03/02 22/03/2016  13.20.03  by  Michael Scheer
*CMZ :  3.00/00 11/03/2013  15.12.11  by  Michael Scheer
*CMZ :  2.70/12 01/03/2013  16.28.24  by  Michael Scheer
*CMZ :  2.66/03 29/04/2010  11.46.31  by  Michael Scheer
*CMZ :  2.51/02 08/10/2009  09.58.11  by  Michael Scheer
*CMZ :  2.16/08 23/10/2000  14.22.46  by  Michael Scheer
*CMZ :  2.13/10 14/04/2000  14.26.49  by  Michael Scheer
*CMZ :  2.13/04 21/01/2000  14.54.46  by  Michael Scheer
*CMZ :  2.13/03 11/01/2000  18.22.28  by  Michael Scheer
*CMZ : 00.00/04 29/04/94  17.56.38  by  Michael Scheer
*CMZ : 00.00/00 28/04/94  16.12.16  by  Michael Scheer
*-- Author : Michael Scheer
      SUBROUTINE POWFOLD
*KEEP,gplhint.
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
*KEND.

*KEEP,spectf90u.
      include 'spectf90u.cmn'
*KEEP,sourcef90u.
      include 'sourcef90u.cmn'
*KEEP,observf90u.
      include 'observf90u.cmn'
*KEND.

C--- CALCULATES FOLDING OF POWER-DENSITY IN PINHOLE WITH ELECTRON PHASE SPACE
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
*KEEP,wfoldf90.
      include 'wfoldf90.cmn'
*KEND.

      double precision, dimension(:,:), allocatable :: poww,powwf
      double precision, dimension(:), allocatable :: powz,powy,powfz,powfy
     &  ,w1,w2,w3,w4,coef

      integer iz,iy,iobsv,isour,i

      allocate(poww(nobsvz,nobsvy))
      allocate(powwf(nobsvz,nobsvy))
      allocate(w1(max(nobsvz,nobsvy)))
      allocate(w2(max(nobsvz,nobsvy)))
      allocate(w3(max(nobsvz,nobsvy)))
      allocate(w4(max(nobsvz,nobsvy)))
      allocate(coef(max(nobsvz,nobsvy)))

      allocate(powz(nobsvz))
      allocate(powy(nobsvy))
      allocate(powfz(nobsvz))
      allocate(powfy(nobsvy))

      specpowf=0.0d0
      specpowtf=0.0d0

      do isour=1,nsource

        iobsv=0
        do iy=1,nobsvy
          do iz=1,nobsvz
            iobsv=iobsv+1
            poww(iz,iy)=SPECPOW(ISOUR+(IOBSV-1)*NSOURCE)
          enddo
        enddo

       do iy=(nobsvy-mobsvy)/2+1,(nobsvy-mobsvy)/2+mobsvy
         powz=poww(1:nobsvz,iy)
         if (dgsigz(isour).gt.0.0d0.and.wsigz(isour).gt.0.0d0) then
           call util_fold_function_gauss(
     &        nobsvz,obsvz,powz,wsigz(isour),dgsigz(1),powfz,
     &        coef,w1,w2,w3,w4)
         else
           write(lungfo,*)"*** Warning in powfold: Dgsigz or wsigz is zero for source",isour
           write(6,*)"*** Warning in powfold: Dgsigz or wsigz is zero for source",isour
           powfz=powz
         endif
         powwf(1:nobsvz,iy)=powfz
       enddo !iy

       do iz=1,nobsvz
         powy=powwf(iz,1:nobsvy)
         if (dgsigy(isour).gt.0.0d0.and.wsigy(isour).gt.0.0d0) then
           call util_fold_function_gauss(
     &        nobsvy,obsvy,powy,wsigy(isour),dgsigy(1),powfy,
     &        coef,w1,w2,w3,w4)
         else
           write(lungfo,*)"*** Warning in powfold: Dgsigy or wsigy is zero for source",isour
           write(6,*)"*** Warning in powfold: Dgsigy or wsigy is zero for source",isour
           powfy=powy
         endif
         powwf(iz,1:nobsvy)=powfy
       enddo !iz

       do iy=(nobsvy-mobsvy)/2+1,(nobsvy-mobsvy)/2+mobsvy
         do iz=(nobsvz-mobsvz)/2+1,(nobsvz-mobsvz)/2+mobsvz
           iobsv=nobsvz*(iy-1)+iz
           SPECPOWf(ISOUR+(IOBSV-1)*NSOURCE)=powwf(iz,iy)
           specpowtf(iobsv)=specpowtf(iobsv)+specpowf(isour+nsource*(iobsv-1))
         enddo
       enddo

      enddo !isour=1,nsource

      return
      end
