*CMZ : 00.00/19 11/08/2015  16.40.22  by  Michael Scheer
*CMZ : 00.00/18 11/08/2015  14.14.14  by  Michael Scheer
*CMZ : 00.00/17 26/07/2015  10.35.14  by  Michael Scheer
*CMZ : 00.00/16 24/07/2015  13.16.38  by  Michael Scheer
*-- Author :    Michael Scheer   21/07/2015
      subroutine util_spline_fit_solve(npar,p,v,ifail)

c +PATCH,//UTIL/FOR
c +DECK,util_spline_fit_solve.

      implicit none

      integer npar,ifail

      double precision
     &  p(npar*6),v(npar),ws(6,6),pt(6,6),vt(6)
     &  ,pt4(4,4)

      integer nival,ival,ip,iv,it,i,k,idebug,i1,i2,ipar

      idebug=0
      if (ifail.ne.0) idebug=ifail

      ifail=-1
      nival=(npar-2)/2

      if (idebug.ne.0) then

        print*,"util_spline_fit: p"
        do i=1,npar
          i1=(i-1)*6+1
          i2=i1+5
          print*,i,i1,sngl(p(i1:i2))
        enddo

        print*,"util_spline_fit: v"
        do i=1,npar
          print*,i,v(i)
        enddo
        print*,"-------------"
      endif

      if (nival.eq.1) then
        vt=v
        pt4=0.0d0
        do i=1,npar
          if (i.lt.3) then
            do k=1,4
              pt4(i,k)=p(k+(i-1)*6)
            enddo
          else if (i.gt.npar-2) then
            do k=1,4
              pt4(i,npar-4+k)=p(2+k+(i-1)*6)
            enddo
          else
            do k=1,6
              if (mod(i,2).eq.1) pt4(i,i-3+k)=p(k+(i-1)*6)
              if (mod(i,2).eq.0) pt4(i,i-4+k)=p(k+(i-1)*6)
            enddo
          endif
        enddo

        if (idebug.ne.0) then
          print*,"-------------"
          print*,"Matrix to solve:"
          do i=1,4
            print*,i,sngl(pt4(i,1:4))
          enddo
          print*,"-------------"
        endif

        call util_lineqnsys_gauss(pt4,npar,vt,ws,ifail)
        v(1:npar)=vt(1:npar)

        if (idebug.ne.0) then
          print*,"-------------"
          print*," Solution:"
          do ipar=1,npar
            print*,ipar,sngl(vt(ipar))
          enddo
          print*,"-------------"
        endif

        return

      else if (nival.eq.2) then
        vt=v
        pt=0.0d0
        do i=1,npar
          if (i.lt.3) then
            do k=1,4
              pt(i,k)=p(k+(i-1)*6)
            enddo
          else if (i.gt.npar-2) then
            do k=1,4
              pt(i,npar-4+k)=p(2+k+(i-1)*6)
            enddo
          else
            do k=1,6
              if (mod(i,2).eq.1) pt(i,i-3+k)=p(k+(i-1)*6)
              if (mod(i,2).eq.0) pt(i,i-4+k)=p(k+(i-1)*6)
            enddo
          endif
        enddo
        if (idebug.ne.0) then
          do i=1,npar
            print*,sngl(pt(i,1:npar)),sngl(vt(i))
          enddo
        endif
        call util_lineqnsys_gauss(pt,npar,vt,ws,ifail)
        v(1:npar)=vt(1:npar)
        if (idebug.ne.0) print*,"v:",v
        return
      endif

      v(1)=v(1)/p(1)
      p(1:4)=p(1:4)/p(1)
      v(2)=v(2)/p(7)
      p(7:10)=p(7:10)/p(7)
      v(2)=v(2)-v(1)
      p(7:10)=p(7:10)-p(1:4)
      v(2)=v(2)/p(8)
      p(8:10)=p(8:10)/p(8)

      do ival=2,nival

        if (ival.eq.2) then
c Erste Reihe

          ip=(ival-1)*12+1
          iv=(ival-1)*2+1
          v(iv)=v(iv)/p(ip)
          p(ip:ip+5)=p(ip:ip+5)/p(ip)

          v(iv)=v(iv)-v(iv-2)
          p(ip:ip+5)=p(ip:ip+5)-p(ip-12:ip-12+5)

          ip=ip+1
          v(iv)=v(iv)/p(ip)
          p(ip:ip+4)=p(ip:ip+4)/p(ip)

          v(iv)=v(iv)-v(iv-1)
          p(ip:ip+4)=p(ip:ip+4)-p(ip-6:ip-6+4)

          ip=ip+1
          v(iv)=v(iv)/p(ip)
          p(ip:ip+3)=p(ip:ip+3)/p(ip)

c Zweite Reihe

          ip=(ival-1)*12+1+6
          iv=(ival-1)*2+1+1
          v(iv)=v(iv)/p(ip)
          p(ip:ip+5)=p(ip:ip+5)/p(ip)

          v(iv)=v(iv)-v(iv-2-1)
          p(ip:ip+5)=p(ip:ip+5)-p(ip-12-6:ip-12+5-6)

          ip=ip+1
          v(iv)=v(iv)/p(ip)
          p(ip:ip+4)=p(ip:ip+4)/p(ip)

          v(iv)=v(iv)-v(iv-1-1)
          p(ip:ip+4)=p(ip:ip+4)-p(ip-6-6:ip-6+4-6)

          ip=ip+1
          v(iv)=v(iv)/p(ip)
          p(ip:ip+3)=p(ip:ip+3)/p(ip)

          v(iv)=v(iv)-v(iv-1)
          p(ip:ip+4)=p(ip:ip+4)-p(ip-6:ip-6+4)

          ip=ip+1
          v(iv)=v(iv)/p(ip)
          p(ip:ip+2)=p(ip:ip+2)/p(ip)

        else !if (ival.eq.3) then

c Erste Reihe

          ip=(ival-1)*12+1
          iv=(ival-1)*2+1
          v(iv)=v(iv)/p(ip)
          p(ip:ip+5)=p(ip:ip+5)/p(ip)

          v(iv)=v(iv)-v(iv-2)
          p(ip:ip+5-2)=p(ip:ip+5-2)-p(ip-12+2:ip-12+5-2)

          ip=ip+1
c          print*,ip,p(ip:ip+4)
          v(iv)=v(iv)/p(ip)
          p(ip:ip+4)=p(ip:ip+4)/p(ip)

          v(iv)=v(iv)-v(iv-1)
          p(ip:ip+2)=p(ip:ip+2)-p(ip-6+2:ip-6+4+2)

          ip=ip+1
c          print*,ip,p(ip:ip+4)
          v(iv)=v(iv)/p(ip)
          p(ip:ip+3)=p(ip:ip+3)/p(ip)

c Zweite Reihe

          ip=(ival-1)*12+1+6
          iv=(ival-1)*2+1+1
          v(iv)=v(iv)/p(ip)
          p(ip:ip+5)=p(ip:ip+5)/p(ip)

          v(iv)=v(iv)-v(iv-2-1)
          p(ip:ip+5-2)=p(ip:ip+5-2)-p(ip-18+2:ip-18+5-2)

          ip=ip+1
          v(iv)=v(iv)/p(ip)
          p(ip:ip+4)=p(ip:ip+4)/p(ip)

c          print*,iv,v(iv-2),v(iv)
c          print*,ip,p(ip:ip+1),p(ip-12+2:ip-12+2+1)

          v(iv)=v(iv)-v(iv-2)
          p(ip:ip+4)=p(ip:ip+4)-p(ip-12+2:ip-12+2+4)

          ip=ip+1
          v(iv)=v(iv)/p(ip)
          p(ip:ip+3)=p(ip:ip+3)/p(ip)

          v(iv)=v(iv)-v(iv-1)
          p(ip:ip+4)=p(ip:ip+4)-p(ip-6:ip-6+4)

          ip=ip+1
          v(iv)=v(iv)/p(ip)
          p(ip:ip+2)=p(ip:ip+2)/p(ip)

        endif !ival
      enddo

c Vorletzte Reihe

      ip=(npar-2)*6+3
      iv=npar-1

      v(iv)=v(iv)/p(ip)
      p(ip:ip+3)=p(ip:ip+3)/p(ip)

      v(iv)=v(iv)-v(iv-2)
      p(ip:ip+3)=p(ip:ip+3)-p(ip-12:ip-12+3)

      ip=ip+1
      v(iv)=v(iv)/p(ip)
      p(ip:ip+2)=p(ip:ip+2)/p(ip)

      v(iv)=v(iv)-v(iv-1)
      p(ip:ip+2)=p(ip:ip+2)-p(ip-6:ip-6+2)

      ip=ip+1
      v(iv)=v(iv)/p(ip)
      p(ip:ip+1)=p(ip:ip+1)/p(ip)

c Letzte Reihe

      ip=(npar-1)*6+3
      iv=npar

      v(iv)=v(iv)/p(ip)
      p(ip:ip+3)=p(ip:ip+3)/p(ip)

      v(iv)=v(iv)-v(iv-3)
      p(ip:ip+3)=p(ip:ip+3)-p(ip-18:ip-18+3)

      ip=ip+1
      v(iv)=v(iv)/p(ip)
      p(ip:ip+2)=p(ip:ip+2)/p(ip)

      v(iv)=v(iv)-v(iv-2)
      p(ip:ip+2)=p(ip:ip+2)-p(ip-12:ip-12+2)

      ip=ip+1
      v(iv)=v(iv)/p(ip)
      p(ip:ip+1)=p(ip:ip+1)/p(ip)

      v(iv)=v(iv)-v(iv-1)
      p(ip:ip+2)=p(ip:ip+2)-p(ip-6:ip-6+2)

      ip=ip+1
      v(iv)=v(iv)/p(ip)
      p(ip)=1.0d0

c Jetzt ist die Matrix eine Dreiecksmatrix

      ip=(npar-2)*6+5
      iv=npar-1

      v(iv)=v(iv)/p(ip+1)
      p(ip:ip+1)=p(ip:ip+1)/p(ip+1)

      v(iv)=v(iv)-v(iv+1)
      p(ip:ip+2)=p(ip:ip+2)-p(ip+6:ip+6+2)

      v(iv)=v(iv)/p(ip)
      p(ip)=1.0d0

      do ival=nival,1,-1

        if (ival.ge.2) then

          do it=1,2

            ip=(npar-2-(nival-ival)*2-it)*6+2+3-it
            iv=npar-1-(nival-ival)*2-it
            do k=1,1+it
              v(iv)=v(iv)-v(iv+k)*p(ip+k)
              p(ip+k)=0.0d0
            enddo
          enddo
        else

          do it=1,2

            ip=(npar-2-(nival-ival)*2-it)*6+2+1-it
            iv=npar-1-(nival-ival)*2-it
            do k=1,2+it
              v(iv)=v(iv)-v(iv+k)*p(ip+k)
              p(ip+k)=0.0d0
            enddo

          enddo

        endif !ival

      enddo

      ifail=0

      return
      end
