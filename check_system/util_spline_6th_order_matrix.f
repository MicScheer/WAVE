*CMZ : 00.00/10 28/12/2010  13.38.45  by  Michael Scheer
*CMZ : 00.00/09 21/12/2010  17.21.30  by  Michael Scheer
*CMZ :          22/12/2010  10.11.25  by  Michael Scheer
*CMZ : 00.00/09 21/12/2010  17.21.30  by  Michael Scheer
*-- Author :    Michael Scheer   21/12/2010
      subroutine util_spline_6th_order_matrix(nord,amat,vec,ifail)

c +PATCH,//UTIL/FOR
c +DECK,util_spline_6th_order_matrix.

c The data for the first two lines are always expected in amat(1:4,1:2)
c The data for the last two lines are always expected in amat(3:6,nord-1:nord)

      integer nord,ifail,i,k1,k2,k3,k4,k5,k6

      double precision amat(6,nord),vec(nord)

      if (nord.lt.4) then
        ifail=-1
        return
      endif

      if (mod(nord,2).ne.0) then
        ifail=-1
        return
      endif

      do i=1,(nord-2)/4

        k1=i
        k2=k1+1
        k3=k2+1
        k4=k3+1
        k5=k4+1
        k6=k5+1

        vec(k1)=vec(k1)/amat(1,k1)
        amat(1:4,1)=amat(1:4,1)/amat(1,k1)

        if (amat(1,k2).ne.0.0d0) then
          vec(k2)=vec(k2)/amat(1,k2)
          amat(1:4,k2)=amat(1:4,k2)/amat(1,k2)
          vec(k2)=vec(k2)-vec(k1)
          amat(2:4,k2)=amat(2:4,k2)-amat(2:4,k1)
          amat(1,k2)=0.0d0
        endif
        vec(k2)=vec(k2)/amat(2,k2)
        amat(2:4,k2)=amat(2:4,k2)/amat(2,k2)

        if (amat(1,k3).ne.0.0d0) then
          vec(k3)=vec(k3)/amat(1,k3)
          amat(1:6,k3)=amat(1:6,k3)/amat(1,k3)
          vec(k3)=vec(k3)-vec(k1)
          amat(2:4,k3)=amat(2:4,k3)-amat(2:4,k1)
          amat(1,k3)=0.0d0
        endif

        if (amat(2,k3).ne.0.0d0) then
          vec(k3)=vec(k3)/amat(2,k3)
          amat(2:6,k3)=amat(2:6,k3)/amat(2,k3)
          vec(k3)=vec(k3)-vec(k2)
          amat(3:4,k3)=amat(3:4,k3)-amat(3:4,k2)
          amat(2,k3)=0.0d0
        endif
        vec(k3)=vec(k3)/amat(3,k3)
        amat(3:6,k3)=amat(3:6,k3)/amat(3,k3)

        if (amat(1,k4).ne.0.0d0) then
          vec(k4)=vec(k4)/amat(1,k4)
          amat(1:6,k4)=amat(1:6,k4)/amat(1,k4)
          vec(k4)=vec(k4)-vec(k1)
          amat(2:4,k4)=amat(2:4,k4)-amat(2:4,k1)
          amat(1,k4)=0.0d0
        endif

        if (amat(2,k4).ne.0.0d0) then
          vec(k4)=vec(k4)/amat(2,k4)
          amat(2:6,k4)=amat(2:6,k4)/amat(2,k4)
          vec(k4)=vec(k4)-vec(k2)
          amat(3:4,k4)=amat(3:4,k4)-amat(3:4,k2)
          amat(2,k4)=0.0d0
        endif
        vec(k4)=vec(k4)/amat(3,k4)
        amat(3:6,k4)=amat(3:6,k4)/amat(3,k4)

      enddo !i=1,(nord-2)/4

c last four rows

      if (nord.eq.4) then
c We must shift the data from amat(1:4,nord-3:nord-2), since we get the first
c two of the last four rows there...
        amat(3:6,nord-3)=amat(1:4,nord-3)
        amat(1:2,nord-3)=0.0d0
        amat(3:6,nord-2)=amat(1:4,nord-2)
        amat(1:2,nord-2)=0.0d0
      endif

      i=nord-3

      k1=i
      k2=k1+1
      k3=k2+1
      k4=k3+1

      vec(k1)=vec(k1)/amat(3,k1)
      amat(3:6,k1)=amat(3:6,k1)/amat(3,k1)

      if (amat(3,k2).ne.0.0d0) then
        vec(k2)=vec(k2)/amat(3,k2)
        amat(3:6,k2)=amat(3:6,k2)/amat(3,k2)
        vec(k2)=vec(k2)-vec(k1)
        amat(3:6,k2)=amat(3:6,k2)-amat(3:6,k1)
      endif
      vec(k2)=vec(k2)/amat(4,k2)
      amat(4:6,k2)=amat(4:6,k2)/amat(4,k2)

      if (amat(3,k3).ne.0.0d0) then
        vec(k3)=vec(k3)/amat(3,k3)
        amat(3:6,k3)=amat(3:6,k3)/amat(3,k3)
        vec(k3)=vec(k3)-vec(k1)
        amat(3:6,k3)=amat(3:6,k3)-amat(3:6,k1)
      endif
      if (amat(4,k3).ne.0.0d0) then
        vec(k3)=vec(k3)/amat(4,k3)
        amat(4:6,k3)=amat(4:6,k3)/amat(4,k3)
        vec(k3)=vec(k3)-vec(k2)
        amat(4:6,k3)=amat(4:6,k3)-amat(4:6,k2)
      endif
      vec(k3)=vec(k3)/amat(5,k3)
      amat(5:6,k3)=amat(5:6,k3)/amat(5,k3)

      if (amat(3,k4).ne.0.0d0) then
        vec(k4)=vec(k4)/amat(3,k4)
        amat(3:6,k4)=amat(3:6,k4)/amat(3,k4)
        vec(k4)=vec(k4)-vec(k1)
        amat(3:6,k4)=amat(3:6,k4)-amat(3:6,k1)
      endif

      if (amat(4,k4).ne.0.0d0) then
        vec(k4)=vec(k4)/amat(4,k4)
        amat(4:6,k4)=amat(4:6,k4)/amat(4,k4)
        vec(k4)=vec(k4)-vec(k2)
        amat(4:6,k4)=amat(4:6,k4)-amat(4:6,k2)
      endif

      if (amat(5,k4).ne.0.0d0) then
        vec(k4)=vec(k4)/amat(5,k4)
        amat(5:6,k4)=amat(5:6,k4)/amat(5,k4)
        vec(k4)=vec(k4)-vec(k3)
        amat(5:6,k4)=amat(5:6,k4)-amat(5:6,k3)
      endif

      vec(k4)=vec(k4)/amat(6,k4)
      amat(6,k4)=1.0d0

c now we have a triangle

      if (amat(6,k3).ne.0.0d0) then
        vec(k3)=vec(k3)/amat(6,k3)
        amat(3:6,k3)=amat(3:6,k3)/amat(6,k3)
        vec(k3)=vec(k3)-vec(k4)
        amat(6,k3)=amat(6,k3)-amat(6,k4)
        vec(k3)=vec(k3)/amat(5,k3)
        amat(5,k3)=1.0d0
      endif

      if (amat(6,k2).ne.0.0d0) then
        vec(k2)=vec(k2)/amat(6,k2)
        amat(4:6,k2)=amat(4:6,k2)/amat(6,k2)
        vec(k2)=vec(k2)-vec(k4)
        amat(4:5,k2)=amat(4:5,k2)-amat(4:5,k4)
        amat(6,k2)=0.0d0
      endif

      if (amat(5,k2).ne.0.0d0) then
        vec(k2)=vec(k2)/amat(5,k2)
        amat(4:5,k2)=amat(4:5,k2)/amat(5,k2)
        vec(k2)=vec(k2)-vec(k3)
        amat(4,k2)=amat(4,k2)-amat(4,k3)
        amat(5,k2)=0.0d0
        vec(k2)=vec(k2)/amat(4,k2)
        amat(4,k2)=1.0d0
      endif

      if (amat(6,k1).ne.0.0d0) then
        vec(k1)=vec(k1)/amat(6,k1)
        amat(3:5,k1)=amat(3:5,k1)/amat(6,k1)
        vec(k1)=vec(k1)-vec(k4)
        amat(6,k1)=0.0d0
      endif

      if (amat(5,k1).ne.0.0d0) then
        vec(k1)=vec(k1)/amat(5,k1)
        amat(3:5,k1)=amat(3:5,k1)/amat(5,k1)
        vec(k1)=vec(k1)-vec(k3)
        amat(4,k1)=amat(4,k1)-amat(4,k3)
        amat(5,k1)=0.0d0
      endif

      if (amat(4,k1).ne.0.0d0) then
        vec(k1)=vec(k1)/amat(4,k1)
        amat(3,k1)=amat(3,k1)/amat(4,k1)
        vec(k1)=vec(k1)-vec(k2)
        amat(4,k1)=0.0d0
        vec(k1)=vec(k1)/amat(3,k1)
        amat(3,k1)=1.0d0
      endif

c all remaining upper rows

      do i=nord-5,1,-2

        k1=i
        k2=k1+1
        k3=k2+1
        k4=k3+1

        vec(k1)=vec(k1)/amat(1,k1)
        amat(1:4,k1)=amat(1:4,k1)/amat(1,k1)

        if (amat(1,k2).ne.0.0d0) then
          vec(k2)=vec(k2)/amat(1,k2)
          amat(1:4,k2)=amat(1:4,k2)/amat(1,k2)
          vec(k2)=vec(k2)-vec(k1)
          amat(1:4,k2)=amat(1:4,k2)-amat(1:4,k1)
        endif
        vec(k2)=vec(k2)/amat(2,k2)
        amat(2:4,k2)=amat(2:4,k2)/amat(2,k2)

c now we have a triangle

        if (amat(4,k2).ne.0.0d0) then
          vec(k2)=vec(k2)/amat(4,k2)
          amat(2:4,k2)=amat(2:4,k2)/amat(4,k2)
          vec(k2)=vec(k2)-vec(k4)
          amat(2:3,k2)=amat(2:3,k2)-amat(2:3,k4)
          amat(4,k2)=0.0d0
        endif

        if (amat(3,k2).ne.0.0d0) then
          vec(k2)=vec(k2)/amat(3,k2)
          amat(2:3,k2)=amat(2:3,k2)/amat(3,k2)
          vec(k2)=vec(k2)-vec(k3)
          amat(2,k2)=amat(2,k2)-amat(2,k3)
          amat(3,k2)=0.0d0
          vec(k2)=vec(k2)/amat(2,k2)
          amat(2,k2)=1.0d0
        endif

        if (amat(4,k1).ne.0.0d0) then
          vec(k1)=vec(k1)/amat(4,k1)
          amat(1:3,k1)=amat(1:3,k1)/amat(4,k1)
          vec(k1)=vec(k1)-vec(k4)
          amat(4,k1)=0.0d0
        endif

        if (amat(3,k1).ne.0.0d0) then
          vec(k1)=vec(k1)/amat(3,k1)
          amat(1:3,k1)=amat(1:3,k1)/amat(3,k1)
          vec(k1)=vec(k1)-vec(k3)
          amat(2,k1)=amat(2,k1)-amat(2,k3)
          amat(3,k1)=0.0d0
        endif

        if (amat(2,k1).ne.0.0d0) then
          vec(k1)=vec(k1)/amat(2,k1)
          amat(1,k1)=amat(1,k1)/amat(2,k1)
          vec(k1)=vec(k1)-vec(k2)
          amat(2,k1)=0.0d0
          vec(k1)=vec(k1)/amat(1,k1)
          amat(1,k1)=1.0d0
        endif

      enddo !i=nord-2,1,-4

      return
      end
