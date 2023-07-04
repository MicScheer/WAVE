*CMZ : 00.00/02 13/04/2004  08.08.28  by  Michael Scheer
*-- Author :    Michael Scheer   10/04/2004
      function util_int_p0expiwp2(a0,p0,p1,p2,w,tl,th)

c returns int(a0*exp(I*w*(p0+p1*t+p2*t**2t)),t=tl...th)

      external function util_zerf

      complex*16 util_int_p0expiwp2,util_zerf,zi,z,ewp0,rwp2i,zh,zl,z1

      double precision rpi,pi,a0,p0,p1,p2,w,tl,th

      parameter (pi=3.141592653589793D0)
      parameter (rpi=1.77245385090552D0)
      parameter (zi=(0.0d0,1.0d0))

      if (p2.eq.0.0d0) then

        if (p1.eq.0.0d0) then
          util_int_p0expiwp2=a0*exp(zi*w*p0)*(th-tl)
        else
          zh=exp(zi*th*w*p1)
          zl=exp(zi*tl*w*p1)
          util_int_p0expiwp2=-zi/w*a0/p1*exp(zi*w*p0)*(zh-zl)
        endif

      else

        rwp2i=sqrt(-zi*w*p2)
c        ewp0=exp(zi*w*p0)

        z1=a0/2.0d0*rpi/rwp2i*exp(zi*w*p0)*exp(-zi*w*p1*p1/p2/4.0d0)

        zh=util_zerf((p1/2.0d0+th*p2)*rwp2i/p2)
        zl=util_zerf((p1/2.0d0+tl*p2)*rwp2i/p2)

        util_int_p0expiwp2=z1*(zh-zl)

      endif !p2.eq.0

      return
      end
