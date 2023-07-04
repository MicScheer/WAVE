*CMZ : 00.00/02 13/04/2004  13.28.54  by  Michael Scheer
*-- Author :    Michael Scheer   10/04/2004
      function util_int_p3expiwp2(a0,a1,a2,a3,p0,p1,p2,w,tl,th)

      external function util_zerf

      complex*16 util_int_p3expiwp2,util_zerf,zi,
     &  ewp0,rwp2i,ziw,
     &  az0,az1,az2,az3,az3h,az3l,azh,azl,erfh,erfl,
     &  c1,c4,c8,c16,e1,e2,e2h,e2l,e3h,e3l,azhl

      double precision pi,a0,a1,a2,a3,p0,p1,p2,w,tl,th,t,
     &  rpi,wp1,wp2,thp1,tlp1,thp2,tlp2

      parameter (pi=3.141592653589793D0)
      parameter (rpi=1.77245385090552D0)
      parameter (zi=(0.0d0,1.0d0))

      wp1=w*p1
      wp2=w*p2
      thp1=th*p1
      thp2=th*p2
      tlp1=tl*p1
      tlp2=tl*p2
      ziw=zi*w

      ewp0=exp(ziw*p0)
      rwp2i=sqrt(-zi*wp2)

      if (p2.eq.0.0d0) then

        if (p1.eq.0.0d0) then

          util_int_p3expiwp2=ewp0*
     &      (
     &      th*
     &      (a0+th*
     &      (a1/2.0d0+th*
     &      (a2/3.0d0+a3/4.0d0*th)))-
     &      tl*
     &      (a0+tl*
     &      (a1/2.0d0+tl*
     &      (a2/3.0d0+a3/4.0d0*tl)))
     &     )


        else !(p1.eq.0.0d0)

          e1=exp(zi*th*wp1)/ziw
          e2=exp(zi*tl*wp1)/ziw

c az0:=int(a0*exp(i*w*p),t);
          az0=a0/p1*(e1-e2)

c az1:=int(a1*t*exp(i*w*p),t);
          az1=a1/p1*
     &      (
     &      (th - 1.0d0/ziw/p1)*e1-
     &      (tl - 1.0d0/ziw/p1)*e2
     &      )

c az2:=int(a2*t^2*exp(i*w*p),t);
          az2=a2/p1*
     &      (
     &      e1*(2.0d0/ziw**2/p1**2 - 2.0d0/ziw*th/p1+th**2) -
     &      e2*(2.0d0/ziw**2/p1**2 - 2.0d0/ziw*tl/p1+tl**2)
     &      )

c az3:=int(a3*t^3*exp(i*w*p),t);
          az3=a3/p1*(
     &      e1*
     &      (6.0d0/ziw**2*th/p1**2 -
     &      6.0d0/ziw**3/p1**3 + th**3 -
     &      3.0d0/ziw*th**2/p1) -
     &      e2*
     &      (6.0d0/ziw**2*tl/p1**2 -
     &      6.0d0/ziw**3/p1**3 + tl**3 -
     &      3.0d0/ziw*tl**2/p1))

          util_int_p3expiwp2=ewp0*(az0+az1+az2+az3)

        endif !(p1.eq.0.0d0)

      else !(p2.eq.0.0d0)

        erfh=util_zerf((p1/2.0d0+ thp2)/p2*rwp2i)
        erfl=util_zerf((p1/2.0d0+ tlp2)/p2*rwp2i)

        e2h=rwp2i*exp(ziw*th*(p1+thp2))
        e2l=rwp2i*exp(ziw*tl*(p1+tlp2))

        e1=exp(-1.0d0/4.0d0*zi*wp1*p1/p2)

        c1=rpi*e1
        c4=1.0d0/(4.0d0*wp2*rwp2i)
        c8=c4/(2.0d0*p2)
        c16=c8/(2.0d0*wp2)

        e3h=-4.0d0*a2*c8*thp2*zi*e2h
        e3l=-4.0d0*a2*c8*tlp2*zi*e2l

c  az3:=int(t**3*exp(zi*w*ph),t);
        az3h=a3*c16*
     &    (
     &    +8.0d0*p2*e2h
     &    +(4.0d0*p1*thp2
     &    -2.0d0*p1**2
     &    -8.0d0*thp2*thp2)*ziw*e2h
     &    -(6.0d0*zi*p2+wp1*p1)*wp1*c1*erfh
     &    )

        az3l=a3*c16*
     &    (
     &    +8.0d0*p2*e2l
     &    +(4.0d0*p1*tlp2
     &    -2.0d0*p1**2
     &    -8.0d0*tlp2*tlp2)*ziw*e2l
     &    -(6.0d0*zi*p2+wp1*p1)*wp1*c1*erfl
     &    )

        azhl=
     &    (a0/2.0d0*c1/rwp2i+(2.0d0*zi*p2+wp1*p1)*a2*c8*c1-a1*c4*c1*wp1)
     &                                                     * (erfh-erfl)
     &    +(2.0d0*a2*c8*p1*zi-2.0d0*a1*c4*zi)              * (e2h-e2l)
     &    +                                                  (e3h-e3l)
     &    +                                                  (az3h-az3l)

        util_int_p3expiwp2=ewp0*azhl

      endif !(p2.eq.0.0d0)

      return
      end
