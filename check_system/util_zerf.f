*CMZ : 00.00/02 10/04/2004  15.46.21  by  Michael Scheer
*-- Author :    Michael Scheer   10/04/2004
      function util_zerf(z)

      external function wwerf

      complex*16 util_zerf,z,wwerf,zone,zi
      parameter (zone=(1.0d0,0.0d0))
      parameter (zi=(0.0d0,1.0d0))

      util_zerf=zone-wwerf(zi*z)*exp(-z*z)

      return
      end
