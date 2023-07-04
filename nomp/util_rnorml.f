*CMZ :  4.00/07 05/05/2020  21.22.42  by  Michael Scheer
*CMZ :  4.00/05 29/11/2019  12.01.22  by  Michael Scheer
*-- Author :    Michael Scheer   28/08/2014

c Based on RNORML of CERNLIB, call to RANMAR substituted

cmsh Generated with: cpp -E -DCERNLIB_DOUBLE -DCERNLIB_UNIX rnorml.F

# 1 "rnorml.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/usr/include/stdc-predef.h" 1 3 4
# 1 "<command-line>" 2
# 1 "rnorml.F"
*
* $Id: rnorml.F,v 1.1.1.1 1996/04/01 15:02:55 mclareni Exp $
*
* $Log: rnorml.F,v $
* Revision 1.1.1.1 1996/04/01 15:02:55 mclareni
* Mathlib gen
*
*
# 1 "/usr/include/gen/pilot.h" 1 3 4
# 10 "rnorml.F" 2
      SUBROUTINE util_RNORML(DEVIAS,NDEV)
C Generator of a vector of independent Gaussian-distributed
C (pseudo-)random numbers, of mean zero and variance one,
C making use of a uniform pseudo-random generator (RANMAR).
C The algorithm for converting uniform numbers to Gaussian
C is that of "Ratio of Uniforms with Quadratic Bounds." The
C method is in principle exact (apart from rounding errors),
C and is based on the variant published by Joseph Leva in
C ACM TOMS vol. 18(1992), page 449 for the method and 454 for
C the Fortran algorithm (ACM No. 712).
C It requires at least 2 and on average 2.74 uniform deviates
C per Gaussian (normal) deviate.
C WARNING -- The uniform generator should not produce exact zeroes,
C since the pair (0.0, 0.5) provokes a floating point exception.
      SAVE S, T, A, B, R1, R2
      DIMENSION U(2), DEVIAS(*)
      DATA S, T, A, B / 0.449871, -0.386595, 0.19600, 0.25472/
      DATA R1, R2/ 0.27597, 0.27846/
C generate pair of uniform deviates
      DO 200 IDEV = 1, NDEV
cmsh   50 CALL RANMAR(U,2)
   50 CALL util_random(2,U)
      V = 1.7156 * (U(2) - 0.5)
      X = U(1) - S
      Y = ABS(V) - T
      Q = X**2 + Y*(A*Y - B*X)
C accept P if inside inner ellipse
      IF (Q .LT. R1) GO TO 100
C reject P if outside outer ellipse
      IF (Q .GT. R2) GO TO 50
C reject P if outside acceptance region
      IF (V**2 .GT. -4.0 *ALOG(U(1)) *U(1)**2) GO TO 50
C ratio of P's coordinates is normal deviate
  100 DEVIAT = V/U(1)
cmsh  200 DEVIAS(IDEV) = DEVIAT
      DEVIAS(IDEV) = DEVIAT
 200  continue
      RETURN
      END
