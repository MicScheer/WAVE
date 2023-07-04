*CMZ : 00.00/15 09/10/2012  13.45.25  by  Michael Scheer
*-- Author :    Michael Scheer   05/10/2012
      subroutine util_beta_function_drift(
     &  s0,beta0,gamma0,
     &  s1,beta1,betap1,alpha1,gamma1,phase1,
     &  s2,beta2,betap2,alpha2,gamma2,phase2)

      real*8 s1,beta1,betap1,alpha1,gamma1,phase1,phase2,
     &  s2,beta2,betap2,alpha2,gamma2,beta0,gamma0,s0

c Calculates beta(s2) etc. from beta(s1) and betap(s1)
c beta(s)=beta0+(s-s0)**2/beta(0)

      alpha1=-betap1/2.0d0
      gamma1=(1.0d0+alpha1**2)/beta1
      s0=s1+alpha1/gamma1
      beta0=1.0d0/gamma1
      gamma0=1.0d0/beta0
      phase1=atan((s1-s0)/beta0)

      beta2=beta0+(s2-s0)**2/beta0
      betap2=2.0d0*(s2-s0)/beta0
      alpha2=-(s2-s0)/beta0
      gamma2=(1.0d0+alpha2**2)/beta2
      phase2=atan((s2-s0)/beta0)

      return
      end
