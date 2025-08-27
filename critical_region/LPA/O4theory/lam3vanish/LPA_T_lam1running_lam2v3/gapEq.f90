subroutine gapEq(n, x, fvec)

  implicit none

  integer n
  real(16) x(n), fvec(n)                     !x(1)
  real(16) rho,sigma
  real(16) lam1,lam2,lam3,lam4,lam5,lam6,lam7,c,kappa

  common /gapPara/ lam1,lam2,lam3,lam4,lam5,lam6,lam7,c,kappa


  sigma=x(1)

  !fvec(1)=lam1 - c/(Sqrt(2.Q+0)*Sqrt(rho)) + lam2*(-kappa + rho) + &
  !        (lam3*(-kappa + rho)**2)/2.Q+0
  fvec(1)=lam1*sigma + lam2*sigma*(-kappa + sigma**2/2) + (lam3*sigma*(-kappa + sigma**2/2)**2)/2.Q+0

end
