subroutine initial(a)
!make the initialization

  implicit none

  integer(8) gridnum,i
  !number of phi
  integer Nflow,Vmax
  !number of flow equations
  parameter (Nflow=800,Vmax=200)
  real(16) drho,rho(Vmax),yflow(Nflow),a(Vmax),f(Vmax)
  real(16) kappa(Vmax)
  real(16) lami1,lami2
  real(16) h,cbar
  !Yukawa coupling and explicit chiral symmetry breaking 

  common /gridrho/ gridnum,drho,rho
  common /iniconst/ h,cbar

  h=6.50Q+00
  !Yukawa coupling
  cbar=1.75Q+6
!  cbar=1.Q+00
!  cbar=121.Q+00**3.Q+00
  !explicit chiral symmetry breaking term
  !in unit of MeV^3

  lami1=0.Q+00
  lami2=45.Q+00

  kappa=rho
  yflow=0.Q+00
  yflow(1:gridnum)=lami1*kappa+lami2*kappa**2/2.Q+00
  yflow((1+gridnum):(2*gridnum))=lami2*kappa
  yflow((2*gridnum+1):(3*gridnum))=lami2
  f=yflow(1:gridnum)

  call chebft(a,gridnum,gridnum,f)

end
