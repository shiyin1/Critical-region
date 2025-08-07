subroutine initial(yflow,kappa0)
!make the initialization

  implicit none

  integer gridnum,i
  !number of phi
  integer Nflow,Vmax,maxim
  !number of flow equations
  parameter (Nflow=600,Vmax=600)
  real(16) drho,rho(Vmax),yflow(Nflow),f(Vmax)
  real(16) rhomin,rhomax,rhoacc
  real(16) kappa(Vmax),kappa0
  real(16) lam1,lam2,lam3,lami
  real(16) Z_pi,Z_sigma,Zt_pi,Zt_sigma
  real(16) h,cbar
  !Yukawa coupling and explicit chiral symmetry breaking 

  common /gridrho/ gridnum,drho,rho,rhomin,rhomax,rhoacc
  common /iniconst/ h,cbar

  h=6.4Q+00
  !Yukawa coupling
!  cbar=0.Q+00
!  cbar=121.Q+00**3.Q+00
  !explicit chiral symmetry breaking term
  !in unit of MeV^3

  lam1=2.Q+00
  lam2=20.Q+00
  lam3=2.319Q+00*10Q-138
  Z_pi=1.Q+00
  Z_sigma=1.Q+00
  Zt_pi=1.Q+00
  Zt_sigma=1.Q+00

  kappa=rho
  yflow=0.Q+00
  yflow(1:gridnum)=lam1*(kappa-kappa0)
  yflow(gridnum+1)=Z_pi
  yflow(gridnum+2)=Zt_pi
!  write(*,*) f
!  stop


end
