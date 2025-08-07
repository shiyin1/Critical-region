subroutine initial(a)
!make the initialization

  implicit none

  integer(8) gridnum,i,anum
  !number of phi
  integer Nflow,Vmax,maxim
  !number of flow equations
  parameter (Nflow=800,Vmax=400)
  real(16) drho,rho(Vmax),yflow(Nflow),a(Vmax),f(Vmax),rhomin,rhomax
  real(16) kappa(Vmax)
  real(16) lami1,lami2,lami3
  real(16) h,cbar
  !Yukawa coupling and explicit chiral symmetry breaking 

  common /gridrho/ gridnum,anum,rho,rhomin,rhomax,maxim
  common /iniconst/ h,cbar

  h=6.5Q+00
  !Yukawa coupling
 cbar=1.73e06
!  cbar=121.Q+00**3.Q+00
!   c=1.7D-3*(1000.)**3
!  c=0.Q+00
  !explicit chiral symmetry breaking term
  !in unit of MeV^3

  lami1=495.Q+00**2
  lami2=5.Q+00
  lami3=3.33123Q+00*10Q-82 !anum=21

  kappa=rho
  yflow=0.Q+00
  yflow(1:gridnum)=lami1*kappa+lami2*kappa**2/2.Q+00+lami3*kappa**20.Q+00
  yflow((1+gridnum):(2*gridnum))=lami2*kappa
  yflow((2*gridnum+1):(3*gridnum))=lami2
  f=yflow(1:gridnum)
!  write(*,*) f
!  stop
  call chebft(a,gridnum,anum,f)

end
