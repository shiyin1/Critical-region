subroutine initial(Nflow,yflow,kappa)
!make the initialization

  implicit none
  integer Nflow
  real(16) yflow(Nflow) !sigma is the fixed expansion point at UV
  integer N_str(4) !store the structure of functions of ODE
  integer Nv,Nh,Nz,Nck
  real(16) lam0,lam1,lam2,lam3,lam4,lam5,lam6,lam7
  real(16) h
  real(16) Zphi,Zpsi
  real(16) c,kappa
  real(16) pi,hc
  parameter(pi=3.141592653589793238462643383279Q+0)
  parameter(hc=197.33Q+0)
  real(16) lambda,nu,cl

  common /strucFun/ N_str
  common /break/ cl


  Nv=N_str(1)
  Nh=N_str(2)
  Nz=N_str(3)
  Nck=N_str(4)

  lambda=44.Q+0
  nu=(0.Q+0/hc)**2

  h=6.5Q+0
!expansion coefficients of Yukawa coupling

  c=1.Q-5!1.75Q-3*(1000.Q+0/hc)**3  !explicit chiral symmetry breaking term, in unit of fm**(-3)

  !kappa=0.Q+00

  lam0=(kappa**2*lambda)/2.Q+0+nu*kappa
  lam1=kappa*lambda+nu
  lam2=lambda
  lam3=0.Q+0
  lam4=0.Q+0
  lam5=0.Q+0
  lam6=0.Q+0
  lam7=0.Q+0
!expansion coefficients of effective potential V

  Zphi=1.Q+0 !meson wave function renormalization
  Zpsi=1.Q+0 !quark wave function renormalization


  yflow(1)=lam1
  yflow(2)=lam2
  yflow(3)=lam3
  yflow(4)=lam4
  yflow(5)=lam5
  yflow(6)=lam6
  yflow(7)=lam7
  yflow(Nv+1)=lam0
  yflow((Nv+1)+1)=h
  yflow((Nv+1)+(Nh+1)+1)=Zphi
  yflow((Nv+1)+(Nh+1)+2)=Zpsi
  yflow((Nv+1)+(Nh+1)+Nz+1)=c
  yflow((Nv+1)+(Nh+1)+Nz+2)=kappa

end





