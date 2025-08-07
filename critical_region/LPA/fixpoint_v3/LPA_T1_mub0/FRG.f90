subroutine FRG(kappa_UV_i,kappa_UV,rho0,mPion,mSigma,mf,Vall,fpi,h,Zphi,Zpsi,c,kappa)
!This subroutine solve FRG flow equations with fixed expansion point, search for phyiscal point, find
!appropriate expansion point. Inputing a guess kappa_UV_i, outputing kappa_UV and other physical variables

  implicit none

  integer Nv,Nh
! Nv: order of Tylor expansion for effective potential V, Nh: order of Yukawa coupling h
  parameter(Nv=7)
  parameter(Nh=0)
  integer Nz !number of wave function renormalizations
  parameter(Nz=2)
  integer Nck
  parameter(Nck=2)
  integer Nflow !number of flow equations
  parameter(Nflow=(Nv+1)+(Nh+1)+Nz+Nck)
  real(16) yflow(Nflow)
!dependent variables in flow equations
  integer N_str(4) !store the structure of functions of ODE
  real(16) pi,hc
  parameter(pi=3.141592653589793238462643383279Q+0)
  parameter(hc=197.33Q+0)
  real(16) kappa_UV_i,kappa_UV
  real(16) T,mu
  real(16) k_UV,k_IR,t_UV,t_IR
  external derivs,rkqs
  real(16) eps_ode,h1,hmin !variables in subroutine odeint
  integer nok,nbad !variables in subroutine odeint
  INTEGER kmax,kount !variables in common block of subroutine odeint
  INTEGER KMAXX,NMAX
  PARAMETER (NMAX=50,KMAXX=2000)
  real(16) dxsav,xp(KMAXX),yp(NMAX,KMAXX) !variables in common block of subroutine odeint
  real(16) rho0,mPion,mSigma,mf,Vall
  real(16) fpi,h,Zphi,Zpsi,c,kappa
  integer i
  integer iT,iv
  real(16) l_com,lb_com


  common /strucFun/ N_str
  common /Tmu/ T,mu
  common /kRange/k_UV,k_IR,t_UV,t_IR
  common /odeContr/ eps_ode,h1,hmin
  COMMON /path/ kmax,kount,dxsav,xp,yp
  common /polyakov_com/ l_com,lb_com
  common /iTiv/ iT,iv


  N_str(1)=Nv
  N_str(2)=Nh
  N_str(3)=Nz
  N_str(4)=Nck

  k_UV=700.Q+0/hc !in unit of fm**(-1)
  k_IR=0.01Q+0/hc   !in unit of fm**(-1)
  t_UV=0.Q+0
  t_IR=log(k_IR/k_UV)

  eps_ode=1.Q-8
  h1=t_IR/200.Q+0
  hmin=0.Q+0
  kmax=KMAXX
  dxsav=t_IR/10000.Q+0
!  dxsav=t_IR/200.Q+0

  kappa_UV=kappa_UV_i
  call expaPoint(Nflow,yflow,kappa_UV)
  call initial(Nflow,yflow,kappa_UV)
  call odeint(yflow,Nflow,t_UV,t_IR,eps_ode,h1,hmin,nok,nbad,derivs,rkqs)
  call phypoint2(Nflow,yflow,rho0,mPion,mSigma,mf,Vall)

  fpi=sqrt(2.Q+0*rho0)
  h=yflow((Nv+1)+1)
  Zphi=yflow((Nv+1)+(Nh+1)+1)
  Zpsi=yflow((Nv+1)+(Nh+1)+2)
  c=yflow((Nv+1)+(Nh+1)+Nz+1)
  kappa=yflow((Nv+1)+(Nh+1)+Nz+2)

end


