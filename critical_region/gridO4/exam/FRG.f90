Program FRG
!this code solve QM model by chebyshev method
!Yong-Rui Chen 2020.1.6

  implicit none
  
  real(16) pi
  parameter (pi=3.14159265358979323846264338327950288419716939937510Q+00)
  real(16) T,Num,d,const
  real(16) rho0up,rho0down
  integer Tnum,munum
  real(16) k_UV,k_IR,t_UV,t_IR
  integer Nflow,Vmax
  !number of flow equations
  parameter (Nflow=600,Vmax=600)
  real(16) aV1(Vmax),bV1(Vmax),cV1(Vmax),dV1(Vmax)
  real(16) aV2(Vmax),bV2(Vmax),cV2(Vmax),dV2(Vmax)
  real(16) aV3(Vmax),bV3(Vmax),cV3(Vmax),dV3(Vmax)
  real(16) dVdrho1(Vmax),dVdrho2(Vmax),dVdrho3(Vmax)
  !Chebyshev coefficient
  integer gridnum,maxim
  !number of rho
  real(16) drho,rho(Vmax),rhomin,rhomax,yflow(Nflow),rhoacc
  external derivs,rkqs
  real(16) eps_ode,h1,hmin 
  !variables in subroutine odeint
  integer nok,nbad 
  !variables in subroutine odeint
  integer kmax,kount
  !variables in common block of subroutine odeint
  integer KMAXX,NMAX
  parameter (NMAX=600,KMAXX=1000)
  real(16) dxsav,xp(KMAXX),yp(NMAX,KMAXX),funp(10,KMAXX)
  !variables in common block of subroutine odeint
  real(16) h,cbar
  !Yukawa coupling and explicit chiral symmetry breaking 
  integer i,imin,j,k,j1
  real(16) dydx(Nflow),y(Nflow)
  real(16) rho0

  common /Tmu/ T,Num,d,const
  common /gridrho/ gridnum,drho,rho,rhomin,rhomax,rhoacc
  common /odeContr/ eps_ode,h1,hmin
  common /path/ kmax,kount,dxsav,xp,yp,funp
  common /iniconst/ h,cbar
  common /interpol/ aV1,bV1,cV1,dV1
  
 
!  T=real(Tnum,kind=16)
!  mu=0.1Q+00*real(munum,kind=8)+245.Q+00

  T=150.Q+00
  Num=4.Q+00
  d=3.Q+00
  const=1.Q+00/(6.Q+00*pi**2.Q+00) 

  k_UV=1000.Q+00
  k_IR=0.00001Q+00

  t_UV=0.Q+00
  t_IR=log(k_IR/k_UV)

  eps_ode=1.Q-012
  h1=(k_IR-k_UV)/1000.Q+00
  hmin=0.Q+00
  kmax=KMAXX
  dxsav=(k_IR-k_UV)/10000.Q+00

  gridnum=500
!  gridnum<Vmax=100
  
  rhomin=0.Q+00
  rhomax=15.Q+00
  rhoacc=1Q-10
  drho=rhomax/real(gridnum-1,kind=16)
  do i=1,gridnum
  rho(i)=real(i-1,kind=16)*drho
  end do

  rho0=0.0771353129090741525887632114887237474Q+00


  call initial(yflow,rho0)
  call odeint(yflow,gridnum+2,t_UV,t_IR,eps_ode,h1,hmin,nok,nbad,derivs,rkqs)
  call phyphi(yflow,rho0,imin)
!  call Int_Z(K_UV,K_IR,rho_min)

  write(*,*)'T=',T
!  write(*,*)'phi0=',sqrt(2.Q+00*rho(imin))
!  write(*,*)'mqk=',sqrt(h**2.Q+00*rho(imin)/2.Q+00)


  open(unit=50,file='./buffer/xp.dat')
  do i=1,kount
  write(50,*)exp(xp(i))*k_UV
  end do
  close(50)
  open(unit=50,file='./buffer/rho0_k.dat')
  do i=1,kount
  write(50,*)funp(1,i)
  end do
  close(50)
  open(unit=50,file='./buffer/etapi_k.dat')
  do i=1,kount
  write(50,*)funp(2,i)
  end do
  close(50)
  open(unit=50,file='./buffer/etatpi_k.dat')
  do i=1,kount
  write(50,*)funp(3,i)
  end do
  close(50)
  open(unit=50,file='./buffer/mpi2min_k.dat')
  do i=1,kount
  write(50,*)funp(4,i)
  end do
  close(50)
end
