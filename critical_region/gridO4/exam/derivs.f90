subroutine derivs(x,y,dydx)
!Calculating the right hand side of differential equations
!Calculating the Chebyshev coefficient

  implicit none

  real(16) pi
  parameter(pi=3.14159265358979323846264338327950288419716939937510Q+00)
  integer NMAX,Vmax,maxim
  !maximal number of differential equations
  !maximal number of potential
  parameter (NMAX=600,Vmax=600)
  real(16) x,y(NMAX),dydx(NMAX)
  real(16) k_UV
  parameter(k_UV=1000.Q+00)
  real(16) summation,summation1,summation2,summation3
  integer gridnum
  !number of rho
  real(16) drho,rho(Vmax),rhomin,rhomax
  real(16) rho0,rhoacc,rtbis
  real(16) T,Num,d,const
  real(16) k 
  ! IR cutoff in flow equations
  real(16) h,cbar
  !Yukawa coupling and explicit chiral symmetry breaking 
  ! effective potential
  real(16) dVdrho1(Vmax),dVdrho2(Vmax),dVdrho3(Vmax)
  real(16) dVdrho1min,dVdrho2min,dVdrho3min
  real(16) mpi2(Vmax),msg2(Vmax)
  real(16) mpi2min,msg2min
  real(16) lb1pi(Vmax),lb1sg(Vmax)
  ! sqr of mass 
  ! pi, sigma, quark
  real(16) aV1(Vmax),bV1(Vmax),cV1(Vmax),dV1(Vmax)
  real(16) aV2(Vmax),bV2(Vmax),cV2(Vmax),dV2(Vmax)
  real(16) aV3(Vmax),bV3(Vmax),cV3(Vmax),dV3(Vmax)
  real(16) etaphi,Zpi
  real(16) etatpi,etatsigma,Ztpi,Ztsigma
  !calculate accelerate, no sence
  real(16) dVdtdrho0(Vmax),dVdtdrho1(Vmax),dVdtdrho2(Vmax)
  real(16) dZpidt,dZsigmadt
  real(16) dZtpidt,dZtsigmadt
  real(16) funsave(10)
  integer i,j,nn
  parameter(nn=1) 
  logical check
  external gapEq

  common /Tmu/ T,Num,d,const
  common /gridrho/ gridnum,drho,rho,rhomin,rhomax,rhoacc
  common /iniconst/ h,cbar
  common /interpol/ aV1,bV1,cV1,dV1
  common /odecomm/ funsave

  k=exp(x)*k_UV
  dVdrho1(1:gridnum)=y(1:gridnum)
  Zpi=y(gridnum+1)
  Ztpi=y(gridnum+2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call dfdxO5(dVdrho1,drho,gridnum,dVdrho2,dVdrho3)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call spline(rho,dVdrho1,gridnum,aV1,bV1,cV1,dV1)
  if(dVdrho1(1)>=0.Q+00)then
          rho0=rho(1)
  elseif(dVdrho1(gridnum)<=0.Q+00)then
          rho0=rho(gridnum)
  else
          rho0=rtbis(gapEq,rhomin,rhomax,rhoacc)
  endif

  call spline(rho,dVdrho2,gridnum,aV2,bV2,cV2,dV2)
  call spline(rho,dVdrho3,gridnum,aV3,bV3,cV3,dV3)

  call splint(aV1,bV1,cV1,dV1,rho0,rho,gridnum,dVdrho1min)
  call splint(aV2,bV2,cV2,dV2,rho0,rho,gridnum,dVdrho2min)
  call splint(aV3,bV3,cV3,dV3,rho0,rho,gridnum,dVdrho3min)
  if(abs(dVdrho1min)<1Q-9)then
          dVdrho1min=0.Q+00
  endif
  
  mpi2min=dVdrho1min
  msg2min=dVdrho1min+2.Q+00*rho0*dVdrho2min


!mass of mesons and quarks
  mpi2=dVdrho1
  msg2=dVdrho1+2.Q+00*rho*dVdrho2



  etaphi=4.Q+00*const*rho0*dVdrho2min**2.Q+00/((1.Q+00+mpi2min)**2.Q+00*(1.Q+00+msg2min)**2.Q+00)
 

  etatpi=(-2.Q+00*dVdrho2min**2.Q+00*(-5.Q+00 + etaphi)*&
      (6.Q+00 + mpi2min**2.Q+00 + 6.Q+00*msg2min + msg2min**2.Q+00 + mpi2min*(6.Q+00 + 4.Q+00*msg2min))*&
      rho0)/(15.Q+00*(1.Q+00 + mpi2min)**2.Q+00*(1.Q+00 + msg2min)**2.Q+00*&
      (2.Q+00 + mpi2min + msg2min)**2.Q+00*Pi**2.Q+00)


  lb1pi=dVdrho2/(1.Q+00+mpi2)**2.Q+00
  lb1sg=(3.Q+00*dVdrho2+2.Q+00*rho*dVdrho3)/(1.Q+00+msg2)**2.Q+00

  dVdtdrho1=(-2.Q+00+etaphi)*dVdrho1+(d-2.Q+00+etaphi)*rho*dVdrho2-const*(1.Q+00-etaphi/(d+2.Q+00))*(lb1sg+(Num-1)*lb1pi)


  dZpidt=-etaphi*Zpi
  dZtpidt=-etatpi*Ztpi



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  dydx=0.Q+00
  dydx(1:gridnum)=dVdtdrho1(1:gridnum)
  dydx(gridnum+1)=dZpidt
  dydx(gridnum+2)=dZtpidt


  funsave(1)=rho0
  funsave(2)=etaphi
  funsave(3)=etatpi
  funsave(4)=mpi2min
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  write(*,*)'k=',k
!  open(unit=50,file='./buffer/etapi.dat')
!  write(50,*)etapi
!  close(50)
!  open(unit=50,file='./buffer/etasigma.dat')
!  write(50,*)etasigma
!  close(50)
!  if(k<200)then
!  stop
end
