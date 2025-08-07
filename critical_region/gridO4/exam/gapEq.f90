subroutine gapEq(x,fvec) 

  implicit none

  integer n,Vmax,gridnum
  parameter(Vmax=600)
  real(16) x, fvec ,aV1(Vmax),bV1(Vmax),cV1(Vmax),dV1(Vmax)                    !x(1)
  real(16) drho,rho(Vmax),rhomin,rhomax,rhoacc
  real(16) Num,d,const
  real(16) h,cbar
  real(16) rho0,V0

  common /gridrho/ gridnum,drho,rho,rhomin,rhomax,rhoacc
  common /interpol/ aV1,bV1,cV1,dV1 

  

!  write(*,*)gridnum
!  stop
  call splint(aV1,bV1,cV1,dV1,x,rho,gridnum,fvec)


 !open(unit=50,file='./buffer/rho0.dat',position='append')
 !write(*,*)rho0
 !close(50)
 ! open(unit=50,file='./buffer/rho0.dat',position='append')
 ! write(*,*)V0
 ! close(50)

end 
