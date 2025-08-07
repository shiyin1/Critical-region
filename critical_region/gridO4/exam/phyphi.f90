subroutine phyphi(yflow,rho0,imin)

  implicit none
  integer Nflow,Vmax,maxim
  parameter (Nflow=600,Vmax=600)
  integer gridnum,anum
  real(16) drho,rho(Vmax),yflow(Nflow),rhomin,rhomax
  real(16) rho0,rhoacc,rtbis
  real(16) V0(Vmax),V1(Vmax),V2(Vmax),V3(Vmax)
  real(16) summation1,summation2,summation3
  real(16) dVdrho1_0,dVdrho1_max
  real(16) dVdrho1min,dVdrho2min,dVdrho3min
  real(16) Zpi,Zsigma,Ztpi,Ztsigma
  integer i,imin,j
  real(16) h,cbar
  external gapEq

  common /gridrho/ gridnum,drho,rho,rhomin,rhomax,rhoacc
  common /iniconst/ h,cbar
 
  V1(1:gridnum)=yflow(1:gridnum)
  Zpi=yflow(gridnum+1)
  Ztpi=yflow(gridnum+2)

  call dfdxO5(V1,drho,gridnum,V2,V3) 

   if(V1(1)>0.Q+00)then
    imin=0
   elseif(V1(1)<0.Q+00)then
    imin=1
   end if 



!  V0=V0-cbar*sqrt(2.Q+00*rhonew)

  write(*,*)'V1=',V1(1)
  write(*,*)'V2=',V2(1)
  open(unit=50,file='./buffer/rho.dat')
  do i=1,gridnum
  write(50,*)rho(i)
  end do
  close(50)
  open(unit=50,file='./buffer/V1.dat')
  do i=1,gridnum
  write(50,*)V1(i)
  end do
  close(50)
  open(unit=50,file='./buffer/V1_bis.dat',position='append')
  write(50,*) V1(1)
  close(50)
  open(unit=50,file='./buffer/V2.dat',position='append')
  write(50,*) V2(1)
  close(50)
  open(unit=50,file='./buffer/Zpi.dat')
  write(50,*)Zpi
  close(50)
  open(unit=50,file='./buffer/Ztpi.dat')
  write(50,*)Ztpi
  close(50)

end

