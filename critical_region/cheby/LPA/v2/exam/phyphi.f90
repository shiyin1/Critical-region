subroutine phyphi(yflow,V0,imin)

  implicit none
  integer Nflow,Vmax,maxim
  parameter (Nflow=400,Vmax=100)
  integer(8) gridnum,anum
  real(16) drho,rho(Vmax),yflow(Nflow),rhonew(100000),rhomin,rhomax,etaphi(Vmax)
  real(16) V0(100000),V1(100000),V2(100000),a(Vmax),Tn(Vmax,Vmax),b(Vmax),c(Vmax)
  real(16) summation0,summation1,summation2
  integer i,imin,j
  real(16) h,cbar,V2_sig(100000)

  common /gridrho/ gridnum,anum,rho,rhomin,rhomax,maxim
  common /iniconst/ h,cbar
  common /chebyshev/ Tn,b,c
  common /anomalous/ etaphi
  
  a(1:anum)=yflow(1:anum)

  
!  do i=1,gridnum
!    summation0=0.5Q+00*a(1)
!    summation1=0.5Q+00*b(1)
!    summation2=0.5Q+00*c(1)
!    do j=2,gridnum
!    summation0=summation0+a(j)*Tn(j,i)
!    summation1=summation1+b(j)*Tn(j,i)
!    summation2=summation2+c(j)*Tn(j,i)
!    end do
!    V0(i)=summation0
!    V1(i)=summation1
!    V2(i)=summation2
!  end do
  !V versus phi rather than rho

!  V0=V0-cbar*sqrt(2.Q+00*rho)

   maxim=20001
   drho=(rhomax-rhomin)/real(maxim-1,kind=16)
   do i=1,maxim
   rhonew(i)=real(i-1,kind=16)*drho
   end do

   do i=1,maxim
   summation0=0.5Q+00*a(1)
   summation1=0.5Q+00*b(1)
   summation2=0.5Q+00*c(1)
   do j=2,anum
   summation0=summation0+a(j)*cos((j-1)*Acos((2*rhonew(i)-rhomin-rhomax)/(rhomax-rhomin)))
   summation1=summation1+b(j)*cos((j-1)*Acos((2*rhonew(i)-rhomin-rhomax)/(rhomax-rhomin)))
   summation2=summation2+c(j)*cos((j-1)*Acos((2*rhonew(i)-rhomin-rhomax)/(rhomax-rhomin)))
   end do
   V0(i)=summation0
   V1(i)=summation1
   V2(i)=summation2
   end do

  !V0=V0-cbar*sqrt(2.Q+00*rhonew)

  imin=1
  do i=1,maxim
    if(V0(i)<V0(imin))then
      imin=i
    end if
  end do
  write(*,*)'phi0=',sqrt(2.Q+00*rhonew(imin))
  write(*,*)'mqk=',sqrt(h**2.Q+00*rhonew(imin)/2.Q+00)
  write(*,*)'mpi=',sqrt(V1(imin))
  write(*,*)'msg=',sqrt(V1(imin)+2.Q+00*rhonew(imin)*V2(imin))
  open(unit=50,file='./buffer/a.dat')
  do i=1,anum
  write(50,*) a(i)
  end do
  close(50)
  open(unit=50,file='./buffer/b.dat')
  do i=1,anum
  write(50,*) b(i)
  end do
  close(50)
  open(unit=50,file='./buffer/c.dat')
  do i=1,anum
  write(50,*) c(i)
  end do
  close(50)
   open(unit=50,file='./buffer/etaphi.dat')
  do i=1,gridnum
  write(50,*) etaphi(i)
  end do
  close(50)
  open(unit=50,file='./buffer/phi0.dat',position='append')
  write(50,*) sqrt(2.Q+00*rhonew(imin))
  close(50)

end

