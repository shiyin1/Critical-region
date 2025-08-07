subroutine dfdxO5(fx,dx,fxnum,dfdx1,dfdx2)
!To calculate df/dx and d^2f/dx^2 up to O(5)
!Rui Wen and Yong-rui Chen 20190413
      
  implicit none

  integer(8) fxnum
  integer(8) i,NMAX
  parameter(NMAX=600)
  real(16) dx,fxU(NMAX),fx(NMAX),dfdx1(NMAX),dfdx2(NMAX)
  real(16) a1(3),b1(3),c1(3),d1(3),e1(3)
  real(16) a2(3),b2(3),c2(3),d2(3),e2(3)

  a1(1)=-25.Q+00/12.Q+00
  b1(1)=4.Q+00
  c1(1)=-3.Q+00
  d1(1)=4.Q+00/3.Q+00
  e1(1)=-1.Q+00/4.Q+00
 
  a2(1)=35.Q+00/12.Q+00
  b2(1)=-26.Q+00/3.Q+00
  c2(1)=19.Q+00/2.Q+00
  d2(1)=-14.Q+00/3.Q+00
  e2(1)=11.Q+00/12.Q+00

  a1(2)=-1.Q+00/4.Q+00
  b1(2)=-5.Q+00/6.Q+00
  c1(2)=3.Q+00/2.Q+00
  d1(2)=-1.Q+00/2.Q+00
  e1(2)=1.Q+00/12.Q+00

  a2(2)=11.Q+00/12.Q+00
  b2(2)=-5.Q+00/3.Q+00
  c2(2)=1.Q+00/2.Q+00
  d2(2)=1.Q+00/3.Q+00
  e2(2)=-1.Q+00/12.Q+00

  a1(3)=1.Q+00/12.Q+00
  b1(3)=-2.Q+00/3.Q+00
  c1(3)=0.Q+00
  d1(3)=2.Q+00/3.Q+00
  e1(3)=-1.Q+00/12.Q+00

  a2(3)=-1.Q+00/12.Q+00
  b2(3)=4.Q+00/3.Q+00
  c2(3)=-5.Q+00/2.Q+00
  d2(3)=4.Q+00/3.Q+00
  e2(3)=-1.Q+00/12.Q+00

  dfdx1(1)=(a1(1)*fx(1)+b1(1)*fx(2)+c1(1)*fx(3)+d1(1)*fx(4)+e1(1)*fx(5))/dx
  dfdx2(1)=(a2(1)*fx(1)+b2(1)*fx(2)+c2(1)*fx(3)+d2(1)*fx(4)+e2(1)*fx(5))/dx**2
  dfdx1(2)=(a1(2)*fx(1)+b1(2)*fx(2)+c1(2)*fx(3)+d1(2)*fx(4)+e1(2)*fx(5))/dx
  dfdx2(2)=(a2(2)*fx(1)+b2(2)*fx(2)+c2(2)*fx(3)+d2(2)*fx(4)+e2(2)*fx(5))/dx**2
  do i=3,fxnum-2
    dfdx1(i)=(a1(3)*fx(i-2)+b1(3)*fx(i-1)+c1(3)*fx(i)+d1(3)*fx(i+1)+e1(3)*fx(i+2))/dx
    dfdx2(i)=(a2(3)*fx(i-2)+b2(3)*fx(i-1)+c2(3)*fx(i)+d2(3)*fx(i+1)+e2(3)*fx(i+2))/dx**2
  end do
  dfdx1(fxnum-1)=-(a1(2)*fx(fxnum)+b1(2)*fx(fxnum-1)+c1(2)*fx(fxnum-2)+d1(2)*fx(fxnum-3)+e1(2)*fx(fxnum-4))/dx
  dfdx2(fxnum-1)=(a2(2)*fx(fxnum)+b2(2)*fx(fxnum-1)+c2(2)*fx(fxnum-2)+d2(2)*fx(fxnum-3)+e2(2)*fx(fxnum-4))/dx**2
  dfdx1(fxnum)=-(a1(1)*fx(fxnum)+b1(1)*fx(fxnum-1)+c1(1)*fx(fxnum-2)+d1(1)*fx(fxnum-3)+e1(1)*fx(fxnum-4))/dx
  dfdx2(fxnum)=(a2(1)*fx(fxnum)+b2(1)*fx(fxnum-1)+c2(1)*fx(fxnum-2)+d2(1)*fx(fxnum-3)+e2(1)*fx(fxnum-4))/dx**2

end
