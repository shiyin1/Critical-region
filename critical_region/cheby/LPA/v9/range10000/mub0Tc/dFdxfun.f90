subroutine dfdxfun(u1,u2,dx,uxnum,u3,u4)
!To calculate d^3 u/(dx)^3 and d^4 u/(dx)^4 from d^u/dx and d^2 u/(dx)^2
!Rui Wen and Yong-rui Chen 20190413
      
  implicit none

  integer(8) uxnum
  integer(8) i,NMAX
  parameter(NMAX=3000)
  real(16) dx,u1(NMAX),u2(NMAX),u3(NMAX),u4(NMAX)
  real(16) a1(2*NMAX),a2(2*NMAX),a3(2*NMAX),a4(2*NMAX),a5(2*NMAX)
  real(16) b1(2*NMAX)
  real(16) Tem_a,u3u4(2*NMAX)

  b1=0.Q+00
  a1=0.Q+00
  a2=0.Q+00
  a3=0.Q+00
  a4=0.Q+00
  a5=0.Q+00
  u3=0.Q+00
  u4=0.Q+00
  u3u4=0.Q+00

  a1(1)=0.Q+00
  do i=1,uxnum-1
    a1(2*i)  =dx/2.Q+00
    a1(2*i+1)=-dx**3.Q+00/96.Q+00
  end do
  a1(2*uxnum)=1.Q+00

  a2(1)=1.Q+00
  do i=1,uxnum-1
    a2(2*i)  =dx**2.Q+00/8.Q+00
    a2(2*i+1)=-dx**2.Q+00/4.Q+00
  end do
  a2(2*uxnum)=-dx/4.Q+00

  a3(1)=dx/2.Q+00
  do i=1,uxnum-1
    a3(2*i)  =dx/2.Q+00
    a3(2*i+1)=5.Q+00*dx**3.Q+00/96.Q+00
  end do
  a3(2*uxnum)=0.Q+00

  a4(1)=-1.Q+00
  do i=1,uxnum-1
    a4(2*i)  =-dx**2.Q+00/8.Q+00
    a4(2*i+1)=0.Q+00
  end do
  a4(2*uxnum)=0.Q+00

  a5(1)=dx/2.Q+00
  a5(2:2*uxnum)=0.Q+00

  b1(1)=0.Q+00
  do i=1,uxnum-1
    b1(2*i)  =u2(i+1)-u2(i)
    b1(2*i+1)=u1(i+1)-u1(i)-dx*(3.Q+00*u2(i+1)+u2(i))/4.Q+00
  end do
  b1(2*uxnum)=(-6.Q+00*(u1(uxnum)-u1(uxnum-1))+dx*(5.Q+00*u2(uxnum)+u2(uxnum-1)))/2.Q+00/dx**2.Q+00

  do i=2,2*uxnum
    Tem_a= - a1(i)/a2(i-1)
    a2(i)=a2(i)+Tem_a*a3(i-1)
    a3(i)=a3(i)+Tem_a*a4(i-1)
    a4(i)=a4(i)+Tem_a*a5(i-1)
    b1(i)=b1(i)+Tem_a*b1(i-1)
  end do
  
  u3u4=0.Q+00
  u3u4(2*uxnum)=b1(2*uxnum)/a2(2*uxnum)
  u3u4(2*uxnum-1)=(b1(2*uxnum-1)-a3(2*uxnum-1)*u3u4(2*uxnum))/a2(2*uxnum-1)
  do i=2*uxnum-2,2,-1
    u3u4(i)=(b1(i)-a3(i)*u3u4(i+1)-a4(i)*u3u4(i+2))/a2(i)
  end do
  u3u4(1)=(b1(1)-a3(1)*u3u4(2)-a4(1)*u3u4(3)-a5(1)*u3u4(4))/a2(1)

  do i=1,uxnum
   u3(i)=u3u4(2*i-1)
   u4(i)=u3u4(2*i)
  end do

end
