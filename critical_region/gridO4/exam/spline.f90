SUBROUTINE spline(x,y,m,a,b,c,d)
  implicit none

  integer m,n,Nnum,i
  parameter(Nnum=600)
  real(16) x(Nnum),y(Nnum)
  real(16) a(Nnum),b(Nnum),c(Nnum),d(Nnum)
  real(16) h(Nnum),u(Nnum)
  real(16) s(Nnum),z(Nnum),t(Nnum)

  n=m-1
  do i=1,m
    a(i)=y(i)
  end do

  do i=1,n
    h(i)=x(i+1)-x(i)
  end do
  u(1)=0.Q+00
  do i=2,n
    u(i)=3.Q+00*(a(i+1)-a(i))/h(i)-3.Q+00*(a(i)-a(i-1))/h(i-1)
  end do
  s(1)=1.Q+00
  z(1)=0.Q+00
  t(1)=0.Q+00
  do i=2,n
    s(i)=2.Q+00*(x(i+1)-x(i-1))-h(i-1)*t(i-1)
    t(i)=h(i)/s(i)
    z(i)=(u(i)-h(i-1)*z(i-1))/s(i)
  end do
  s(m)=1.Q+00
  z(m)=0.Q+00
  c(m)=0.Q+00
  do i=n,1,-1
    c(i)=z(i)-t(i)*c(i+1)
    b(i)=(a(i+1)-a(i))/h(i)-h(i)*(c(i+1)+2*c(i))/3.Q+00
    d(i)=(c(i+1)-c(i))/(3.Q+00*h(i))
  end do
  end


