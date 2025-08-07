SUBROUTINE splint(a,b,c,d,w,x,m,y)
  implicit none

  INTEGER m,n,p,Nnum,i
  parameter(Nnum=600)
  real(16) a(Nnum),b(Nnum),c(Nnum),d(Nnum),lampi(Nnum)
  real(16) w,y,x(Nnum)


   n=m-1
   p=1
   do i=1,n
     if(w<=x(i+1))then
             exit
     else
             p=p+1
     endif
   end do
   y=a(p)+b(p)*(w-x(p))+c(p)*(w-x(p))**2.Q+00+d(p)*(w-x(p))**3.Q+00
   end





