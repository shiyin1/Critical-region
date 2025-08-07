subroutine Fnb(x,T,nb)          
!boson distribution function
  implicit none

  integer Nnum
  parameter(Nnum=400)
  real(16) x(Nnum),T
  real(16) nb(Nnum)
  real(16) over(Nnum)
  integer i

  over=x/T
  do i=1,Nnum
    if (over(i) > 8.Q+02)then
      nb(i)=0.Q+00
    else
      nb(i)=1.Q+00/(exp(over(i))-1.Q+00)
    end if
  end do
end

subroutine Fnf(x,T,nf)
!fermion distribution function
  implicit none

  integer Nnum
  parameter(Nnum=400)
  real(16) x(Nnum),T
  real(16) nf(Nnum)
  real(16) over(Nnum)
  integer i

  over=x/T
  do i=1,Nnum
    if(over(i) > 8.Q+02)then
      nf(i)=0.Q+00
    else
      nf(i)=1.Q+00/(exp(over(i))+1.Q+00)
    end if
  end do
end
