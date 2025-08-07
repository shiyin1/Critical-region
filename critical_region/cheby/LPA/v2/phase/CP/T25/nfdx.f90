subroutine nbdx(T,nb,dnbdx)
!Calculate the derivatives of bose distribution functions
  implicit none

  integer Nnum
  parameter(Nnum=400)
  real(16) nb(Nnum),dnbdx(Nnum,7)
  real(16) nbd0x(Nnum),nbd1x(Nnum),nbd2x(Nnum),nbd3x(Nnum),nbd4x(Nnum),nbd5x(Nnum),nbd6x(Nnum)
  real(16) T
  !temperature and chemical potential

  nbd0x=nb
  nbd1x=-((nb*(1.Q+00 + nb))/T)
  nbd2x=(nb*(1.Q+00 + nb)*(1.Q+00 + 2.Q+00*nb))/T**2.Q+00
  nbd3x=-((nb*(1.Q+00 + nb)*(1.Q+00 + 6.Q+00*nb*(1.Q+00 + nb)))/T**3.Q+00)
  nbd4x=(nb*(1.Q+00 + nb)*(1.Q+00 + 2.Q+00*nb)*(1.Q+00 + 12.Q+00*nb*(1.Q+00 + nb)))/T**4.Q+00
  nbd5x=-((nb*(1.Q+00 + nb)*(1.Q+00 + 30.Q+00*nb*(1.Q+00 + nb)*(1.Q+00 + 2.Q+00*nb)**2.Q+00))/T**5.Q+00)
  nbd6x=(nb*(1.Q+00 + nb)*(1.Q+00 + 2.Q+00*nb)*(1.Q+00 + 60.Q+00*nb*(1.Q+00 + nb)*(1.Q+00 + 6.Q+00*nb*(1.Q+00 + nb))))/T**6.Q+00

  dnbdx(:,1)=nbd0x
  dnbdx(:,2)=nbd1x
  dnbdx(:,3)=nbd2x
  dnbdx(:,4)=nbd3x
  dnbdx(:,5)=nbd4x
  dnbdx(:,6)=nbd5x
  dnbdx(:,7)=nbd6x

end

subroutine nfdx(T,nf,dnfdx)
!Calculate the derivatives of bose distribution functions
  implicit none

  integer Nnum
  parameter(Nnum=400)
  real(16) nf(Nnum)
  real(16) nfd0x(Nnum),nfd1x(Nnum),nfd2x(Nnum),nfd3x(Nnum),nfd4x(Nnum),nfd5x(Nnum)
  real(16) dnfdx(Nnum,6)
  real(16) T
  !temperature and chemical potential

  nfd0x=nf
  nfd1x=((-1.Q+00 + nf)*nf)/T
  nfd2x=((-1.Q+00 + nf)*nf*(-1.Q+00 + 2.Q+00*nf))/T**2.Q+00
  nfd3x=((-1.Q+00 + nf)*nf*(1.Q+00 - 6.Q+00*nf + 6.Q+00*nf**2.Q+00))/T**3.Q+00
  nfd4x=((-1.Q+00 + nf)*nf*(-1.Q+00 + 2.Q+00*nf)*(1.Q+00 + 12.Q+00*(-1.Q+00 + nf)*nf))/T**4.Q+00
  nfd5x=((-1.Q+00 + nf)*nf*(1.Q+00 + 30.Q+00*(1.Q+00 - 2.Q+00*nf)**2.Q+00*(-1.Q+00 + nf)*nf))/T**5.Q+00

  dnfdx(:,1)=nfd0x
  dnfdx(:,2)=nfd1x
  dnfdx(:,3)=nfd2x
  dnfdx(:,4)=nfd3x
  dnfdx(:,5)=nfd4x
  dnfdx(:,6)=nfd5x
end
