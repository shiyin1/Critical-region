subroutine derivs(x,y,dydx)
!Calculating the right hand side of differential equations under beyondLPA
!Calculating the Chebyshev coefficient
!Yong-Rui Chen change 2020.1.12

  implicit none

  real(16) pi
  parameter(pi=3.1415926535897932Q+00)
  real(16) Nc,Nf
  parameter(Nc=3.Q+00,Nf=2.Q+00)
  integer NMAX,Vmax
  !maximal number of differential equations
  !maximal number of potential
  parameter (NMAX=800,Vmax=2000)
  real(16) x,y(NMAX),dydx(NMAX),a(Vmax),b(Vmax),c(Vmax),dadk(Vmax),Tn(Vmax,Vmax),summation
  real(16) summationbn,summationb0,summationV0
  integer(8) gridnum,maxim,anum
  !number of rho
  real(16) rho(Vmax),rhomin,rhomax
  real(16) T,mu
  real(16) k 
  ! IR cutoff in flow equations
  real(16) h,cbar
  !Yukawa coupling and explicit chiral symmetry breaking 
  real(16) dVdrho0(Vmax)
  ! effective potential
  real(16) dVdrho1(Vmax),dVdrho2(Vmax),dVdrho3(Vmax),dVdrho4(Vmax),dVdrho5(Vmax)
  real(16) mpi2(Vmax),mpi2drho1(Vmax),mpi2drho2(Vmax),mpi2drho3(Vmax)
  real(16) msg2(Vmax),msg2drho1(Vmax),msg2drho2(Vmax),msg2drho3(Vmax)
  real(16) mqk2(Vmax),mqk2drho1(Vmax)
  ! sqr of mass 
  ! pi, sigma, quark
  real(16) Epi(Vmax),Esg(Vmax),Eqk(Vmax)
  real(16) nbpi(Vmax),nbsg(Vmax),nfqk(Vmax),naqk(Vmax)
  real(16) dnbpidx(Vmax,7),dnbsgdx(Vmax,7)
  real(16) nbpidx0(Vmax),nbpidx1(Vmax),nbpidx2(Vmax),nbpidx3(Vmax)
  real(16) nbsgdx0(Vmax),nbsgdx1(Vmax),nbsgdx2(Vmax),nbsgdx3(Vmax)
  real(16) dnfdx(Vmax,6),dnadx(Vmax,6)
  real(16) nfqkdx0(Vmax),nfqkdx1(Vmax),nfqkdx2(Vmax),nfqkdx3(Vmax)
  real(16) naqkdx0(Vmax),naqkdx1(Vmax),naqkdx2(Vmax),naqkdx3(Vmax)
  ! pi, sigma, quark
  real(16) etaphi(Vmax),etaphidrho1(Vmax),etaphidrho2(Vmax),etaphidrho3(Vmax)
  real(16) etaq(Vmax)
  real(16) Fall(Vmax,5)
  real(16) BBall(Vmax,10)
  real(16) B2B2(Vmax),B3B2(Vmax),B4B2(Vmax),B5B2(Vmax)
  real(16) B2B3(Vmax),B3B3(Vmax),B4B3(Vmax)
  real(16) B2B4(Vmax),B3B4(Vmax),B2B5(Vmax)
  real(16) F2(Vmax),F3(Vmax),F4(Vmax),F5(Vmax),F6(Vmax)
  real(16) lb0pi(Vmax),lb1pi(Vmax),lb2pi(Vmax),lb3pi(Vmax)
  real(16) lb0sg(Vmax),lb1sg(Vmax),lb2sg(Vmax),lb3sg(Vmax)
  real(16) lf0qk(Vmax),lf1qk(Vmax),lf2qk(Vmax),lf3qk(Vmax)
  real(16) k34pi2
  real(16) flowbn,flowb0,flowV0
  !calculate accelerate, no sence
  real(16) dVdkdrho0(Vmax),dVdkdrho1(Vmax),dVdkdrho2(Vmax),dVdkdrho3(Vmax)
  integer i,j,n 
  real(16) funsave(51)

  common /Tmu/ T,mu
  common /gridrho/ gridnum,anum,rho,rhomin,rhomax,maxim
  common /iniconst/ h,cbar
  common /chebyshev/ Tn,b,c
  common /anomalous/ etaphi
  common /sav/ funsave

  k=x
  a(1:anum)=y(1:anum)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call chder(rhomin,rhomax,a,b,anum)
  call chder(rhomin,rhomax,b,c,anum)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  open(unit=50,file='./buffer/a.dat')
!  do i=1,gridnum
!  write(50,*) a(i)
!  end do
!  close(50)
!  open(unit=50,file='./buffer/b.dat')
!  do i=1,gridnum
!  write(50,*) b(i)
!  end do
!  close(50)
!  open(unit=50,file='./buffer/c.dat')
!  do i=1,gridnum
!  write(50,*) c(i)
!  end do
!  close(50)

  do i=1,gridnum
    summation=0.5Q+00*b(1)
    do j=2,anum
    summation=summation+b(j)*Tn(j,i)
!    write(*,*)'summation=', summation
    end do
    dVdrho1(i)=summation
!    write(*,*)'dVdrho1=',dVdrho1(i)
  end do

  do i=1,gridnum
    summation=0.5Q+00*c(1)
    do j=2,anum
    summation=summation+c(j)*Tn(j,i)
    end do
    dVdrho2(i)=summation
!    write(*,*)dVdrho2(i)
  end do
!mass of mesons and quarks
  mpi2=dVdrho1

  msg2=dVdrho1+2.Q+00*rho*dVdrho2

  mqk2=h**2.Q+00*rho/2.Q+00
  mqk2drho1=h**2.Q+00/2.Q+00

  Epi=Sqrt(k**2.Q+00 + mpi2)
  Esg=Sqrt(k**2.Q+00 + msg2)
  Eqk=Sqrt(k**2.Q+00 + mqk2)

  call Fnb(Epi,T,nbpi)
  call nbdx(T,nbpi,dnbpidx)
  nbpidx0=dnbpidx(:,1)

  call Fnb(Esg,T,nbsg)
  call nbdx(T,nbsg,dnbsgdx)
  nbsgdx0=dnbsgdx(:,1)

  call Fnf(Eqk-mu,T,nfqk)
  call nfdx(T,nfqk,dnfdx)
  nfqkdx0=dnfdx(:,1)

  call Fnf(Eqk+mu,T,naqk)
  call nfdx(T,naqk,dnadx)
  naqkdx0=dnadx(:,1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  etaq=0.Q+00
  etaphi=0.Q+00!((-4.Q+00*(-2.Q+00 + etaq)*F3 + (-3.Q+00 + 2.Q+00*etaq)*F2)*h**2.Q+00*Nc +              &
      !(4.Q+00*B2B2*dVdrho2**2.Q+00*rho)/k**2.Q+00)/(6.Q+00*Pi**2.Q+00)
!  write(*,*)'etaphi=',etaphi(1:gridnum)
  etaphidrho1=0.Q+00
  etaphidrho2=0.Q+00
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  lb0pi=(2.Q+00*(1.Q+00 - etaphi/5.Q+00)*k*(0.5Q+00 + nbpidx0))/(3.Q+00*Sqrt(k**2.Q+00 + mpi2))
!  lb1pi=(-2.Q+00*etaphidrho1*k*(0.5Q+00 + nbpidx0))/(15.Q+00*Sqrt(k**2.Q+00 + mpi2)) -                &
!    ((1.Q+00 - etaphi/5.Q+00)*k*mpi2drho1*(0.5Q+00 + nbpidx0))/(3.Q+00*(k**2.Q+00 + mpi2)**1.5Q+00) + &
!    ((1.Q+00 - etaphi/5.Q+00)*k*mpi2drho1*nbpidx1)/(3.Q+00*(k**2.Q+00 + mpi2))
!  lb2pi=(-2.Q+00*etaphidrho2*k*(0.5Q+00 + nbpidx0))/(15.Q+00*Sqrt(k**2.Q+00 + mpi2)) +                &
!    (2.Q+00*etaphidrho1*k*mpi2drho1*(0.5Q+00 + nbpidx0))/(15.Q+00*(k**2.Q+00 + mpi2)**1.5Q+00) +      &
!    ((1.Q+00 - etaphi/5.Q+00)*k*mpi2drho1**2.Q+00*(0.5Q+00 + nbpidx0))/                               &
!     (2.Q+00*(k**2.Q+00 + mpi2)**2.5Q+00) -                                                           &
!    ((1.Q+00 - etaphi/5.Q+00)*k*mpi2drho2*(0.5Q+00 + nbpidx0))/(3.Q+00*(k**2.Q+00 + mpi2)**1.5Q+00) - &
!    (2.Q+00*etaphidrho1*k*mpi2drho1*nbpidx1)/(15.Q+00*(k**2.Q+00 + mpi2)) -                           &
!    ((1.Q+00 - etaphi/5.Q+00)*k*mpi2drho1**2.Q+00*nbpidx1)/(2.Q+00*(k**2.Q+00 + mpi2)**2.Q+00) +      &
!    ((1.Q+00 - etaphi/5.Q+00)*k*mpi2drho2*nbpidx1)/(3.Q+00*(k**2.Q+00 + mpi2)) +                      &
!    ((1.Q+00 - etaphi/5.Q+00)*k*mpi2drho1**2.Q+00*nbpidx2)/(6.Q+00*(k**2.Q+00 + mpi2)**1.5Q+00)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  lb0sg=(2.Q+00*(1.Q+00 - etaphi/5.Q+00)*k*(0.5Q+00 + nbsgdx0))/(3.Q+00*Sqrt(k**2.Q+00 + msg2))
!  lb1sg=(-2.Q+00*etaphidrho1*k*(0.5Q+00 + nbsgdx0))/(15.Q+00*Sqrt(k**2.Q+00 + msg2)) -                &
!    ((1.Q+00 - etaphi/5.Q+00)*k*msg2drho1*(0.5Q+00 + nbsgdx0))/(3.Q+00*(k**2.Q+00 + msg2)**1.5Q+00) + &
!    ((1.Q+00 - etaphi/5.Q+00)*k*msg2drho1*nbsgdx1)/(3.Q+00*(k**2.Q+00 + msg2))
!  lb2sg=(-2.Q+00*etaphidrho2*k*(0.5Q+00 + nbsgdx0))/(15.Q+00*Sqrt(k**2.Q+00 + msg2)) +                &
!    (2.Q+00*etaphidrho1*k*msg2drho1*(0.5Q+00 + nbsgdx0))/(15.Q+00*(k**2.Q+00 + msg2)**1.5Q+00) +      &
!    ((1.Q+00 - etaphi/5.Q+00)*k*msg2drho1**2.Q+00*(0.5Q+00 + nbsgdx0))/                               &
!     (2.Q+00*(k**2.Q+00 + msg2)**2.5Q+00) -                                                           &
!    ((1.Q+00 - etaphi/5.Q+00)*k*msg2drho2*(0.5Q+00 + nbsgdx0))/(3.Q+00*(k**2.Q+00 + msg2)**1.5Q+00) - &
!    (2.Q+00*etaphidrho1*k*msg2drho1*nbsgdx1)/(15.Q+00*(k**2.Q+00 + msg2)) -                           &
!    ((1.Q+00 - etaphi/5.Q+00)*k*msg2drho1**2.Q+00*nbsgdx1)/(2.Q+00*(k**2.Q+00 + msg2)**2.Q+00) +      &
!    ((1.Q+00 - etaphi/5.Q+00)*k*msg2drho2*nbsgdx1)/(3.Q+00*(k**2.Q+00 + msg2)) +                      &
!    ((1.Q+00 - etaphi/5.Q+00)*k*msg2drho1**2.Q+00*nbsgdx2)/(6.Q+00*(k**2.Q+00 + msg2)**1.5Q+00)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  lf0qk=((1.Q+00 - etaq/4.Q+00)*k*(1.Q+00 - naqk - nfqk))/(3.Q+00*Sqrt(k**2.Q+00 + mqk2))
!  lf1qk=-((1.Q+00 - etaq/4.Q+00)*k*mqk2drho1*(1.Q+00 - naqk - nfqk))/                                 &
!     (6.Q+00*(k**2.Q+00 + mqk2)**1.5Q+00) +                                                           &
!    ((1.Q+00 - etaq/4.Q+00)*k*(-(mqk2drho1*naqkdx1)/(2.Q+00*Sqrt(k**2.Q+00 + mqk2)) -                 &
!         (mqk2drho1*nfqkdx1)/(2.Q+00*Sqrt(k**2.Q+00 + mqk2))))/(3.Q+00*Sqrt(k**2.Q+00 + mqk2))
!  lf2qk=((1.Q+00 - etaq/4.Q+00)*k*mqk2drho1**2.Q+00*(1.Q+00 - naqk - nfqk))/                          &
!     (4.Q+00*(k**2.Q+00 + mqk2)**2.5Q+00) -                                                           &
!    ((1.Q+00 - etaq/4.Q+00)*k*mqk2drho1*                                                              &
!       (-(mqk2drho1*naqkdx1)/(2.Q+00*Sqrt(k**2.Q+00 + mqk2)) -                                        &
!         (mqk2drho1*nfqkdx1)/(2.Q+00*Sqrt(k**2.Q+00 + mqk2))))/                                       &
!     (3.Q+00*(k**2.Q+00 + mqk2)**1.5Q+00) +                                                           &
!    ((1.Q+00 - etaq/4.Q+00)*k*((mqk2drho1**2.Q+00*naqkdx1)/(4.Q+00*(k**2.Q+00 + mqk2)**1.5Q+00) -     &
!         (mqk2drho1**2.Q+00*naqkdx2)/(4.Q+00*(k**2.Q+00 + mqk2)) +                                    &
!         (mqk2drho1**2.Q+00*nfqkdx1)/(4.Q+00*(k**2.Q+00 + mqk2)**1.5Q+00) -                           &
!         (mqk2drho1**2.Q+00*nfqkdx2)/(4.Q+00*(k**2.Q+00 + mqk2))))/(3.Q+00*Sqrt(k**2.Q+00 + mqk2))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  k34pi2=k**3.Q+00/4.Q+00/pi/pi
  dVdkdrho0=k34pi2*(3.Q+00*lb0pi+lb0sg-4.Q+00*Nf*Nc*lf0qk)

  do i=1,anum
     summation=0.Q+00
     do j=1,gridnum
        summation=summation+dVdkdrho0(j)*Tn(i,j)
     end do
     dadk(i)=2.Q+00/gridnum*summation
  end do

  dydx=0.Q+00
  dydx(1:anum)=dadk(1:anum)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  funsave=b(1:51)
  write(*,*)'k=',k

end
