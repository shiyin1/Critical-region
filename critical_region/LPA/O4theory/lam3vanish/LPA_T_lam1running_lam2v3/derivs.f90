subroutine derivs(x,y,dydx)
!Calculating the right hand side of differential equations

  implicit none

  integer NMAX !maximal number of differential equations
  parameter(NMAX=50)
  real(16) x,y(NMAX),dydx(NMAX)
  integer N_str(4) !store the structure of functions of ODE
  integer Nv,Nh,Nz,Nck
  real(16) k ! IR cutoff in flow equations
  real(16) lam0,lam1,lam2,lam3,lam4,lam5,lam6,lam7
  real(16) h
  real(16) Zphi,Zpsi
  real(16) c,kappa
  real(16) pi,hc
  parameter(pi=3.141592653589793238462643383279Q+0)
  parameter(hc=197.33Q+0)
  real(16) etaphi,etapsi
!meson and quark anomanous dimension
  real(16) Nc,Nf
  parameter(Nc=0.Q+0,Nf=2.Q+0)
  real(16) v3
  parameter(v3=1.Q+0/(2.Q+0*pi**2))
  real(16) rho !phi_a**2/2
  real(16) Fnb,Fnf0,Fnf1,Fnf2
  external Fnb,Fnf0,Fnf1,Fnf2
  real(16) zb,zf !distinguish the transverse and longituidanl wave function renormalization
  real(16) p0,p0c !temporal compontent of external momentum
  real(16) T,mu
  real(16) mu0
  real(16) l,lb !polyakov loop
  real(16) k_UV,k_IR,t_UV,t_IR
  real(16) mp2,ms2,mf2,mp2d1rho,mp2d2rho,mp2d3rho,mp2d4rho,mp2d5rho,ms2d1rho,ms2d2rho,ms2d3rho,ms2d4rho,ms2d5rho,mf2d1rho
  real(16) nb,nbd0x,nbd1x,nbd2x,nbd3x,nbd4x,nbd5x
  real(16) nbPion,nbd1xPion,nbd2xPion,nbd3xPion,nbd4xPion,nbd5xPion,nbSigma,nbd1xSigma,nbd2xSigma, &
          nbd3xSigma,nbd4xSigma,nbd5xSigma
  real(16) nff,nfd1xf,nfd2xf,nfd3xf,nfd4xf,nfd5xf,nfa,nfd1xa,nfd2xa,nfd3xa,nfd4xa,nfd5xa
  real(16) nf0,nf1,nf2,nfd0x,nfd1x,nfd2x,nfd3x,nfd4x,nfd5x
  real(16) xff,xfa
  real(16) f2a,f3a,b2b2PS,b2f1aPion,b2f1aPionI,b2f1aSigma,b2f1aSigmaI,b1f2Pion,b1f2PionI,b1f2Sigma,b1f2SigmaI
  complex(16) b2f1aPionC,b2f1aSigmaC,b1f2PionC,b1f2SigmaC
  real(16) dr0dtV,dr1dtV,dr2dtV,dr3dtV,dr4dtV,dr5dtV
  real(16) dlam0dt,dlam1dt,dlam2dt,dlam3dt,dlam4dt,dlam5dt,dlam6dt,dlam7dt
  real(16) L11Pion,L11Sigma
  real(16) dth,dhdt
  real(16) dZphidt,dZpsidt,dcdt,dkappadt
  real(16) l_com,lb_com


  common /strucFun/ N_str
  common /Tmu/ T,mu
  common /kRange/k_UV,k_IR,t_UV,t_IR
  common /polyakov_com/ l_com,lb_com


  Nv=N_str(1)
  Nh=N_str(2)
  Nz=N_str(3)
  Nck=N_str(4)

  k=k_UV*exp(x)
  lam1=y(1)
  lam2=y(2)
  lam3=y(3)
  lam4=y(4)
  lam5=y(5)
  lam6=y(6)
  lam7=y(7)
  lam0=y(Nv+1)
  h=y((Nv+1)+1)
  Zphi=y((Nv+1)+(Nh+1)+1)
  Zpsi=y((Nv+1)+(Nh+1)+2)
  c=y((Nv+1)+(Nh+1)+Nz+1)
  kappa=y((Nv+1)+(Nh+1)+Nz+2)


  rho=kappa !calculations are performed at expansion point kappa

  zb=1.Q+0
  zf=1.Q+0

  mu0=0.Q+0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!mass and their derivatives
  mp2=lam1/k**2

  ms2=(lam1 + 2*lam2*rho)/k**2

  mf2=(h**2*rho)/(k**2*Nf)

  mp2d1rho=lam2/k**2

  mp2d2rho=lam3/k**2

  mp2d3rho=lam4/k**2

  mp2d4rho=lam5/k**2

  mp2d5rho=lam6/k**2

  ms2d1rho=(3*lam2 + 2*lam3*rho)/k**2

  ms2d2rho=(5*lam3 + 2*lam4*rho)/k**2

  ms2d3rho=(7*lam4 + 2*lam5*rho)/k**2

  ms2d4rho=(9*lam5 + 2*lam6*rho)/k**2

  ms2d5rho=(11*lam6 + 2*lam7*rho)/k**2

  mf2d1rho=h**2/(k**2*Nf)




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  nb=Fnb((k*Sqrt(1 + mp2))/Sqrt(zb),T)
  call nbdx(nb,nbd0x,nbd1x,nbd2x,nbd3x,nbd4x,nbd5x)
  nbPion=nbd0x
  nbd1xPion=nbd1x
  nbd2xPion=nbd2x
  nbd3xPion=nbd3x
  nbd4xPion=nbd4x
  nbd5xPion=nbd5x

  nb=Fnb((k*Sqrt(1 + ms2))/Sqrt(zb),T)
  call nbdx(nb,nbd0x,nbd1x,nbd2x,nbd3x,nbd4x,nbd5x)
  nbSigma=nbd0x
  nbd1xSigma=nbd1x
  nbd2xSigma=nbd2x
  nbd3xSigma=nbd3x
  nbd4xSigma=nbd4x
  nbd5xSigma=nbd5x


  xff=-mu + (k*Sqrt(1 + mf2))/zf
  xfa=mu + (k*Sqrt(1 + mf2))/zf

  l=l_com
  lb=lb_com

  nf0=Fnf0(xff,T,l,lb)
  nf1=Fnf1(xff,T,l,lb)
  nf2=Fnf2(xff,T,l,lb)
  call nfdx(l,lb,nf0,nf1,nf2,nfd0x,nfd1x,nfd2x,nfd3x,nfd4x,nfd5x)
  nff=nfd0x
  nfd1xf=nfd1x
  nfd2xf=nfd2x
  nfd3xf=nfd3x
  nfd4xf=nfd4x
  nfd5xf=nfd5x



  nf0=Fnf0(xfa,T,lb,l)
  nf1=Fnf1(xfa,T,lb,l)
  nf2=Fnf2(xfa,T,lb,l)
  call nfdx(lb,l,nf0,nf1,nf2,nfd0x,nfd1x,nfd2x,nfd3x,nfd4x,nfd5x)
  nfa=nfd0x
  nfd1xa=nfd1x
  nfd2xa=nfd2x
  nfd3xa=nfd3x
  nfd4xa=nfd4x
  nfd5xa=nfd5x

  etaphi=0.
  etapsi=0.

  dZphidt=-etaphi*Zphi
  dZpsidt=-etapsi*Zpsi

  dcdt=(1.Q+0/2.Q+0)*etaphi*c
  dkappadt=-etaphi*kappa

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  dr0dtV=(k**4*v3*((2*(1 - etaphi/5.Q+0)*(0.5Q+0 + nbSigma))/(3.Q+0*Sqrt(1 + ms2)*Sqrt(zb)) +   &
      (2*(1 - etaphi/5.Q+0)*(0.5Q+0 + nbPion)*(-1 + Nf**2))/                          &
       (3.Q+0*Sqrt(1 + mp2)*Sqrt(zb))))/2.Q+0

  dr1dtV=(k**4*v3*(((1 - etaphi/5.Q+0)*k*ms2d1rho*nbd1xSigma)/(3.Q+0*(1 + ms2)*zb) +    &
      ((1 - etaphi/5.Q+0)*k*mp2d1rho*nbd1xPion*(-1 + Nf**2))/                      &
       (3.Q+0*(1 + mp2)*zb) - ((1 - etaphi/5.Q+0)*ms2d1rho*(0.5Q+0 + nbSigma))/          &
       (3.Q+0*(1 + ms2)**1.5Q+0*Sqrt(zb)) -                                           &
      ((1 - etaphi/5.Q+0)*mp2d1rho*(0.5Q+0 + nbPion)*(-1 + Nf**2))/                   &
       (3.Q+0*(1 + mp2)**1.5Q+0*Sqrt(zb))))/2.Q+0

  dr2dtV=(k**4*v3*(((1 - etaphi/5.Q+0)*k**2*ms2d1rho**2*nbd2xSigma)/               &
       (6.Q+0*(1 + ms2)**1.5Q+0*zb**1.5Q+0) +                                            &
      ((1 - etaphi/5.Q+0)*k**2*mp2d1rho**2*nbd2xPion*(-1 + Nf**2))/                &
       (6.Q+0*(1 + mp2)**1.5Q+0*zb**1.5Q+0) -                                            &
      ((1 - etaphi/5.Q+0)*k*ms2d1rho**2*nbd1xSigma)/(2.Q+0*(1 + ms2)**2*zb) +         &
      ((1 - etaphi/5.Q+0)*k*ms2d2rho*nbd1xSigma)/(3.Q+0*(1 + ms2)*zb) -               &
      ((1 - etaphi/5.Q+0)*k*mp2d1rho**2*nbd1xPion*(-1 + Nf**2))/                   &
       (2.Q+0*(1 + mp2)**2*zb) + ((1 - etaphi/5.Q+0)*k*mp2d2rho*nbd1xPion*            &
         (-1 + Nf**2))/(3.Q+0*(1 + mp2)*zb) +                                      &
      ((1 - etaphi/5.Q+0)*ms2d1rho**2*(0.5Q+0 + nbSigma))/                            &
       (2.Q+0*(1 + ms2)**2.5Q+0*Sqrt(zb)) -                                           &
      ((1 - etaphi/5.Q+0)*ms2d2rho*(0.5Q+0 + nbSigma))/                               &
       (3.Q+0*(1 + ms2)**1.5Q+0*Sqrt(zb)) +                                           &
      ((1 - etaphi/5.Q+0)*mp2d1rho**2*(0.5Q+0 + nbPion)*(-1 + Nf**2))/                &
       (2.Q+0*(1 + mp2)**2.5Q+0*Sqrt(zb)) -                                           &
      ((1 - etaphi/5.Q+0)*mp2d2rho*(0.5Q+0 + nbPion)*(-1 + Nf**2))/                   &
       (3.Q+0*(1 + mp2)**1.5Q+0*Sqrt(zb))))/2.Q+0

  dr3dtV=(k**4*v3*(((1 - etaphi/5.Q+0)*k**3*ms2d1rho**3*nbd3xSigma)/               &
       (12.Q+0*(1 + ms2)**2*zb**2) +                                               &
      ((1 - etaphi/5.Q+0)*k**3*mp2d1rho**3*nbd3xPion*(-1 + Nf**2))/                &
       (12.Q+0*(1 + mp2)**2*zb**2) -                                               &
      ((1 - etaphi/5.Q+0)*k**2*ms2d1rho**3*nbd2xSigma)/                            &
       (2.Q+0*(1 + ms2)**2.5Q+0*zb**1.5Q+0) +                                            &
      ((1 - etaphi/5.Q+0)*k**2*ms2d1rho*ms2d2rho*nbd2xSigma)/                      &
       (2.Q+0*(1 + ms2)**1.5Q+0*zb**1.5Q+0) -                                            &
      ((1 - etaphi/5.Q+0)*k**2*mp2d1rho**3*nbd2xPion*(-1 + Nf**2))/                &
       (2.Q+0*(1 + mp2)**2.5Q+0*zb**1.5Q+0) +                                            &
      ((1 - etaphi/5.Q+0)*k**2*mp2d1rho*mp2d2rho*nbd2xPion*(-1 + Nf**2))/          &
       (2.Q+0*(1 + mp2)**1.5Q+0*zb**1.5Q+0) +                                            &
      (5*(1 - etaphi/5.Q+0)*k*ms2d1rho**3*nbd1xSigma)/(4.Q+0*(1 + ms2)**3*zb) -       &
      (3*(1 - etaphi/5.Q+0)*k*ms2d1rho*ms2d2rho*nbd1xSigma)/                       &
       (2.Q+0*(1 + ms2)**2*zb) + ((1 - etaphi/5.Q+0)*k*ms2d3rho*nbd1xSigma)/          &
       (3.Q+0*(1 + ms2)*zb) + (5*(1 - etaphi/5.Q+0)*k*mp2d1rho**3*nbd1xPion*          &
         (-1 + Nf**2))/(4.Q+0*(1 + mp2)**3*zb) -                                   &
      (3*(1 - etaphi/5.Q+0)*k*mp2d1rho*mp2d2rho*nbd1xPion*(-1 + Nf**2))/           &
       (2.Q+0*(1 + mp2)**2*zb) + ((1 - etaphi/5.Q+0)*k*mp2d3rho*nbd1xPion*            &
         (-1 + Nf**2))/(3.Q+0*(1 + mp2)*zb) -                                      &
      (5*(1 - etaphi/5.Q+0)*ms2d1rho**3*(0.5Q+0 + nbSigma))/                          &
       (4.Q+0*(1 + ms2)**3.5Q+0*Sqrt(zb)) +                                           &
      (3*(1 - etaphi/5.Q+0)*ms2d1rho*ms2d2rho*(0.5Q+0 + nbSigma))/                    &
       (2.Q+0*(1 + ms2)**2.5Q+0*Sqrt(zb)) -                                           &
      ((1 - etaphi/5.Q+0)*ms2d3rho*(0.5Q+0 + nbSigma))/                               &
       (3.Q+0*(1 + ms2)**1.5Q+0*Sqrt(zb)) -                                           &
      (5*(1 - etaphi/5.Q+0)*mp2d1rho**3*(0.5Q+0 + nbPion)*(-1 + Nf**2))/              &
       (4.Q+0*(1 + mp2)**3.5Q+0*Sqrt(zb)) +                                           &
      (3*(1 - etaphi/5.Q+0)*mp2d1rho*mp2d2rho*(0.5Q+0 + nbPion)*(-1 + Nf**2))/        &
       (2.Q+0*(1 + mp2)**2.5Q+0*Sqrt(zb)) -                                           &
      ((1 - etaphi/5.Q+0)*mp2d3rho*(0.5Q+0 + nbPion)*(-1 + Nf**2))/                   &
       (3.Q+0*(1 + mp2)**1.5Q+0*Sqrt(zb))))/2.Q+0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  dkappadt=-etaphi*kappa

  dlam0dt=dr0dtV
  dlam1dt=dr1dtV
  dlam2dt=dr2dtV
  dlam3dt=dr3dtV
  dlam4dt=0.Q+0
  dlam5dt=0.Q+0
  dlam6dt=0.Q+0
  dlam7dt=0.Q+0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  dydx(1)=dlam1dt
  dydx(2)=dlam2dt
  dydx(3)=dlam3dt
  dydx(4)=dlam4dt
  dydx(5)=dlam5dt
  dydx(6)=dlam6dt
  dydx(7)=dlam7dt
  dydx(Nv+1)=dlam0dt
  dydx((Nv+1)+1)=0.
  dydx((Nv+1)+(Nh+1)+1)=0.
  dydx((Nv+1)+(Nh+1)+2)=0.
  dydx((Nv+1)+(Nh+1)+Nz+1)=0.
  dydx((Nv+1)+(Nh+1)+Nz+2)=0.

  write(*,*)'k=',k*197.33
  !goto 100

  open(unit=101,file='./buffer/k1.dat')
  write(101, "(e20.9)")k*hc

  open(unit=102,file='./buffer/lam1.dat')
  write(102, "(e20.9)")lam1

  open(unit=103,file='./buffer/lam2.dat')
  write(103, "(e20.9)")lam2

  open(unit=104,file='./buffer/lam3.dat')
  write(104, "(e20.9)")lam3
100 continue

end



