Program FRG
!this code solve QM model by chebyshev method
!Yong-Rui Chen 2020.1.6

  implicit none
  
  real(16) pi
  parameter (pi=3.141592653589793Q+00)
  real(16) T,mu,Tup,Tdown
  real(16) Tc_save(49),mu_save(49)
  integer Tnum,munum
  real(16) K_UV,K_IR
  integer Nflow,Vmax
  !number of flow equations
  parameter (Nflow=800,Vmax=3000)
  real(16) a(Vmax),b(Vmax),c(Vmax),Tn(Vmax,Vmax),summation
  !Chebyshev coefficient
  integer(8) gridnum,maxim,anum
  !number of rho
  real(16) rho(Vmax),rhomin,rhomax,yflow(Nflow)
  external derivs,rkqs
  real(16) eps_ode,h1,hmin 
  !variables in subroutine odeint
  integer nok,nbad 
  !variables in subroutine odeint
  integer(8) kmax,kount
  !variables in common block of subroutine odeint
  integer KMAXX,NMAX
  parameter (NMAX=800,KMAXX=20000)
  real(16) dxsav,xp(KMAXX),yp(NMAX,KMAXX) 
  !variables in common block of subroutine odeint
  real(16) h,cbar,delta
  !Yukawa coupling and explicit chiral symmetry breaking 
  integer i,imin,j,k,j1
  real(16) dydx(Nflow),y(Nflow)
  real(16) V0(100000)
  real(16) funp(51,KMAXX)


  common /Tmu/ T,mu
  common /gridrho/ gridnum,anum,rho,rhomin,rhomax,maxim
  common /odeContr/ eps_ode,h1,hmin
  common /path/ kmax,kount,dxsav,xp
  common /path1/ yp
  common /path2/ funp
  common /iniconst/ h,cbar
  common /chebyshev/ Tn,b,c

  T=124.6915Q+00
  mu=0.Q+00

  k_UV=700.Q+00
  k_IR=0.01Q+00

  eps_ode=1.Q-12
  h1=(k_IR-k_UV)/5000.Q+00
  hmin=0.Q+00
  kmax=KMAXX
  dxsav=(k_IR-k_UV)/10000.Q+00

  anum=51
  gridnum=2000
  
  rhomin=0.Q+00
  rhomax=10000.Q+00
  do i=1,anum
  do j=1,gridnum
  rho(j)=rhomax/2.Q+00*cos(pi*(j-0.5Q+00)/gridnum)+rhomax/2.Q+00
  Tn(i,j)=cos((pi*(i-1)*(j-0.5Q+00)/gridnum))
  end do
  end do

  call initial(a)
  call odeint(a,anum,k_UV,k_IR,eps_ode,h1,hmin,nok,nbad,derivs,rkqs)
!  call derivs(k_UV,yflow,dydx)
  call phyphi(a,V0,imin)

  write(*,*)'T=',T
  write(*,*)'mu=',mu
!  write(*,*)'phi0=',sqrt(2.Q+00*rho(imin))
!  write(*,*)'mqk=',sqrt(h**2.Q+00*rho(imin)/2.Q+00)

  open(unit=50,file='./buffer/mu.dat')
    write(50,*)mu
  close(50)

  open(unit=50,file='./buffer/T.dat',position='append')
    write(50,*)T
  close(50)

  open(unit=50,file='./buffer/k.dat',position='append')
    write(50,*)xp(1:kount)
  close(50)

  open(unit=50,file='./buffer/bk.dat',position='append')
  do i=1,kount
     do j=1,anum
        write(50,*)funp(j,i)
     end do
  end do
  close(50)

  open(unit=50,file='./buffer/ak.dat',position='append')
  do i=1,kount
     do j=1,anum
        write(50,*)yp(j,i)
     end do
  end do
  close(50)
end
