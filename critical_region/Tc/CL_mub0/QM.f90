program QM

  implicit none

  real(16) pi,hc
  parameter(pi=3.141592653589793238462643383279Q+0)
  parameter(hc=197.33Q+0)
  real(16) T,mu !temperature and chemical potential
  real(16) l,lb !polyakov loop
  real(16) rho0,mPion,mSigma,mf,Vtotal,fpi,h,Zphi,Zpsi,c,kappa
  real(16) sigma_UV,kappa_UV_i,kappa_UV
  real(16) dT,mudelta
  integer i,iTmax
  parameter(iTmax=300)
  real(16) T_MeV,muB_MeV
  real(16) muBi,muB
  integer mmax
  parameter(mmax=10)
  integer iT,iv
  real(16) l_i,lb_i,cl


  common /Tmu/ T,mu
  common /iTiv/ iT,iv



  common /cl_com/ cl

  dT=1.Q-5

  muBi=0.Q+0/hc

  sigma_UV=80.Q+0/hc
  kappa_UV_i=sigma_UV**2/2.Q+0
  l_i=0.396610714974989101424322305634385553Q+00!.Q-10
  lb_i=0.396610714974989101424322305634385553Q+00!1.Q-10

  mudelta=0.

  do i=1,50!iTmax
    cl=1.Q-10+(0.167003Q+00-1.Q-10)/50.*real(i-1,kind=16)!1.Q-10+1.Q-9*real(i-1,kind=16)
!    T_MeV=153.1722Q+00!+dT*real(i,kind=16)
    T_MeV=153.173000Q+0!+dT*real(i,kind=16)
    T=T_MeV/hc

    muB=muBi+mudelta
    mu=muB/3.Q+0
    muB_MeV=muB*hc

      call selfEQ(kappa_UV_i,l_i,lb_i,kappa_UV,l,lb,rho0,mPion,mSigma,mf,Vtotal,fpi,h,Zphi,Zpsi,c,kappa)
      kappa_UV_i=kappa_UV
      l_i=l
      lb_i=lb

    write(*,"('i=', I4)")i
    write(*,"('muB_MeV=', f15.7, t25, 'T_MeV=', f15.7)")muB_MeV,T_MeV
    write(*,"('fpi=', f15.7, t25, 'mPion=', f15.7)")fpi*hc,mPion*hc
    write(*,"('mSigma=', f15.7, t25, 'mf=', f15.7)")mSigma*hc,mf*hc
    write(*,"('kappa_UV=', e21.14)")kappa_UV

    open(unit=51,file='./buffer/TMeV.dat',position='append')
    write(51,*)T_MeV
    close(51)

    open(unit=51,file='./buffer/l.dat',position='append')
    write(51,*)l
    close(51)

    open(unit=51,file='./buffer/lb.dat',position='append')
    write(51,*)lb
    close(51)

    open(unit=51,file='./buffer/mPion.dat',position='append')
    write(51, "(e21.14)")mPion*hc
    close(51)

    open(unit=51,file='./buffer/mSigma.dat',position='append')
    write(51,*)mSigma*hc
    close(51)

    open(unit=51,file='./buffer/mf.dat',position='append')
    write(51,*)mf*hc
    close(51)

    open(unit=51,file='./buffer/h.dat',position='append')
    write(51,*)h
    close(51)

    open(unit=51,file='./buffer/Vtotal.dat',position='append')
    write(51,*)Vtotal
    close(51)

    open(unit=51,file='./buffer/fpi.dat',position='append')
    write(51,*)fpi*hc
    close(51)

    open(unit=51,file='./buffer/Zphi.dat',position='append')
    write(51,*)Zphi
    close(51)

    open(unit=51,file='./buffer/kappaUV.dat',position='append')
    write(51,*)kappa_UV
    close(51)

    open(unit=51,file='./buffer/cl.dat',position='append')
    write(51,*)cl
    close(51)


  end do



end






