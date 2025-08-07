subroutine cal_BB(k,mpi2,msg2,dnbpidx,dnbsgdx,BBall,gridnum)

  implicit none
  real(16) k 
  ! IR cutoff in flow equations
  integer Vmax,i
  !maximal number of potential
  parameter (Vmax=100)
  real(16) mpi2(Vmax),msg2(Vmax)
  real(16) ma2,mb2
  real(16) dnbpidx(Vmax,7),dnbsgdx(Vmax,7)
  real(16) dnbdxa0,dnbdxa1,dnbdxa2,dnbdxa3,dnbdxa4,dnbdxa5,dnbdxa6
  real(16) dnbdxb0,dnbdxb1,dnbdxb2,dnbdxb3,dnbdxb4
  real(16) B2B2,B3B2,B4B2,B5B2
  real(16) B2B3,B3B3,B4B3
  real(16) B2B4,B3B4,B2B5
  real(16) BBall(Vmax,10)
  integer(8) gridnum
  !number of rho

  do i=1,gridnum
    ma2=mpi2(i)
    mb2=msg2(i)
  
    dnbdxa0=dnbpidx(i,1)
    dnbdxa1=dnbpidx(i,2)
    dnbdxa2=dnbpidx(i,3)
    dnbdxa3=dnbpidx(i,4)
    dnbdxa4=dnbpidx(i,5)
    dnbdxa5=dnbpidx(i,6)
    dnbdxa6=dnbpidx(i,7)

    dnbdxb0=dnbsgdx(i,1)
    dnbdxb1=dnbsgdx(i,2)
    dnbdxb2=dnbsgdx(i,3)
    dnbdxb3=dnbsgdx(i,4)
    dnbdxb4=dnbsgdx(i,5)

  if(abs(ma2-mb2)<1.Q-20)then

  b2b2=(k**7.Q+00*(12.Q+00*dnbdxa2*(k**2.Q+00 + ma2) -                                                         &
        2.Q+00*dnbdxa3*(k**2.Q+00 + ma2)**1.5Q+00 +                                                            &
        15.Q+00*(1.Q+00 + 2.Q+00*dnbdxa0 - 2.Q+00*dnbdxa1*Sqrt(k**2.Q+00 + ma2))))/                            &
    (16.Q+00*(k**2.Q+00 + ma2)**3.5Q+00)
  b3b2= -(k**9.Q+00*(2.Q+00*(k**2.Q+00 + ma2)*(-3.Q+00*dnbdxa2 + 3.Q+00*dnbdxa3*Sqrt(k**2.Q+00 + ma2) -        &
            dnbdxa4*(k**2.Q+00 + ma2)) -                                                                       &
         7.Q+00*(12.Q+00*dnbdxa2*(k**2.Q+00 + ma2) - 2.Q+00*dnbdxa3*(k**2.Q+00 + ma2)**1.5Q+00 +               &
            15.Q+00*(1.Q+00 + 2.Q+00*dnbdxa0 - 2.Q+00*dnbdxa1*Sqrt(k**2.Q+00 + ma2)))))/                       &
    (64.Q+00*(k**2.Q+00 + ma2)**4.5Q+00)
  b4b2= -(k**11.Q+00*(-945.Q+00 - 1890.Q+00*dnbdxa0 - 840.Q+00*dnbdxa2*k**2.Q+00 - 30.Q+00*dnbdxa4*k**4.Q+00 - &
         840.Q+00*dnbdxa2*ma2 - 60.Q+00*dnbdxa4*k**2.Q+00*ma2 - 30.Q+00*dnbdxa4*ma2**2.Q+00 +                      &
         1890.Q+00*dnbdxa1*Sqrt(k**2.Q+00 + ma2) +                                        &
         210.Q+00*dnbdxa3*k**2.Q+00*Sqrt(k**2.Q+00 + ma2) +                                    &
         2.Q+00*dnbdxa5*k**4.Q+00*Sqrt(k**2.Q+00 + ma2) +                                      &
         210.Q+00*dnbdxa3*ma2*Sqrt(k**2.Q+00 + ma2) +                                     &
         4.Q+00*dnbdxa5*k**2.Q+00*ma2*Sqrt(k**2.Q+00 + ma2) +                                  &
         2.Q+00*dnbdxa5*ma2**2.Q+00*Sqrt(k**2.Q+00 + ma2)))/(384.Q+00*(k**2.Q+00 + ma2)**5.5Q+00)
  b5b2=(k**13.Q+00*(10395.Q+00 + 20790.Q+00*dnbdxa0 + 9450.Q+00*dnbdxa2*k**2.Q+00 +                      &
        420.Q+00*dnbdxa4*k**4.Q+00 + 2.Q+00*dnbdxa6*k**6.Q+00 + 9450.Q+00*dnbdxa2*ma2 +                  &
        840.Q+00*dnbdxa4*k**2.Q+00*ma2 + 6.Q+00*dnbdxa6*k**4.Q+00*ma2 + 420.Q+00*dnbdxa4*ma2**2.Q+00 +        &
        6.Q+00*dnbdxa6*k**2.Q+00*ma2**2.Q+00 + 2.Q+00*dnbdxa6*ma2**3.Q+00 -                              &
        20790.Q+00*dnbdxa1*Sqrt(k**2.Q+00 + ma2) -                                        &
        2520.Q+00*dnbdxa3*k**2.Q+00*Sqrt(k**2.Q+00 + ma2) -                                    &
        42.Q+00*dnbdxa5*k**4.Q+00*Sqrt(k**2.Q+00 + ma2) -                                      &
        2520.Q+00*dnbdxa3*ma2*Sqrt(k**2.Q+00 + ma2) -                                     &
        84.Q+00*dnbdxa5*k**2.Q+00*ma2*Sqrt(k**2.Q+00 + ma2) -                                  &
        42.Q+00*dnbdxa5*ma2**2.Q+00*Sqrt(k**2.Q+00 + ma2)))/(3072.Q+00*(k**2.Q+00 + ma2)**6.5Q+00)
  b2b3=-(k**9.Q+00*(2.Q+00*(k**2.Q+00 + ma2)*                                                  &
          (-3.Q+00*dnbdxa2 + 3.Q+00*dnbdxa3*Sqrt(k**2.Q+00 + ma2) - dnbdxa4*(k**2.Q+00 + ma2))      &
- 7.Q+00*(12.Q+00*dnbdxa2*(k**2.Q+00 + ma2) - 2.Q+00*dnbdxa3*(k**2.Q+00 + ma2)**1.5Q+00 +                    &
            15.Q+00*(1.Q+00 + 2.Q+00*dnbdxa0 - 2.Q+00*dnbdxa1*Sqrt(k**2.Q+00 + ma2)))))/                 &
    (64.Q+00*(k**2.Q+00 + ma2)**4.5Q+00)
  b3b3= -(k**11.Q+00*(-945.Q+00 - 1890.Q+00*dnbdxa0 - 840.Q+00*dnbdxa2*k**2.Q+00 - 30.Q+00*dnbdxa4*k**4.Q+00 -     &
         840.Q+00*dnbdxa2*ma2 - 60.Q+00*dnbdxa4*k**2.Q+00*ma2 - 30.Q+00*dnbdxa4*ma2**2.Q+00 +         &
         1890.Q+00*dnbdxa1*Sqrt(k**2.Q+00 + ma2) +                                        &
         210.Q+00*dnbdxa3*k**2.Q+00*Sqrt(k**2.Q+00 + ma2) +                                    &
         2.Q+00*dnbdxa5*k**4.Q+00*Sqrt(k**2.Q+00 + ma2) +                                      &
         210.Q+00*dnbdxa3*ma2*Sqrt(k**2.Q+00 + ma2) +                                     &
         4.Q+00*dnbdxa5*k**2.Q+00*ma2*Sqrt(k**2.Q+00 + ma2) +                                  &
         2.Q+00*dnbdxa5*ma2**2.Q+00*Sqrt(k**2.Q+00 + ma2)))/(256.Q+00*(k**2.Q+00 + ma2)**5.5Q+00)
  b4b3=(k**13.Q+00*(10395.Q+00 + 20790.Q+00*dnbdxa0 + 9450.Q+00*dnbdxa2*k**2.Q+00 +                      &
        420.Q+00*dnbdxa4*k**4.Q+00 + 2.Q+00*dnbdxa6*k**6.Q+00 + 9450.Q+00*dnbdxa2*ma2 +                  &
        840.Q+00*dnbdxa4*k**2.Q+00*ma2 + 6.Q+00*dnbdxa6*k**4.Q+00*ma2 + 420.Q+00*dnbdxa4*ma2**2.Q+00 +        &
        6.Q+00*dnbdxa6*k**2.Q+00*ma2**2.Q+00 + 2.Q+00*dnbdxa6*ma2**3.Q+00 -                              &
        20790.Q+00*dnbdxa1*Sqrt(k**2.Q+00 + ma2) -                                        &
        2520.Q+00*dnbdxa3*k**2.Q+00*Sqrt(k**2.Q+00 + ma2) -                                    &
        42.Q+00*dnbdxa5*k**4.Q+00*Sqrt(k**2.Q+00 + ma2) -                                      &
        2520.Q+00*dnbdxa3*ma2*Sqrt(k**2.Q+00 + ma2) -                                     &
        84.Q+00*dnbdxa5*k**2.Q+00*ma2*Sqrt(k**2.Q+00 + ma2) -                                  &
        42.Q+00*dnbdxa5*ma2**2.Q+00*Sqrt(k**2.Q+00 + ma2)))/(1536.Q+00*(k**2.Q+00 + ma2)**6.5Q+00)
  b2b4=-(k**11.Q+00*(-945.Q+00 - 1890.Q+00*dnbdxa0 - 840.Q+00*dnbdxa2*k**2.Q+00 -                        &
         30.Q+00*dnbdxa4*k**4.Q+00 - 840.Q+00*dnbdxa2*ma2 - 60.Q+00*dnbdxa4*k**2.Q+00*ma2 -              &
         30.Q+00*dnbdxa4*ma2**2.Q+00 + 1890.Q+00*dnbdxa1*Sqrt(k**2.Q+00 + ma2) +                    &
         210.Q+00*dnbdxa3*k**2.Q+00*Sqrt(k**2.Q+00 + ma2) +                                    &
         2.Q+00*dnbdxa5*k**4.Q+00*Sqrt(k**2.Q+00 + ma2) +                                      &
         210.Q+00*dnbdxa3*ma2*Sqrt(k**2.Q+00 + ma2) +                                     &
         4.Q+00*dnbdxa5*k**2.Q+00*ma2*Sqrt(k**2.Q+00 + ma2) +                                  &
         2.Q+00*dnbdxa5*ma2**2.Q+00*Sqrt(k**2.Q+00 + ma2)))/(384.Q+00*(k**2.Q+00 + ma2)**5.5Q+00)
  b3b4=(k**13.Q+00*(10395.Q+00 + 20790.Q+00*dnbdxa0 + 9450.Q+00*dnbdxa2*k**2.Q+00 +                      &
        420.Q+00*dnbdxa4*k**4.Q+00 + 2.Q+00*dnbdxa6*k**6.Q+00 + 9450.Q+00*dnbdxa2*ma2 +                  &
        840.Q+00*dnbdxa4*k**2.Q+00*ma2 + 6.Q+00*dnbdxa6*k**4.Q+00*ma2 + 420.Q+00*dnbdxa4*ma2**2.Q+00 +        &
        6.Q+00*dnbdxa6*k**2.Q+00*ma2**2.Q+00 + 2.Q+00*dnbdxa6*ma2**3.Q+00 -                              &
        20790.Q+00*dnbdxa1*Sqrt(k**2.Q+00 + ma2) -                                        &
        2520.Q+00*dnbdxa3*k**2.Q+00*Sqrt(k**2.Q+00 + ma2) -                                    &
        42.Q+00*dnbdxa5*k**4.Q+00*Sqrt(k**2.Q+00 + ma2) -                                      &
        2520.Q+00*dnbdxa3*ma2*Sqrt(k**2.Q+00 + ma2) -                                     &
        84.Q+00*dnbdxa5*k**2.Q+00*ma2*Sqrt(k**2.Q+00 + ma2) -                                  &
        42.Q+00*dnbdxa5*ma2**2.Q+00*Sqrt(k**2.Q+00 + ma2)))/(1536.Q+00*(k**2.Q+00 + ma2)**6.5Q+00)
  b2b5=(k**13.Q+00*(10395.Q+00 + 20790.Q+00*dnbdxa0 + 9450.Q+00*dnbdxa2*k**2.Q+00 +                      &
        420.Q+00*dnbdxa4*k**4.Q+00 + 2.Q+00*dnbdxa6*k**6.Q+00 + 9450.Q+00*dnbdxa2*ma2 +                  &
        840.Q+00*dnbdxa4*k**2.Q+00*ma2 + 6.Q+00*dnbdxa6*k**4.Q+00*ma2 + 420.Q+00*dnbdxa4*ma2**2.Q+00 +        &
        6.Q+00*dnbdxa6*k**2.Q+00*ma2**2.Q+00 + 2*dnbdxa6*ma2**3.Q+00 -                              &
        20790.Q+00*dnbdxa1*Sqrt(k**2.Q+00 + ma2) -                                        &
        2520.Q+00*dnbdxa3*k**2*Sqrt(k**2.Q+00 + ma2) -                                    &
        42.Q+00*dnbdxa5*k**4.Q+00*Sqrt(k**2.Q+00 + ma2) -                                      &
        2520.Q+00*dnbdxa3*ma2*Sqrt(k**2.Q+00 + ma2) -                                     &
        84.Q+00*dnbdxa5*k**2.Q+00*ma2*Sqrt(k**2.Q+00 + ma2) -                                  &
        42.Q+00*dnbdxa5*ma2**2.Q+00*Sqrt(k**2.Q+00 + ma2)))/(3072.Q+00*(k**2.Q+00 + ma2)**6.5Q+00) 

  else

  b2b2=(k**7.Q+00*((2.Q+00 + 4.Q+00*dnbdxa0)/(Sqrt(k**2.Q+00 + ma2)*(ma2 - mb2)**3.Q+00) +               &
        (0.5Q+00 + dnbdxa0)/((k**2.Q+00 + ma2)**1.5Q+00*(ma2 - mb2)**2.Q+00) -                    &
        dnbdxa1/((k**2.Q+00 + ma2)*(ma2 - mb2)**2.Q+00) +                                 &
        (0.5Q+00 + dnbdxb0)/((ma2 - mb2)**2.Q+00*(k**2.Q+00 + mb2)**1.5Q+00) -                    &
        dnbdxb1/((ma2 - mb2)**2.Q+00*(k**2.Q+00 + mb2)) +                                 &
        (2.Q+00 + 4.Q+00*dnbdxb0)/(Sqrt(k**2.Q+00 + mb2)*(-ma2 + mb2)**3.Q+00)))/2.Q+00
  b3b2=-(k**9.Q+00*((-6.Q+00*(2.Q+00 + 4.Q+00*dnbdxa0))/(Sqrt(k**2.Q+00 + ma2)*(ma2 - mb2)**4.Q+00) -         &
         (4.Q+00*(1.Q+00 + 2.Q+00*dnbdxa0))/((k**2.Q+00 + ma2)**1.5Q+00*(ma2 - mb2)**3.Q+00) +               &
         (8.Q+00*dnbdxa1)/((k**2.Q+00 + ma2)*(ma2 - mb2)**3.Q+00) -                            &
         (3.Q+00*(0.5Q+00 + dnbdxa0))/((k**2.Q+00 + ma2)**2.5Q+00*(ma2 - mb2)**2.Q+00) +               &
         (3.Q+00*dnbdxa1)/((k**2.Q+00 + ma2)**2.Q+00*(ma2 - mb2)**2.Q+00) -                         &
         dnbdxa2/((k**2.Q+00 + ma2)**1.5Q+00*(ma2 - mb2)**2.Q+00) +                           &
         (4.Q+00*dnbdxb1)/((ma2 - mb2)**3.Q+00*(k**2.Q+00 + mb2)) +                            &
         (6.Q+00*(2.Q+00 + 4.Q+00*dnbdxb0))/((ma2 - mb2)**4.Q+00*Sqrt(k**2.Q+00 + mb2)) +                &
         (2.Q+00 + 4.Q+00*dnbdxb0)/((k**2.Q+00 + mb2)**1.5Q+00*(-ma2 + mb2)**3.Q+00)))/8.Q+00
  b4b2=(k**11.Q+00*((48.Q+00*(1.Q+00 + 2.Q+00*dnbdxa0))/(Sqrt(k**2.Q+00 + ma2)*(ma2 - mb2)**5.Q+00) +         &
        (18.Q+00*(1.Q+00 + 2.Q+00*dnbdxa0))/((k**2.Q+00 + ma2)**1.5Q+00*(ma2 - mb2)**4.Q+00) -               &
        (36.Q+00*dnbdxa1)/((k**2.Q+00 + ma2)*(ma2 - mb2)**4.Q+00) +                            &
        (6.Q+00*(1.Q+00 + 2.Q+00*dnbdxa0))/((k**2.Q+00 + ma2)**2.5Q+00*(ma2 - mb2)**3.Q+00) +                &
        (3.Q+00 + 6.Q+00*dnbdxa0)/((k**2.Q+00 + ma2)**2.5Q+00*(ma2 - mb2)**3.Q+00) -                    &
        (18.Q+00*dnbdxa1)/((k**2.Q+00 + ma2)**2.Q+00*(ma2 - mb2)**3.Q+00) +                         &
        (6.Q+00*dnbdxa2)/((k**2.Q+00 + ma2)**1.5Q+00*(ma2 - mb2)**3.Q+00) +                        &
        (15.Q+00*(1.Q+00 + 2.Q+00*dnbdxa0))/(4.Q+00*(k**2.Q+00 + ma2)**3.5Q+00*(ma2 - mb2)**2.Q+00) -            &
        (15.Q+00*dnbdxa1)/(2.Q+00*(k**2.Q+00 + ma2)**3.Q+00*(ma2 - mb2)**2.Q+00) +                      &
        (3.Q+00*dnbdxa2)/((k**2.Q+00 + ma2)**2.5Q+00*(ma2 - mb2)**2.Q+00) -                        &
        dnbdxa3/(2.Q+00*(k**2.Q+00 + ma2)**2.Q+00*(ma2 - mb2)**2.Q+00) +                           &
        (6.Q+00*(1.Q+00 + 2.Q+00*dnbdxb0))/((ma2 - mb2)**4.Q+00*(k**2.Q+00 + mb2)**1.5Q+00) -                &
        (12.Q+00*dnbdxb1)/((ma2 - mb2)**4.Q+00*(k**2.Q+00 + mb2)) +                            &
        (48.Q+00*(1.Q+00 + 2.Q+00*dnbdxb0))/(Sqrt(k**2.Q+00 + mb2)*(-ma2 + mb2)**5.Q+00)))/24.Q+00
  b5b2=-(k**13.Q+00*((-240.Q+00*(1.Q+00 + 2.Q+00*dnbdxa0))/                                         &
          (Sqrt(k**2.Q+00 + ma2)*(ma2 - mb2)**6.Q+00) -                                   &
         (96.Q+00*(1.Q+00 + 2.Q+00*dnbdxa0))/((k**2.Q+00 + ma2)**1.5Q+00*(ma2 - mb2)**5.Q+00) +              &
         (192.Q+00*dnbdxa1)/((k**2.Q+00 + ma2)*(ma2 - mb2)**5.Q+00) -                          &
         (54.Q+00*(1.Q+00 + 2.Q+00*dnbdxa0))/((k**2.Q+00 + ma2)**2.5Q+00*(ma2 - mb2)**4.Q+00) +              &
         (108.Q+00*dnbdxa1)/((k**2.Q+00 + ma2)**2.Q+00*(ma2 - mb2)**4.Q+00) -                       &
         (36.Q+00*dnbdxa2)/((k**2.Q+00 + ma2)**1.5Q+00*(ma2 - mb2)**4.Q+00) -                      &
         (45.Q+00*(1.Q+00 + 2.Q+00*dnbdxa0))/(2.Q+00*(k**2.Q+00 + ma2)**3.5Q+00*(ma2 - mb2)**3.Q+00) -           &
         (5.Q+00*(3.Q+00 + 6.Q+00*dnbdxa0))/(2.Q+00*(k**2.Q+00 + ma2)**3.5Q+00*(ma2 - mb2)**3.Q+00) +            &
         (60.Q+00*dnbdxa1)/((k**2.Q+00 + ma2)**3.Q+00*(ma2 - mb2)**3.Q+00) -                        &
         (24.Q+00*dnbdxa2)/((k**2.Q+00 + ma2)**2.5Q+00*(ma2 - mb2)**3.Q+00) +                      &
         (4.Q+00*dnbdxa3)/((k**2.Q+00 + ma2)**2.Q+00*(ma2 - mb2)**3.Q+00) -                         &
         (105.Q+00*(1.Q+00 + 2.Q+00*dnbdxa0))/(8.Q+00*(k**2.Q+00 + ma2)**4.5Q+00*(ma2 - mb2)**2.Q+00) +          &
         (105.Q+00*dnbdxa1)/(4.Q+00*(k**2.Q+00 + ma2)**4.Q+00*(ma2 - mb2)**2.Q+00) -                    &
         (45.Q+00*dnbdxa2)/(4.Q+00*(k**2.Q+00 + ma2)**3.5Q+00*(ma2 - mb2)**2.Q+00) +                   &
         (5.Q+00*dnbdxa3)/(2.Q+00*(k**2.Q+00 + ma2)**3.Q+00*(ma2 - mb2)**2.Q+00) -                      &
         dnbdxa4/(4.Q+00*(k**2.Q+00 + ma2)**2.5Q+00*(ma2 - mb2)**2.Q+00) +                        &
         (48.Q+00*dnbdxb1)/((ma2 - mb2)**5.Q+00*(k**2.Q+00 + mb2)) +                           &
         (240.Q+00*(1.Q+00 + 2.Q+00*dnbdxb0))/((ma2 - mb2)**6.Q+00*Sqrt(k**2.Q+00 + mb2)) +              &
         (24.Q+00*(1.Q+00 + 2.Q+00*dnbdxb0))/((k**2.Q+00 + mb2)**1.5Q+00*(-ma2 + mb2)**5.Q+00)))/96.Q+00
  b2b3=-(k**9.Q+00*((6.Q+00*(2.Q+00 + 4.Q+00*dnbdxa0))/(Sqrt(k**2.Q+00 + ma2)*(ma2 - mb2)**4.Q+00) +          &
         (2.Q+00 + 4.Q+00*dnbdxa0)/((k**2.Q+00 + ma2)**1.5Q+00*(ma2 - mb2)**3.Q+00) -                   &
         (4.Q+00*dnbdxa1)/((k**2.Q+00 + ma2)*(ma2 - mb2)**3.Q+00) -                            &
         (3.Q+00*(0.5Q+00 + dnbdxb0))/((ma2 - mb2)**2.Q+00*(k**2.Q+00 + mb2)**2.5Q+00) +               &
         (3.Q+00*dnbdxb1)/((ma2 - mb2)**2.Q+00*(k**2.Q+00 + mb2)**2.Q+00) +                         &
         (4.Q+00*(0.5Q+00 + dnbdxb0))/((ma2 - mb2)**3.Q+00*(k**2.Q+00 + mb2)**1.5Q+00) -               &
         dnbdxb2/((ma2 - mb2)**2.Q+00*(k**2.Q+00 + mb2)**1.5Q+00) -                           &
         (6.Q+00*(2.Q+00 + 4.Q+00*dnbdxb0))/((ma2 - mb2)**4.Q+00*Sqrt(k**2.Q+00 + mb2)) -                &
         (2.Q+00*(1.Q+00 + 2.Q+00*dnbdxb0))/((k**2.Q+00 + mb2)**1.5Q+00*(-ma2 + mb2)**3.Q+00) +              &
         (8.Q+00*dnbdxb1)/((k**2.Q+00 + mb2)*(-ma2 + mb2)**3.Q+00)))/8.Q+00
  b3b3=(k**11.Q+00*((-48.Q+00*(1.Q+00 + 2.Q+00*dnbdxa0))/                                           &
         (Sqrt(k**2.Q+00 + ma2)*(ma2 - mb2)**5.Q+00) -                                    &
        (6.Q+00*(2.Q+00 + 4.Q+00*dnbdxa0))/((k**2.Q+00 + ma2)**1.5Q+00*(ma2 - mb2)**4.Q+00) +                &
        (24.Q+00*dnbdxa1)/((k**2.Q+00 + ma2)*(ma2 - mb2)**4.Q+00) -                            &
        (3.Q+00*(1.Q+00 + 2.Q+00*dnbdxa0))/((k**2.Q+00 + ma2)**2.5Q+00*(ma2 - mb2)**3.Q+00) +                &
        (6.Q+00*dnbdxa1)/((k**2.Q+00 + ma2)**2.Q+00*(ma2 - mb2)**3.Q+00) -                          &
        (2.Q+00*dnbdxa2)/((k**2.Q+00 + ma2)**1.5Q+00*(ma2 - mb2)**3.Q+00) +                        &
        (6.Q+00*(0.5Q+00 + dnbdxb0))/((ma2 - mb2)**3.Q+00*(k**2.Q+00 + mb2)**2.5Q+00) -                &
        (12.Q+00*(1.Q+00 + 2.Q+00*dnbdxb0))/((ma2 - mb2)**4.Q+00*(k**2.Q+00 + mb2)**1.5Q+00) +               &
        (2.Q+00*dnbdxb2)/((ma2 - mb2)**3.Q+00*(k**2.Q+00 + mb2)**1.5Q+00) +                        &
        (24.Q+00*dnbdxb1)/((ma2 - mb2)**4.Q+00*(k**2.Q+00 + mb2)) +                            &
        (24.Q+00*(2.Q+00 + 4.Q+00*dnbdxb0))/((ma2 - mb2)**5.Q+00*Sqrt(k**2.Q+00 + mb2)) +                &
        (6.Q+00*dnbdxb1)/((k**2.Q+00 + mb2)**2.Q+00*(-ma2 + mb2)**3.Q+00)))/16.Q+00
  b4b3=-(k**13.Q+00*((240.Q+00*(1.Q+00 + 2.Q+00*dnbdxa0))/                                          &
          (Sqrt(k**2.Q+00 + ma2)*(ma2 - mb2)**6.Q+00) +                                   &
         (72.Q+00*(1.Q+00 + 2.Q+00*dnbdxa0))/((k**2.Q+00 + ma2)**1.5Q+00*(ma2 - mb2)**5.Q+00) -              &
         (144.Q+00*dnbdxa1)/((k**2.Q+00 + ma2)*(ma2 - mb2)**5.Q+00) +                          &
         (9.Q+00*(1.Q+00 + 2.Q+00*dnbdxa0))/((k**2.Q+00 + ma2)**2.5Q+00*(ma2 - mb2)**4.Q+00) +               &
         (9.Q+00*(2.Q+00 + 4.Q+00*dnbdxa0))/((k**2.Q+00 + ma2)**2.5Q+00*(ma2 - mb2)**4.Q+00) -               &
         (54.Q+00*dnbdxa1)/((k**2.Q+00 + ma2)**2.Q+00*(ma2 - mb2)**4.Q+00) +                        &
         (18.Q+00*dnbdxa2)/((k**2.Q+00 + ma2)**1.5Q+00*(ma2 - mb2)**4.Q+00) +                      &
         (15.Q+00*(1.Q+00 + 2.Q+00*dnbdxa0))/(2.Q+00*(k**2.Q+00 + ma2)**3.5Q+00*(ma2 - mb2)**3.Q+00) -           &
         (15.Q+00*dnbdxa1)/((k**2.Q+00 + ma2)**3.Q+00*(ma2 - mb2)**3.Q+00) +                        &
         (6.Q+00*dnbdxa2)/((k**2.Q+00 + ma2)**2.5Q+00*(ma2 - mb2)**3.Q+00) -                       &
         dnbdxa3/((k**2.Q+00 + ma2)**2.Q+00*(ma2 - mb2)**3.Q+00) -                             &
         (9.Q+00*(1.Q+00 + 2.Q+00*dnbdxb0))/((ma2 - mb2)**4.Q+00*(k**2.Q+00 + mb2)**2.5Q+00) +               &
         (18.Q+00*dnbdxb1)/((ma2 - mb2)**4.Q+00*(k**2.Q+00 + mb2)**2.Q+00) +                        &
         (48.Q+00*(1.Q+00 + 2.Q+00*dnbdxb0))/((ma2 - mb2)**5.Q+00*(k**2.Q+00 + mb2)**1.5Q+00) -              &
         (6.Q+00*dnbdxb2)/((ma2 - mb2)**4.Q+00*(k**2.Q+00 + mb2)**1.5Q+00) -                       &
         (240.Q+00*(1.Q+00 + 2.Q+00*dnbdxb0))/((ma2 - mb2)**6.Q+00*Sqrt(k**2.Q+00 + mb2)) +              &
         (96.Q+00*dnbdxb1)/((k**2.Q+00 + mb2)*(-ma2 + mb2)**5.Q+00)))/48.Q+00
  b2b4=(k**11.Q+00*(15.Q+00*(0.5Q+00 + dnbdxb0)*(k**2.Q+00 + ma2)**1.5Q+00*(ma2 - mb2)**3.Q+00 -            &
        15.Q+00*dnbdxb1*(k**2.Q+00 + ma2)**1.5Q+00*(ma2 - mb2)**3.Q+00*Sqrt(k**2.Q+00 + mb2) -          &
        18.Q+00*(1.Q+00 + 2.Q+00*dnbdxb0)*(k**2.Q+00 + ma2)**1.5Q+00*(ma2 - mb2)**2.Q+00*(k**2.Q+00 + mb2) +      &
        6.Q+00*dnbdxb2*(k**2.Q+00 + ma2)**1.5Q+00*(ma2 - mb2)**3.Q+00*(k**2.Q+00 + mb2) +               &
        36.Q+00*dnbdxb1*(k**2.Q+00 + ma2)**1.5Q+00*(ma2 - mb2)**2.Q+00*(k**2.Q+00 + mb2)**1.5Q+00 -         &
        dnbdxb3*(k**2.Q+00 + ma2)**1.5Q+00*(ma2 - mb2)**3.Q+00*(k**2.Q+00 + mb2)**1.5Q+00 +            &
        24.Q+00*(1.Q+00 + 2.Q+00*dnbdxb0)*(k**2.Q+00 + ma2)**1.5Q+00*(ma2 - mb2)*(k**2.Q+00 + mb2)**2.Q+00 +      &
        6.Q+00*(2.Q+00 + 4.Q+00*dnbdxb0)*(k**2.Q+00 + ma2)**1.5Q+00*(ma2 - mb2)*(k**2.Q+00 + mb2)**2.Q+00 -       &
        12.Q+00*dnbdxb2*(k**2.Q+00 + ma2)**1.5Q+00*(ma2 - mb2)**2.Q+00*(k**2.Q+00 + mb2)**2.Q+00 -           &
        72.Q+00*dnbdxb1*(k**2.Q+00 + ma2)**1.5Q+00*(ma2 - mb2)*(k**2.Q+00 + mb2)**2.5Q+00 -            &
        96.Q+00*(1.Q+00 + 2.Q+00*dnbdxb0)*(k**2.Q+00 + ma2)**1.5Q+00*(k**2.Q+00 + mb2)**3.Q+00 +                  &
        96.Q+00*(1.Q+00 + 2.Q+00*dnbdxa0)*(k**2.Q+00 + ma2)*(k**2.Q+00 + mb2)**3.5Q+00 +                     &
        6.Q+00*(2.Q+00 + 4.Q+00*dnbdxa0)*(ma2 - mb2)*(k**2.Q+00 + mb2)**3.5Q+00 -                       &
        24.Q+00*dnbdxa1*Sqrt(k**2.Q+00 + ma2)*(ma2 - mb2)*(k**2.Q+00 + mb2)**3.5Q+00))/            &
    (48.Q+00*(k**2.Q+00 + ma2)**1.5Q+00*(ma2 - mb2)**5.Q+00*(k**2.Q+00 + mb2)**3.5Q+00)
  b3b4=-(k**13.Q+00*(-10.Q+00*(k**2.Q+00 + ma2)*                                               &
          (15.Q+00*(0.5Q+00 + dnbdxb0)*(k**2.Q+00 + ma2)**1.5Q+00*(ma2 - mb2)**3.Q+00 -                &
            15.Q+00*dnbdxb1*(k**2.Q+00 + ma2)**1.5Q+00*(ma2 - mb2)**3.Q+00*                        &
             Sqrt(k**2.Q+00 + mb2) -                                                 &
            18.Q+00*(1.Q+00 + 2.Q+00*dnbdxb0)*(k**2.Q+00 + ma2)**1.5Q+00*(ma2 - mb2)**2.Q+00*                &
             (k**2.Q+00 + mb2) + 6.Q+00*dnbdxb2*(k**2.Q+00 + ma2)**1.5Q+00*(ma2 - mb2)**3.Q+00*         &
             (k**2.Q+00 + mb2) + 36.Q+00*dnbdxb1*(k**2.Q+00 + ma2)**1.5Q+00*(ma2 - mb2)**2.Q+00*        &
             (k**2.Q+00 + mb2)**1.5Q+00 -                                                &
            dnbdxb3*(k**2.Q+00 + ma2)**1.5Q+00*(ma2 - mb2)**3.Q+00*(k**2.Q+00 + mb2)**1.5Q+00 +        &
            24.Q+00*(1.Q+00 + 2.Q+00*dnbdxb0)*(k**2.Q+00 + ma2)**1.5Q+00*(ma2 - mb2)*                   &
             (k**2.Q+00 + mb2)**2.Q+00 +                                                  &
            6.Q+00*(2.Q+00 + 4.Q+00*dnbdxb0)*(k**2.Q+00 + ma2)**1.5Q+00*(ma2 - mb2)*                    &
             (k**2.Q+00 + mb2)**2.Q+00 -                                                  &
            12.Q+00*dnbdxb2*(k**2.Q+00 + ma2)**1.5Q+00*(ma2 - mb2)**2.Q+00*(k**2.Q+00 + mb2)**2.Q+00 -       &
            72.Q+00*dnbdxb1*(k**2.Q+00 + ma2)**1.5Q+00*(ma2 - mb2)*(k**2.Q+00 + mb2)**2.5Q+00 -        &
            96.Q+00*(1.Q+00 + 2.Q+00*dnbdxb0)*(k**2.Q+00 + ma2)**1.5Q+00*(k**2.Q+00 + mb2)**3.Q+00 +              &
            96.Q+00*(1.Q+00 + 2.Q+00*dnbdxa0)*(k**2.Q+00 + ma2)*(k**2.Q+00 + mb2)**3.5Q+00 +                 &
            6.Q+00*(2.Q+00 + 4.Q+00*dnbdxa0)*(ma2 - mb2)*(k**2.Q+00 + mb2)**3.5Q+00 -                   &
            24.Q+00*dnbdxa1*Sqrt(k**2.Q+00 + ma2)*(ma2 - mb2)*(k**2.Q+00 + mb2)**3.5Q+00) -        &
         3.Q+00*(ma2 - mb2)*(15.Q+00*(0.5Q+00 + dnbdxb0)*(k**2.Q+00 + ma2)**1.5Q+00*                   &
             (ma2 - mb2)**3.Q+00 -                                                   &
            15.Q+00*dnbdxb1*(k**2.Q+00 + ma2)**1.5Q+00*(ma2 - mb2)**3.Q+00*                        &
             Sqrt(k**2.Q+00 + mb2) -                                                 &
            18.Q+00*(1.Q+00 + 2.Q+00*dnbdxb0)*(k**2.Q+00 + ma2)**1.5Q+00*(ma2 - mb2)**2.Q+00*                &
             (k**2.Q+00 + mb2) + 6.Q+00*dnbdxb2*(k**2.Q+00 + ma2)**1.5Q+00*(ma2 - mb2)**3.Q+00*         &
             (k**2.Q+00 + mb2) + 36.Q+00*dnbdxb1*(k**2.Q+00 + ma2)**1.5Q+00*(ma2 - mb2)**2.Q+00*        &
             (k**2.Q+00 + mb2)**1.5Q+00 -                                                &
            dnbdxb3*(k**2.Q+00 + ma2)**1.5Q+00*(ma2 - mb2)**3.Q+00*(k**2.Q+00 + mb2)**1.5Q+00 +        &
            24.Q+00*(1.Q+00 + 2.Q+00*dnbdxb0)*(k**2.Q+00 + ma2)**1.5Q+00*(ma2 - mb2)*                   &
             (k**2.Q+00 + mb2)**2.Q+00 +                                                  &
            6.Q+00*(2.Q+00 + 4.Q+00*dnbdxb0)*(k**2.Q+00 + ma2)**1.5Q+00*(ma2 - mb2)*                    &
             (k**2.Q+00 + mb2)**2.Q+00 -                                                  &
            12.Q+00*dnbdxb2*(k**2.Q+00 + ma2)**1.5Q+00*(ma2 - mb2)**2.Q+00*(k**2.Q+00 + mb2)**2.Q+00 -       &
            72.Q+00*dnbdxb1*(k**2.Q+00 + ma2)**1.5Q+00*(ma2 - mb2)*(k**2.Q+00 + mb2)**2.5Q+00 -        &
            96.Q+00*(1.Q+00 + 2.Q+00*dnbdxb0)*(k**2.Q+00 + ma2)**1.5Q+00*(k**2.Q+00 + mb2)**3.Q+00 +              &
            96.Q+00*(1.Q+00 + 2.Q+00*dnbdxa0)*(k**2.Q+00 + ma2)*(k**2.Q+00 + mb2)**3.5Q+00 +                 &
            6.Q+00*(2.Q+00 + 4.Q+00*dnbdxa0)*(ma2 - mb2)*(k**2.Q+00 + mb2)**3.5Q+00 -                   &
            24.Q+00*dnbdxa1*Sqrt(k**2.Q+00 + ma2)*(ma2 - mb2)*(k**2.Q+00 + mb2)**3.5Q+00) +        &
         3.Q+00*(k**2.Q+00 + ma2)*(ma2 - mb2)*                                            &
          (15.Q+00*(1.Q+00 + 2.Q+00*dnbdxb0)*(k**2.Q+00 + ma2)**1.5Q+00*(ma2 - mb2)**2.Q+00 +                &
            15.Q+00*(0.5Q+00 + dnbdxb0)*Sqrt(k**2.Q+00 + ma2)*(ma2 - mb2)**3.Q+00 -                &
            30.Q+00*dnbdxb1*(k**2.Q+00 + ma2)**1.5Q+00*(ma2 - mb2)**2.Q+00*Sqrt(k**2.Q+00 + mb2) -      &
            15.Q+00*dnbdxb1*Sqrt(k**2.Q+00 + ma2)*(ma2 - mb2)**3.Q+00*Sqrt(k**2.Q+00 + mb2) -       &
            24.Q+00*(1.Q+00 + 2.Q+00*dnbdxb0)*(k**2.Q+00 + ma2)**1.5Q+00*(ma2 - mb2)*                   &
             (k**2.Q+00 + mb2) - 18.Q+00*(1.Q+00 + 2.Q+00*dnbdxb0)*Sqrt(k**2.Q+00 + ma2)*                &
             (ma2 - mb2)**2.Q+00*(k**2.Q+00 + mb2) +                                      &
            12.Q+00*dnbdxb2*(k**2.Q+00 + ma2)**1.5Q+00*(ma2 - mb2)**2.Q+00*(k**2.Q+00 + mb2) +          &
            6.Q+00*dnbdxb2*Sqrt(k**2.Q+00 + ma2)*(ma2 - mb2)**3.Q+00*(k**2.Q+00 + mb2) +            &
            48.Q+00*dnbdxb1*(k**2.Q+00 + ma2)**1.5Q+00*(ma2 - mb2)*(k**2.Q+00 + mb2)**1.5Q+00 +        &
            36.Q+00*dnbdxb1*Sqrt(k**2.Q+00 + ma2)*(ma2 - mb2)**2.Q+00*(k**2.Q+00 + mb2)**1.5Q+00 -      &
            2.Q+00*dnbdxb3*(k**2.Q+00 + ma2)**1.5Q+00*(ma2 - mb2)**2.Q+00*(k**2.Q+00 + mb2)**1.5Q+00 -      &
            dnbdxb3*Sqrt(k**2.Q+00 + ma2)*(ma2 - mb2)**3.Q+00*(k**2.Q+00+ mb2)**1.5Q+00 +         &
            24.Q+00*(1.Q+00 + 2.Q+00*dnbdxb0)*(k**2.Q+00 + ma2)**1.5Q+00*(k**2.Q+00 + mb2)**2.Q+00 +              &
            24.Q+00*(1.Q+00 + 2.Q+00*dnbdxb0)*Sqrt(k**2.Q+00 + ma2)*(ma2 - mb2)*                    &
             (k**2.Q+00 + mb2)**2.Q+00 +                                                  &
            6.Q+00*(2.Q+00 + 4.Q+00*dnbdxb0)*Sqrt(k**2.Q+00 + ma2)*(ma2 - mb2)*                     &
             (k**2.Q+00 + mb2)**2.Q+00 -                                                  &
            16.Q+00*dnbdxb2*(k**2.Q+00 + ma2)**1.5Q+00*(ma2 - mb2)*(k**2.Q+00 + mb2)**2.Q+00 -          &
            12.Q+00*dnbdxb2*Sqrt(k**2.Q+00 + ma2)*(ma2 - mb2)**2.Q+00*(k**2.Q+00 + mb2)**2.Q+00 -        &
            48.Q+00*dnbdxb1*(k**2.Q+00 + ma2)**1.5Q+00*(k**2.Q+00 + mb2)**2.5Q+00 -                    &
            72.Q+00*dnbdxb1*Sqrt(k**2.Q+00 + ma2)*(ma2 - mb2)*(k**2.Q+00 + mb2)**2.5Q+00 -         &
            96.Q+00*(1.Q+00 + 2.Q+00*dnbdxb0)*Sqrt(k**2.Q+00 + ma2)*(k**2.Q+00 + mb2)**3.Q+00 +               &
            72.Q+00*(1.Q+00 + 2.Q+00*dnbdxa0)*(k**2.Q+00 + mb2)**3.5Q+00 +                              &
            48.Q+00*dnbdxa1*Sqrt(k**2.Q+00 + ma2)*(k**2.Q+00 + mb2)**3.5Q+00 +                     &
            8.Q+00*dnbdxa2*(k**2.Q+00 + mb2)**3.5Q+00*(-ma2 + mb2))))/                        &
    (192.Q+00*(k**2.Q+00 + ma2)**2.5Q+00*(ma2 - mb2)**6.Q+00*(k**2.Q+00 + mb2)**3.5Q+00)
  b2b5=-(k**13.Q+00*(2.Q+00*(ma2 - mb2)*(k**2.Q+00 + mb2)*                                     &
          (-45.Q+00*(0.5Q+00 + dnbdxb0)*(k**2.Q+00 + ma2)**1.5Q+00*(ma2 - mb2)**2.Q+00 -               &
            18.Q+00*(1.Q+00 + 2.Q+00*dnbdxb0)*(k**2.Q+00 + ma2)**1.5Q+00*(ma2 - mb2)**2.Q+00 -               &
            (3.Q+00*dnbdxb2*(k**2.Q+00 + ma2)**1.5Q+00*(ma2 - mb2)**3.Q+00)/2.Q+00 +                   &
            81.Q+00*dnbdxb1*(k**2.Q+00 + ma2)**1.5Q+00*(ma2 - mb2)**2.Q+00*Sqrt(k**2.Q+00 + mb2) +      &
            (3.Q+00*dnbdxb3*(k**2.Q+00 + ma2)**1.5Q+00*(ma2 - mb2)**3.Q+00*Sqrt(k**2.Q+00 + mb2))/      &
             2.Q+00 + 108.Q+00*(1.Q+00 + 2.Q+00*dnbdxb0)*(k**2.Q+00 + ma2)**1.5Q+00*(ma2 - mb2)*            &
             (k**2.Q+00 + mb2) - 24.Q+00*dnbdxb2*(k**2.Q+00 + ma2)**1.5Q+00*(ma2 - mb2)**2.Q+00*        &
             (k**2.Q+00 + mb2) - (dnbdxb4*(k**2.Q+00 + ma2)**1.5Q+00*(ma2 - mb2)**3.Q+00*          &
               (k**2.Q+00 + mb2))/2.Q+00 -                                               &
            216.Q+00*dnbdxb1*(k**2.Q+00 + ma2)**1.5Q+00*(ma2 - mb2)*(k**2.Q+00 + mb2)**1.5Q+00 -       &
            3.Q+00*dnbdxb3*(k**2.Q+00 + ma2)**1.5Q+00*(ma2 - mb2)**2.Q+00*(k**2.Q+00 + mb2)**1.5Q+00 -      &
            312.Q+00*(1.Q+00 + 2.Q+00*dnbdxb0)*(k**2.Q+00 + ma2)**1.5Q+00*(k**2.Q+00 + mb2)**2.Q+00 -             &
            6.Q+00*(2.Q+00 + 4.Q+00*dnbdxb0)*(k**2.Q+00 + ma2)**1.5Q+00*(k**2.Q+00 + mb2)**2.Q+00 -               &
            12.Q+00*dnbdxb2*(k**2.Q+00 + ma2)**1.5Q+00*(ma2 - mb2)*(k**2.Q+00 + mb2)**2.Q+00 +          &
            336.Q+00*(1.Q+00 + 2.Q+00*dnbdxa0)*(k**2.Q+00 + ma2)*(k**2.Q+00 + mb2)**2.5Q+00 -                &
            24.Q+00*dnbdxb1*(k**2.Q+00 + ma2)**1.5Q+00*(k**2.Q+00 + mb2)**2.5Q+00 +                    &
            21.Q+00*(2.Q+00 + 4.Q+00*dnbdxa0)*(ma2 - mb2)*(k**2.Q+00 + mb2)**2.5Q+00 -                  &
            84.Q+00*dnbdxa1*Sqrt(k**2.Q+00 + ma2)*(ma2 - mb2)*(k**2.Q+00 + mb2)**2.5Q+00 -         &
            6.Q+00*(2.Q+00 + 4.Q+00*dnbdxa0)*(k**2.Q+00 + mb2)**3.5Q+00 +                               &
            24.Q+00*dnbdxa1*Sqrt(k**2.Q+00 + ma2)*(k**2.Q+00 + mb2)**3.5Q+00) -                    &
         7.Q+00*(ma2 - mb2)*(15.Q+00*(0.5Q+00 + dnbdxb0)*(k**2.Q+00 + ma2)**1.5Q+00*                   &
             (ma2 - mb2)**3.Q+00 - 15.Q+00*dnbdxb1*(k**2.Q+00 + ma2)**1.5Q+00*(ma2 - mb2)**3.Q+00*      &
             Sqrt(k**2.Q+00 + mb2) -                                                 &
            18.Q+00*(1.Q+00 + 2.Q+00*dnbdxb0)*(k**2.Q+00 + ma2)**1.5Q+00*(ma2 - mb2)**2.Q+00*                &
             (k**2.Q+00 + mb2) + 6.Q+00*dnbdxb2*(k**2.Q+00 + ma2)**1.5Q+00*(ma2 - mb2)**3.Q+00*         &
             (k**2.Q+00 + mb2) + 36.Q+00*dnbdxb1*(k**2.Q+00 + ma2)**1.5Q+00*(ma2 - mb2)**2.Q+00*        &
             (k**2.Q+00 + mb2)**1.5Q+00 -                                                &
            dnbdxb3*(k**2.Q+00 + ma2)**1.5Q+00*(ma2 - mb2)**3.Q+00*(k**2.Q+00 + mb2)**1.5Q+00 +        &
            24.Q+00*(1.Q+00 + 2.Q+00*dnbdxb0)*(k**2.Q+00 + ma2)**1.5Q+00*(ma2 - mb2)*                   &
             (k**2.Q+00 + mb2)**2.Q+00 +                                                  &
            6.Q+00*(2.Q+00 + 4.Q+00*dnbdxb0)*(k**2.Q+00 + ma2)**1.5Q+00*(ma2 - mb2)*                    &
             (k**2.Q+00 + mb2)**2.Q+00 -                                                  &
            12.Q+00*dnbdxb2*(k**2.Q+00 + ma2)**1.5Q+00*(ma2 - mb2)**2.Q+00*(k**2.Q+00 + mb2)**2.Q+00 -       &
            72.Q+00*dnbdxb1*(k**2.Q+00 + ma2)**1.5Q+00*(ma2 - mb2)*(k**2.Q+00 + mb2)**2.5Q+00 -        &
            96.Q+00*(1.Q+00 + 2.Q+00*dnbdxb0)*(k**2.Q+00 + ma2)**1.5Q+00*(k**2.Q+00 + mb2)**3.Q+00 +              &
            96.Q+00*(1.Q+00 + 2.Q+00*dnbdxa0)*(k**2.Q+00 + ma2)*(k**2.Q+00 + mb2)**3.5Q+00 +                 &
            6.Q+00*(2.Q+00 + 4.Q+00*dnbdxa0)*(ma2 - mb2)*(k**2.Q+00 + mb2)**3.5Q+00 -                   &
            24.Q+00*dnbdxa1*Sqrt(k**2.Q+00 + ma2)*(ma2 - mb2)*(k**2.Q+00 + mb2)**3.5Q+00) +        &
         10.Q+00*(k**2.Q+00 + mb2)*(15.Q+00*(0.5Q+00 + dnbdxb0)*(k**2.Q+00 + ma2)**1.5Q+00*                 &
             (ma2 - mb2)**3.Q+00 - 15.Q+00*dnbdxb1*(k**2.Q+00 + ma2)**1.5Q+00*(ma2 - mb2)**3.Q+00*      &
             Sqrt(k**2.Q+00 + mb2) -                                                 &
            18.Q+00*(1.Q+00 + 2.Q+00*dnbdxb0)*(k**2.Q+00 + ma2)**1.5Q+00*(ma2 - mb2)**2.Q+00*                &
             (k**2.Q+00 + mb2) + 6.Q+00*dnbdxb2*(k**2.Q+00 + ma2)**1.5Q+00*(ma2 - mb2)**3.Q+00*         &
             (k**2.Q+00 + mb2) + 36.Q+00*dnbdxb1*(k**2.Q+00 + ma2)**1.5Q+00*(ma2 - mb2)**2.Q+00*        &
             (k**2.Q+00 + mb2)**1.5Q+00 -                                                &
            dnbdxb3*(k**2.Q+00 + ma2)**1.5Q+00*(ma2 - mb2)**3.Q+00*(k**2.Q+00 + mb2)**1.5Q+00 +        &
            24.Q+00*(1.Q+00 + 2.Q+00*dnbdxb0)*(k**2.Q+00 + ma2)**1.5Q+00*(ma2 - mb2)*                   &
             (k**2.Q+00 + mb2)**2.Q+00 +                                                  &
            6.Q+00*(2.Q+00 + 4.Q+00*dnbdxb0)*(k**2.Q+00 + ma2)**1.5Q+00*(ma2 - mb2)*                    &
             (k**2.Q+00 + mb2)**2.Q+00 -                                                  &
            12.Q+00*dnbdxb2*(k**2.Q+00 + ma2)**1.5Q+00*(ma2 - mb2)**2.Q+00*(k**2.Q+00 + mb2)**2.Q+00 -       &
            72.Q+00*dnbdxb1*(k**2.Q+00 + ma2)**1.5Q+00*(ma2 - mb2)*(k**2.Q+00 + mb2)**2.5Q+00 -        &
            96.Q+00*(1.Q+00 + 2.Q+00*dnbdxb0)*(k**2.Q+00 + ma2)**1.5Q+00*(k**2.Q+00 + mb2)**3.Q+00 +              &
            96.Q+00*(1.Q+00 + 2.Q+00*dnbdxa0)*(k**2.Q+00 + ma2)*(k**2.Q+00 + mb2)**3.5Q+00 +                 &
            6.Q+00*(2.Q+00 + 4.Q+00*dnbdxa0)*(ma2 - mb2)*(k**2.Q+00+ mb2)**3.5Q+00 -                   &
            24.Q+00*dnbdxa1*Sqrt(k**2.Q+00 + ma2)*(ma2 - mb2)*(k**2.Q+00 + mb2)**3.5Q+00)))/       &
    (384.Q+00*(k**2.Q+00 + ma2)**1.5Q+00*(ma2 - mb2)**6.Q+00*(k**2.Q+00 + mb2)**4.5Q+00)  

  end if

  BBall(i,1)=B2B2
  BBall(i,2)=B3B2
  BBall(i,3)=B4B2
  BBall(i,4)=B5B2
  BBall(i,5)=B2B3
  BBall(i,6)=B3B3
  BBall(i,7)=B4B3
  BBall(i,8)=B2B4
  BBall(i,9)=B3B4
  BBall(i,10)=B2B5

  end do

end

subroutine cal_F(k,mf2,dnfdx,dnadx,Fall)

  implicit none
  real(16) k 
  ! IR cutoff in flow equations
  integer Vmax
  !maximal number of potential
  parameter (Vmax=100)
  real(16) mf2(Vmax)
  real(16) nfa(Vmax),nfd1xa(Vmax),nfd2xa(Vmax),nfd3xa(Vmax)
  real(16) nfd4xa(Vmax),nfd5xa(Vmax)
  real(16) nff(Vmax),nfd1xf(Vmax),nfd2xf(Vmax),nfd3xf(Vmax)
  real(16) nfd4xf(Vmax),nfd5xf(Vmax)
  real(16) f2(Vmax),f3(Vmax),f4(Vmax),f5(Vmax),f6(Vmax)
  real(16) dnfdx(Vmax,6),dnadx(Vmax,6),Fall(Vmax,5)

  nff=dnfdx(:,1)
  nfd1xf=dnfdx(:,2)
  nfd2xf=dnfdx(:,3)
  nfd3xf=dnfdx(:,4)
  nfd4xf=dnfdx(:,5)
  nfd5xf=dnfdx(:,6)

  nfa=dnadx(:,1)
  nfd1xa=dnadx(:,2)
  nfd2xa=dnadx(:,3)
  nfd3xa=dnadx(:,4)
  nfd4xa=dnadx(:,5)
  nfd5xa=dnadx(:,6)

  f2=(k**3.Q+00*(1.Q+00 - nfa + Sqrt(k**2.Q+00 + mf2)*nfd1xa +                                           &
        Sqrt(k**2.Q+00 + mf2)*nfd1xf - nff))/(4.Q+00*(k**2.Q+00 + mf2)**1.5Q+00)
  f3=-(k**5.Q+00*(-3.Q+00 + 3.Q+00*nfa - 3.Q+00*Sqrt(k**2.Q+00 + mf2)*nfd1xa -                           &
         3.Q+00*Sqrt(k**2.Q+00 + mf2)*nfd1xf + k**2.Q+00*nfd2xa + mf2*nfd2xa +                           &
         k**2.Q+00*nfd2xf + mf2*nfd2xf + 3.Q+00*nff))/(16.Q+00*(k**2.Q+00 + mf2)**2.5Q+00)
  f4=(k**7.Q+00*(15.Q+00 - 15.Q+00*nfa + 15.Q+00*Sqrt(k**2.Q+00 + mf2)*nfd1xa +                          &
        15.Q+00*Sqrt(k**2.Q+00 + mf2)*nfd1xf - 6.Q+00*k**2.Q+00*nfd2xa - 6.Q+00*mf2*nfd2xa -             &
        6.Q+00*k**2.Q+00*nfd2xf - 6.Q+00*mf2*nfd2xf + k**2.Q+00*Sqrt(k**2.Q+00 + mf2)*nfd3xa +           &
        mf2*Sqrt(k**2.Q+00 + mf2)*nfd3xa + k**2.Q+00*Sqrt(k**2.Q+00 + mf2)*nfd3xf +                      &
        mf2*Sqrt(k**2.Q+00 + mf2)*nfd3xf - 15.Q+00*nff))/(96.Q+00*(k**2.Q+00 + mf2)**3.5Q+00)
  f5=-(k**9.Q+00*(-105.Q+00 + 105.Q+00*nfa - 105.Q+00*Sqrt(k**2.Q+00 + mf2)*nfd1xa -                     &
         105.Q+00*Sqrt(k**2.Q+00 + mf2)*nfd1xf + 45.Q+00*k**2.Q+00*nfd2xa + 45.Q+00*mf2*nfd2xa +         &
         45.Q+00*k**2.Q+00*nfd2xf + 45.Q+00*mf2*nfd2xf -                                                 &
         10.Q+00*k**2.Q+00*Sqrt(k**2.Q+00 + mf2)*nfd3xa -                                                &
         10.Q+00*mf2*Sqrt(k**2.Q+00 + mf2)*nfd3xa -                                                      &
         10.Q+00*k**2.Q+00*Sqrt(k**2.Q+00 + mf2)*nfd3xf -                                                &
         10.Q+00*mf2*Sqrt(k**2.Q+00 + mf2)*nfd3xf + k**4.Q+00*nfd4xa +                                   &
         2.Q+00*k**2.Q+00*mf2*nfd4xa + mf2**2.Q+00*nfd4xa + k**4.Q+00*nfd4xf +                           &
         2.Q+00*k**2.Q+00*mf2*nfd4xf + mf2**2.Q+00*nfd4xf + 105.Q+00*nff))/                              &
    (768.Q+00*(k**2.Q+00 + mf2)**4.5Q+00)
  f6=  (k**11.Q+00*(945.Q+00 - 945.Q+00*nfa + 945.Q+00*Sqrt(k**2.Q+00 + mf2)*nfd1xa +                    &
        945.Q+00*Sqrt(k**2.Q+00 + mf2)*nfd1xf - 420.Q+00*k**2.Q+00*nfd2xa - 420.Q+00*mf2*nfd2xa -        &
        420.Q+00*k**2.Q+00*nfd2xf - 420.Q+00*mf2*nfd2xf +                                                &
        105.Q+00*k**2.Q+00*Sqrt(k**2.Q+00 + mf2)*nfd3xa +                                                &
        105.Q+00*mf2*Sqrt(k**2.Q+00 + mf2)*nfd3xa +                                                      &
        105.Q+00*k**2.Q+00*Sqrt(k**2.Q+00 + mf2)*nfd3xf +                                                &
        105.Q+00*mf2*Sqrt(k**2.Q+00 + mf2)*nfd3xf - 15.Q+00*k**4.Q+00*nfd4xa -                           &
        30.Q+00*k**2.Q+00*mf2*nfd4xa - 15.Q+00*mf2**2.Q+00*nfd4xa - 15.Q+00*k**4.Q+00*nfd4xf -           &
        30.Q+00*k**2.Q+00*mf2*nfd4xf - 15.Q+00*mf2**2.Q+00*nfd4xf +                                      &
        k**4.Q+00*Sqrt(k**2.Q+00 + mf2)*nfd5xa + 2.Q+00*k**2.Q+00*mf2*Sqrt(k**2.Q+00 + mf2)*nfd5xa +     &
        mf2**2.Q+00*Sqrt(k**2.Q+00 + mf2)*nfd5xa + k**4.Q+00*Sqrt(k**2.Q+00 + mf2)*nfd5xf +              &
        2.Q+00*k**2.Q+00*mf2*Sqrt(k**2.Q+00 + mf2)*nfd5xf +                                              &
        mf2**2.Q+00*Sqrt(k**2.Q+00 + mf2)*nfd5xf - 945.Q+00*nff))/(7680.Q+00*(k**2.Q+00 + mf2)**5.5Q+00)

 Fall(:,1)=f2
 Fall(:,2)=f3
 Fall(:,3)=f4
 Fall(:,4)=f5
 Fall(:,5)=f6

end
