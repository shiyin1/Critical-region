      FUNCTION rtbis(func,x1,x2,xacc)
      INTEGER JMAX
      REAL(16) rtbis,x1,x2,xacc,func
!      EXTERNAL func
      PARAMETER (JMAX=80)
      INTEGER j
      REAL(16) dx,f,fmid,xmid

      call gapEq(x2,fmid)
      call gapEq(x1,f)
!      fmid=func(x2)
!      f=func(x1)
      if(f*fmid.ge.0.Q+00) pause 'root must be bracketed in rtbis'
      if(f.lt.0.Q+00)then
        rtbis=x1
        dx=x2-x1
      else
        rtbis=x2
        dx=x1-x2
      endif
      do 11 j=1,JMAX
        dx=dx*0.5Q+00
        xmid=rtbis+dx
        call gapEq(xmid,fmid)
!        fmid=func(xmid)
        if(fmid.le.0.Q+00)rtbis=xmid
        if(abs(dx).lt.xacc .or. fmid.eq.0.Q+00) return
11    continue
      pause 'too many bisections in rtbis'
      END
