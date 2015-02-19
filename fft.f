      subroutine fft3s(f0,f,mm,rl,t,mode,mode1)

!     f0=0.0
!     f la fonction à transformer (tableau)
!     mm la puissance à laquelle on met 2, la dimension du tableau (1024==>mm=10)
!     rl la distance maximale dans le tableau des distances
!     t la transforlée de la fonction
!     mode=1 pour la TF
!          -1 pour la TF-1
!     mode1=0
!
!     FFT of a 3-dimensional centrosymmetric function
!
         implicit real*8(a-h,o-z)
         dimension f(1),t(1)
         double precision f0
         common /fdata/ n,n2,ff,rttwo,nc,ns,nd,nr
         common /fwork/ w(50000)
!
!     dimension of w should be .ge. 3*n
!     fixed data preparation
!
         data pi/3.1415926535897931d0/,m/-1/
!
      md = mm+1
!
      call setft(md,m)
!
      n = 2**mm
      n1=n+1
      rn = dfloat(n)
      dr = rl/rn
      rl3 = rl*rl*rl
      c = 4.d0*rl3*ff/rn
!
      if(mode.eq.-1) c=1.0d0/c
      if(mode1.eq.1) c=c*2.0d0*dr*pi*pi
      if(mode1.eq.2) c=c*2.0d0*pi*pi*dr*rl/pi
!
!     copy input vector into working storage
!
      do 204 i=1,n1
         w(nd+i) = c*ff*f(i)
         if(mode1.ne.2) w(nd+i)=w(nd+i)*(i-1)
 204  continue
!
      if (mode1.eq.2) go to 1000
!
      call rsa(m)
!
      w(nr+n1)=0.0d0
      go to 2000
!
1000  call rca(m)
!
!     copy into output vector
!
2000  do 205 i=2,n1
         t(i) = w(nr+i)/dfloat(i-1)
         if(mode1.ne.0) t(i)=w(nr+i)
 205  continue
!
      return
      end
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine setft(mm,m)
!
!     prepare fixed data for fast fourier transform
!
         implicit real*8(a-h,o-z)
         common /fdata/ n,n2,f,rttwo,nc,ns,ndata,nres
         common /fwork/ w(50000)
!
      if(mm.lt.0)return
      if(mm.eq.m)return
!
      m = mm
      n = 2**m
      n2 = n/2
      f = 1.0d0/dfloat(n)
      da = 6.2831853071796d0*f
      rttwo = dsqrt(2.0d0)
      f = dsqrt(f)
      n1 = (n2-1)/2
      nc = 0
      ns = n1
      ndata = n1+n1
      if(n1.le.0)return
!
      do 10 i=1,n1
         a = dfloat(i)*da
         w(i) = dcos(a)
 10      w(n1+i) = dsin(a)
!
      return
      end
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine rsa(m)
!
!     real odd analysis or synthesis
!
         implicit real*8(a-h,o-z)
         common /fdata/ n,n2,f,rttwo,nc,ns,ndata,nres
         common /fwork/ w(50000)
!
      nq = n2
      nib = ndata +1
      nob = ndata+n2+2
!
      do 200 l = 1,m
         nim = nib+nq
         nie = nim+nq
         nqh = nq/2
         nql = (nq-1)/2
         nom = nob+nql+1
!
!     tprime = 0
!
         if(nql.eq.0)goto 20
!
         do 10 it = 1,nql
            w(nob+it)=w(nib+it)-w(nim-it)
 10         w(nom+it) = w(nib+it)+w(nim-it)
!
 20      if(nqh.ne.nql)w(nom) = 2.0d0*w(nib+nqh)
!
!     0.lt.tprime.lt.np/2.
!
         no1 = nq
         no2 = n2-nq
         if(no1-no2)40,80,120
 40      ni=no1+no1
         w(nob+no1) = w(nib+ni)+w(nim+ni)
         w(nob+no2) = -w(nib+ni)+w(nim+ni)
         if(nql.eq.0)goto 60
         cc = w(nc+no1)
         ss = w(ns+no1)
!
         do 50 it = 1,nql
            re = cc*w(nie+ni-it)-ss*w(nim+ni-it)
            ai = ss*w(nie+ni-it)+cc*w(nim+ni-it)
            w(nom+no1+it) = w(nim+ni+it)-re
            w(nom+no2+it) = w(nim+ni+it)+re
            w(nob+no1+it) = w(nib+ni+it)+ai
 50         w(nob+no2+it) = -w(nib+ni+it)+ai
!
 60      if(nqh.eq.nql)goto 70
         nr = no1/2
         w(nom+no1) = 2.0d0*(w(ns+nr)*w(nim+ni+nqh)+w(nc+nr)*
     &        w(nib+ni+nqh))
         w(nom+no2) = 2.0d0*(w(nc+nr)*w(nim+ni+nqh)-w(ns+nr)*
     &        w(nib+ni+nqh))
 70      no1 = no1+nq
         no2 = no2-nq
         if(no1-no2)40,80,120
 80      w(nob+no1)=w(nim)
         if(nql.eq.0)goto 100
!
         do 90 it = 1,nql
            w(nom+no1+it) = w(nim+it)
 90         w(nob+no1+it) = w(nie-it)
!
 100     if(nqh.ne.nql)w(nom+no1) = rttwo*w(nim+nqh)
 120     nt = nib
         nib = nob
         nob = nt
         nq = nqh
 200  continue
!
      nres = nib-1
      return
      end
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      subroutine rca(m)
!
!     real even analysis or synthesis.
!
         implicit real*8(a-h,o-z)
         common/fdata/ n,n2,f,rttwo,nc,ns,ndata,nres
         common/fwork/ w(50000)
!
      nq=n2
      nib=ndata+1
      nob=ndata+n2+2
!
      do 200 l=1,m
         nim=nib+nq
         nie=nim+nq
         nqh=nq/2
         nql=(nq-1)/2
         nom=nob+nql+1
!
!     tprime=0
!
         w(nob)=w(nib)+w(nim)
         w(nob+n2)=w(nib)-w(nim)
!
         if(nql.eq.0)goto 20
!
         do 10 it = 1,nql
            w(nob+it)=w(nib+it)+w(nim-it)
 10         w(nom+it) = w(nib+it)-w(nim-it)
!
 20      if(nqh.ne.nql)w(nom) = 2.0d0*w(nib+nqh)
!
!     0.lt.tprime.lt.np/2.
!
         no1 = nq
         no2 = n2-nq
         if(no1-no2)40,80,120
 40      ni=no1+no1
         w(nob+no1) = w(nib+ni)+w(nim+ni)
         w(nob+no2) = w(nib+ni)-w(nim+ni)
         if(nql.eq.0)goto 60
         cc = w(nc+no1)
         ss = w(ns+no1)
!
         do 50 it = 1,nql
            re = cc*w(nim+ni-it)-ss*w(nie+ni-it)
            ai = ss*w(nim+ni-it)+cc*w(nie+ni-it)
            w(nob+no1+it) = w(nib+ni+it)+re
            w(nob+no2+it) = w(nib+ni+it)-re
            w(nom+no1+it) = w(nim+ni+it)-ai
 50         w(nom+no2+it) = -w(nim+ni+it)-ai
!
 60      if(nqh.eq.nql)goto 70
         nr = no1/2
         w(nom+no1)=2.d0*(w(nc+nr)*w(nib+ni+nqh)-w(ns+nr)*w(nim+ni+nqh))
         w(nom+no2)=2.d0*(w(ns+nr)*w(nib+ni+nqh)+w(nc+nr)*w(nim+ni+nqh))
 70      no1 = no1+nq
         no2 = no2-nq
         if(no1-no2)40,80,120
!
!     tprime=np/2
!
 80      w(nob+no1)=w(nib+n2)
         if(nql.eq.0)goto 100
!
         do 90 it = 1,nql
            w(nob+no1+it) = w(nim+it)
 90         w(nom+no1+it) = -w(nie-it)
!
 100     if(nqh.ne.nql)w(nom+no1) = rttwo*w(nim+nqh)
 120     nt = nib
         nib = nob
         nob = nt
         nq = nqh
 200  continue
!
      nres = nib-1
      return
      end
! 
