subroutine lang (density)

!use ifport

implicit none

double precision::density

integer::k,l

integer::ia,ip,iq,ik
double precision::pr,qr,kr

double precision::a,b

integer,parameter::t_size=100
integer,parameter::a_size=299

double precision::factv1,factv2
double precision,dimension(0:a_size)::factd1
double precision,parameter::pi=4.0d0*datan(1.0d0)
double precision::conv=1.0d0 !Convergence of the calculation
double precision,parameter::prec=1.0d-12 !Precision of the calculation

double precision::hk
double precision::hq
double precision::mult

integer::kk
double precision::iter

double precision,dimension(0:a_size)::q
double precision,dimension(0:a_size)::ccq
double precision,dimension(0:a_size)::hcq
double precision,dimension(0:a_size)::hdq
double precision,dimension(0:a_size)::scq
double precision,dimension(0:a_size)::sdq

double precision,dimension(0:a_size)::tau=0.0d0

double precision,dimension(1:2*t_size,0:a_size)::phi=0.0d0
double precision,dimension(1:2*t_size,0:a_size)::dphi=0.0d0

!For the convergence of phi and m in the subiteration
double precision,dimension(0:a_size)::phik=0.0d0
double precision,dimension(0:a_size)::phik2=0.0d0
double precision::conv_phik=1.0d0
double precision::prec_phik=1.0d-12
double precision,dimension(0:a_size)::mk=0.0d0

double precision,dimension(1:2*t_size,0:a_size)::m=0.0d0
double precision,dimension(1:2*t_size,0:a_size)::dm=0.0d0

double precision,dimension(0:a_size)::ck=0.0d0
double precision,dimension(0:a_size)::dk=0.0d0

double precision,dimension(0:a_size,0:a_size,0:a_size)::v1,v2

integer::wits
!logical::sys=.false.


character(len=17)::phi_file
character(len=15)::m_file

call fileman('ccq.dat',7,11,1)
call fileman('hcq.dat',7,12,1)
call fileman('hdq.dat',7,13,1)
do iq=0,a_size
  read(11,*) a,b
  read(11,*) q(iq),ccq(iq)
  read(12,*) a,b
  read(12,*) q(iq),hcq(iq)
  read(13,*) a,b
  read(13,*) q(iq),hdq(iq)
end do
call fileman('ccq.dat',7,11,0)
call fileman('hcq.dat',7,12,0)
call fileman('hdq.dat',7,13,0)

do ip=0,a_size
  scq(ip)=1.0d0+density*hcq(ip)
  sdq(ip)=density*hdq(ip)
  tau(ip)=(160.0d0*scq(ip))/(q(ip)**2)
end do

hk=1.0d-4*tau(wits(tau,a_size))

hq=q(2)-q(1) !h, the pitch
factv2=density*(hq**3)/(32.0d0*(pi**2))
factv1=(hq**3)/(16.0d0*(pi**2))

do ip=0,a_size
  pr=dble(ip)+0.5d0
  do ik=0,a_size
    kr=dble(ik)+0.5d0
    do iq=0,a_size
      qr=dble(iq)+0.5d0

      if ((pr .ge. (dabs(qr-kr)+0.5d0)) .and. (pr .le. (qr+kr-0.5d0))) then

        v2(iq,ik,ip)=factv2*scq(iq)*scq(ik)*scq(ip)*&
        (kr*pr/(qr**5))*(((kr**2)+(qr**2)-(pr**2))*ccq(ik)+((pr**2)+(qr**2)-(kr**2))*ccq(ip))**2

        v1(iq,ik,ip)=factv1*scq(iq)*scq(ik)*hdq(ip)*&
        (kr*pr/(qr**5))*(((kr**2)+(qr**2)-(pr**2))*ccq(ik)*density+((pr**2)+(qr**2)-(kr**2)))**2

      end if

    end do
  end do
end do





!First calculation of phi and its derivative
!=================================================
do ia=0,a_size
  do k=1,t_size
    phi(k,ia)=1.0d0-(hk*dble(k))/tau(ia)
  end do
end do

do ia=0,a_size
  dphi(1,ia)=0.5d0*(1.0d0+phi(1,ia))
  do k=2,t_size
    dphi(k,ia)=0.5d0*(phi(k-1,ia)+phi(k,ia))*(hk*dble(k)-hk*dble(k-1))
  end do
end do

!First calculation of m and its derivative
!=================================================
do ip=0,a_size
  do ik=0,a_size
    do iq=0,a_size
      do k=0,t_size
        m(k,iq)=m(k,iq)+v2(iq,ik,ip)*phi(k,ip)*phi(k,ik)+v1(iq,ik,ip)*phi(k,ik)
      end do
    end do
  end do
end do

do ia=0,a_size
  dm(1,ia)=0.5d0*(1.0d0+m(1,ia))
  do k=2,t_size
    dm(k,ia)=0.5d0*(m(k-1,ia)+m(k,ia))*(hk*dble(k)-hk*dble(k-1))
  end do
end do

!Create a folder to stock all the phi files
!sys=systemqq('mkdir -p lang/')
call system ('mkdir lang/')


write (6,*) ""
write (6,'(a47)') "          time           ||         convergence"

!Write the first values
do ia=0,a_size
  write (phi_file,'(a10,i3.3,a4)') './lang/phi',ia,'.dat'
  write (m_file,'(a8,i3.3,a4)') './lang/m',ia,'.dat'
  call fileman(phi_file,len(phi_file),10+ia,1)
  call fileman(m_file,len(m_file),11+ia,1)
  write (10+ia,*) hk*1.0d0,phi(1,ia)
  write (11+ia,*) hk*1.0d0,m(1,ia)
  do k=2,t_size
    write (10+ia,*) hk*dble(k),phi(k,ia),dphi(k,ia)
    write (11+ia,*) hk*dble(k),m(k,ia),dm(k,ia)
  end do
  call fileman(phi_file,len(phi_file),10+ia,0)
  call fileman(m_file,len(m_file),11+ia,0)
end do

k=t_size+1
iter=dble(t_size)+1.0d0
mult=1

!ITERATION ON THE CONVERGENCE OF PHI(K)
do while (conv .gt. prec)

!  write (6,*) "mult",mult

  !ITERATION INCREMENTING K
  do while (k .le. t_size*2)

    kk=aint(dble(k)/2.0d0)

!    write (6,*) "iter=",iter
!    write (6,*) "k=",k
!    write (6,*) "kk=",kk

    !Calculate C
    !=================================================
    do ia=0,a_size
      ck(ia)=m(kk,ia)*phi(k-kk,ia)-m(k-1,ia)*dphi(1,ia)-phi(k-1,ia)*dm(1,ia)
    end do
!    write (6,*) "ck=",ck

    do l=2,kk
      do ia=0,a_size
        ck(ia)=ck(ia)+(phi(k-l+1,ia)-phi(k-l,ia))*dm(l,ia)
      end do
    end do
!    write (6,*) "ck=",ck

    do l=2,k-kk
      do ia=0,a_size
        ck(ia)=ck(ia)+(m(k-l+1,ia)-m(k-l,ia))*dphi(l,ia)
      end do
    end do
!    write (6,*) "ck=",ck

    !Calculate D
    !=================================================
    do ia=0,a_size
      factd1(ia)=tau(ia)/(2.0d0*hk*mult)
      dk(ia)=ck(ia)-factd1(ia)*(4.0d0*phi(k-1,ia)-phi(k-2,ia))
    end do

!    write (6,*) "dk=",dk
    do ia=0,a_size
      phik(ia)=phi(k-1,ia)
    end do

    !ITERATE TO CONVERGE PHI AND M
    conv_phik=1.0d0
    do while (conv_phik .gt. prec_phik)

      do ip=0,a_size
        do ik=0,a_size
          do iq=0,a_size
            mk(iq)=mk(iq)+v2(iq,ik,ip)*phik(ip)*phik(ik)+v1(iq,ik,ip)*phik(ik)
          end do
        end do
      end do

      do ia=0,a_size
        phik2(ia)=(mk(ia)*(1.0d0-dphi(1,ia))-dk(ia))/(1.0d0+dm(1,ia)+3.0d0*factd1(ia))
      end do

      conv_phik=0.0d0
      do ia=0,a_size
        conv_phik=conv_phik+dabs(phik2(ia)-phik(ia))
        mk(ia)=0.0d0
        phik(ia)=phik2(ia)
      end do

!      write (6,*) "conv_phik = ",conv_phik

    end do

    do ia=0,a_size
      phi(k,ia)=phik(ia)
      dphi(k,ia)=0.5d0*(phi(k-1,ia)+phi(k,ia))*(hk*mult)
    end do

    do ia=0,a_size
      m(k,ia)=0.0d0
    end do

    do ip=0,a_size
      do ik=0,a_size
        do iq=0,a_size
          m(k,iq)=m(k,iq)+v2(iq,ik,ip)*phi(k,ip)*phi(k,ik)+v1(iq,ik,ip)*phi(k,ik)
        end do
      end do
    end do

    do ia=0,a_size
      dm(k,ia)=0.5d0*(m(k-1,ia)+m(k,ia))*(hk*mult)
    end do


    do ia=0,a_size
      write (phi_file,'(a10,i3.3,a4)') './lang/phi',ia,'.dat'
      write (m_file,'(a8,i3.3,a4)') './lang/m',ia,'.dat'
      call fileman(phi_file,len(phi_file),10+ia,2)
      call fileman(m_file,len(m_file),11+ia,2)
      write (10+ia,*) hk*dble(iter),phi(k,ia),dphi(k,ia)
      write (11+ia,*) hk*dble(iter),m(k,ia),dm(k,ia)
      call fileman(phi_file,len(phi_file),10+ia,0)
    end do

!    write (11,*) hk*iter,phi(k)
!    write (12,*) hk*iter,m(k)

    conv=0.0d0
    do ia=0,a_size
      conv=conv+abs(phi(k,ia)-phi(k-1,ia))
    end do

    write (6,'(es24.16,a4,es24.16)') hk*dble(iter)," || ",conv

    k=k+1
    iter=iter+mult

  end do

  !Compress the calculation space
  do ia=0,a_size
    do k=1,t_size
      phi(k,ia)=phi(2*k,ia)
      m(k,ia)=m(2*k,ia)
    end do
  end do

  !Calculate the derivatives
  do ia=0,a_size
    do k=1,t_size/2
      dphi(k,ia)=0.5d0*(dphi(2*k-1,ia)+dphi(2*k,ia))
      dm(k,ia)=0.5d0*(dm(2*k-1,ia)+dm(2*k,ia))
    end do
  end do

  do ia=0,a_size
    do k=t_size/2+1,t_size
      dphi(k,ia)=(phi(k,ia)+4.0d0*phi(k-1,ia)+phi(k-2,ia))/6.0d0
      dm(k,ia)=(m(k,ia)+4.0d0*m(k-1,ia)+m(k-2,ia))/6.0d0
    end do
  end do

  k=t_size+1
  mult=2.0d0*mult

end do

do ia=0,a_size
  close(10+ia)
end do

!close(11)
!close(12)

end subroutine
