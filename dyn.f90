subroutine dyn (density)

implicit none

double precision::density

integer::k,l

integer::ia,ip,iq,ik
double precision::pr,qr,kr

double precision::a,b

integer,parameter::t_size=100
integer,parameter::a_size=299

double precision::factv1,factv2,factw1,factw2
double precision,dimension(0:a_size)::factd1
double precision,dimension(0:a_size)::factd1_s
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

double precision,parameter::tmic=160.0d0
double precision,dimension(0:a_size)::tau=0.0d0
double precision,dimension(0:a_size)::tau_s=0.0d0

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
double precision::m0=0.0d0

double precision,dimension(0:a_size)::ck=0.0d0
double precision,dimension(0:a_size)::dk=0.0d0

double precision,dimension(0:a_size,0:a_size,0:a_size)::v1,v2

!For the calculation of the self dynamics
double precision,dimension(1:2*t_size,0:a_size)::phi_s=0.0d0
double precision,dimension(1:2*t_size,0:a_size)::dphi_s=0.0d0

double precision,dimension(0:a_size)::phik_s=0.0d0
double precision,dimension(0:a_size)::phik_s2=0.0d0

double precision,dimension(0:a_size)::mk_s=0.0d0

double precision,dimension(1:2*t_size,0:a_size)::m_s=0.0d0
double precision,dimension(1:2*t_size,0:a_size)::dm_s=0.0d0
double precision::m_s0=0.0d0

double precision,dimension(0:a_size)::ck_s=0.0d0
double precision,dimension(0:a_size)::dk_s=0.0d0

double precision,dimension(0:a_size,0:a_size,0:a_size)::v1_s,v2_s

!For the calculation of the mean square displacement
double precision,dimension(0:a_size)::w1,w2
double precision::d0=1.0d0/tmic

double precision::ck_msd=0.0d0
double precision::dk_msd=0.0d0

double precision,dimension(1:2*t_size)::m_msd=0.0d0
double precision,dimension(1:2*t_size)::dm_msd=0.0d0
double precision::m_msd0=0.0d0

double precision,dimension(1:2*t_size)::dr2=0.0d0
double precision,dimension(1:2*t_size)::ddr2=0.0d0



integer::wits

character(len=20)::phi_file
character(len=20)::phi_s_file
character(len=20)::dr2_file

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

do iq=0,a_size
  scq(iq)=1.0d0+density*hcq(iq)
  sdq(iq)=density*hdq(iq)

  tau(iq)=(tmic*scq(iq))/(q(iq)**2)
  tau_s(iq)=tmic/(q(iq)**2)
end do

hk=1.0d-4*minval(tau,dim=1)/dble(t_size)

hq=q(2)-q(1) !h, the pitch

factv2=density*(hq**3)/(32.0d0*(pi**2))
factv1=(hq**3)/(16.0d0*(pi**2))

factw2=(hq*density)/(12.0d0*pi**2)
factw1=hq/(12.0d0*pi**2)

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

        v2_s(iq,ik,ip)=2.0d0*factv2*(ccq(ip)**2)*scq(ip)*(kr*pr/(qr**5))*((pr**2)+(qr**2)-(kr**2))**2

        v1_s(iq,ik,ip)=factv1*hdq(ip)*(kr*pr/(qr**5))*((pr**2)+(qr**2)-(kr**2))**2

      end if

    end do
  end do

  w2(ip)=factw2*(q(ip)**4)*(ccq(ip)**2)*scq(ip)

  w1(ip)=factw1*(q(ip)**4)*hdq(ip)

end do


!First calculation of phi and its derivative
!=================================================
do ia=0,a_size
  do k=1,t_size
    phi(k,ia)=1.0d0-(hk*dble(k))/tau(ia)
    phi_s(k,ia)=1.0d0-(hk*dble(k))/tau_s(ia)
  end do
end do

!Calculation of the mean square displacement
do k=1,t_size
  dr2(k)=6.0d0*d0*hk*dble(k)
end do

do ia=0,a_size
  dphi(1,ia)=0.5d0*(1.0d0+phi(1,ia))
  dphi_s(1,ia)=0.5d0*(1.0d0+phi_s(1,ia))
  do k=2,t_size
    dphi(k,ia)=0.5d0*(phi(k-1,ia)+phi(k,ia))
    dphi_s(k,ia)=0.5d0*(phi_s(k-1,ia)+phi_s(k,ia))
  end do
end do

ddr2(1)=0.5d0*dr2(1)
do k=2,t_size
  ddr2(k)=0.5d0*(dr2(k-1)+dr2(k))
end do

!First calculation of m and its derivative
!=================================================
do ip=0,a_size
  do ik=0,a_size
    do iq=0,a_size
      m0=m0+v2(iq,ik,ip)+v1(iq,ik,ip)
      m_s0=m_s0+v2_s(iq,ik,ip)+v1_s(iq,ik,ip)
      do k=1,t_size
        m(k,iq)=m(k,iq)+v2(iq,ik,ip)*phi(k,ip)*phi(k,ik)+v1(iq,ik,ip)*phi(k,ik)

        m_s(k,iq)=m_s(k,iq)+v2_s(iq,ik,ip)*phi_s(k,ik)*phi(k,ip)+v1_s(iq,ik,ip)*phi_s(k,ik)
      end do
    end do
  end do
end do

do ia=0,a_size
  m_msd0=m_msd0+w2(ia)+w1(ia)
  do k=1,t_size
    m_msd(k)=m_msd(k)+w2(ia)*phi_s(k,ia)*phi(k,ia)+w1(ia)*phi_s(k,ia)
  end do
end do

do ia=0,a_size
  dm(1,ia)=0.5d0*(m0+m(1,ia))
  dm_s(1,ia)=0.5d0*(m_s0+m_s(1,ia))
  dm_msd(1)=0.5d0*(m_msd0+m_msd(1))
  do k=2,t_size
    dm(k,ia)=0.5d0*(m(k-1,ia)+m(k,ia))
    dm_s(k,ia)=0.5d0*(m_s(k-1,ia)+m_s(k,ia))
    dm_msd(k)=0.5d0*(m_msd(k-1)+m_msd(k))
  end do
end do

!Create a folder to stock all the phi files
call system ('mkdir dyn/')

write (6,'(a47)') "          time           ||         convergence"

!Write the first values
do ia=0,a_size
  write (phi_file,'(a10,i3.3,a4)') './dyn/phi_',ia,'.dat'
  write (phi_s_file,'(a12,i3.3,a4)') './dyn/phi_s_',ia,'.dat'
  write (dr2_file,'(a13)') './dyn/dr2.dat'

  call fileman(phi_file,len(phi_file),10+ia,1)
  call fileman(phi_s_file,len(phi_s_file),310+ia,1)
  call fileman(dr2_file,len(dr2_file),610,1)

   write (10+ia,*) 0.0d0,1.0d0
   write (310+ia,*) 0.0d0,1.0d0
   write (610+ia,*) 0.0d0,0.0d0
  do k=1,t_size
    write (10+ia,*) hk*dble(k),phi(k,ia),dphi(k,ia)
    write (310+ia,*) hk*dble(k),phi_s(k,ia),dphi_s(k,ia)
    write (610,*) hk*dble(k),dr2(k)
  end do

  call fileman(phi_file,len(phi_file),10+ia,0)
  call fileman(phi_s_file,len(phi_s_file),310+ia,0)
  call fileman(dr2_file,len(dr2_file),610,0)
end do

k=t_size
iter=dble(k)
mult=1.0d0

!ITERATION ON THE CONVERGENCE OF PHI(K)
do while (conv .gt. prec)

  !ITERATION INCREMENTING K
  do while (k .lt. t_size*2)

    k=k+1
    kk=aint(dble(k)/2.0d0)
    iter=iter+mult

    !Calculate C
    !=================================================
    do ia=0,a_size
      ck(ia)=m(kk,ia)*phi(k-kk,ia)-m(k-1,ia)*dphi(1,ia)-phi(k-1,ia)*dm(1,ia)
      ck_s(ia)=m_s(kk,ia)*phi_s(k-kk,ia)-m_s(k-1,ia)*dphi_s(1,ia)-phi_s(k-1,ia)*dm_s(1,ia)
    end do
    do ia=0,a_size
      do l=2,kk
        ck(ia)=ck(ia)+(phi(k-l+1,ia)-phi(k-l,ia))*dm(l,ia)
        ck_s(ia)=ck_s(ia)+(phi_s(k-l+1,ia)-phi_s(k-l,ia))*dm_s(l,ia)
      end do
    end do
    do ia=0,a_size
      do l=2,k-kk
        ck(ia)=ck(ia)+(m(k-l+1,ia)-m(k-l,ia))*dphi(l,ia)
        ck_s(ia)=ck_s(ia)+(m_s(k-l+1,ia)-m_s(k-l,ia))*dphi_s(l,ia)
      end do
    end do

    ck_msd=m_msd(kk)*dr2(k-kk)-m_msd(k-1)*ddr2(1)-dr2(k-1)*dm_msd(1)
    do l=2,kk
      ck_msd=ck_msd+(dr2(k-l+1)-dr2(k-l))*dm_msd(l)
    end do
    do l=2,k-kk
      ck_msd=ck_msd+(m_msd(k-l+1)-m_msd(k-l))*ddr2(l)
    end do

    !Calculate D
    !=================================================
    do ia=0,a_size
      factd1(ia)=tau(ia)/(2.0d0*hk*mult)
      factd1_s(ia)=tau_s(ia)/(2.0d0*hk*mult)
      dk(ia)=ck(ia)-factd1(ia)*(4.0d0*phi(k-1,ia)-phi(k-2,ia))
      dk_s(ia)=ck_s(ia)-factd1_s(ia)*(4.0d0*phi_s(k-1,ia)-phi_s(k-2,ia))

      !take the last value of phi for a good guess
      phik(ia)=phi(k-1,ia)
      phik_s(ia)=phi_s(k-1,ia)
    end do

    !ITERATE TO CONVERGE PHI AND 
    !=================================================
    conv_phik=1.0d0
    do while (conv_phik .gt. prec_phik)

      do ia=0,a_size
        mk(ia)=0.0d0
        mk_s(ia)=0.0d0
      end do

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

      do ip=0,a_size
        do ik=0,a_size
          do iq=0,a_size
            mk_s(iq)=mk_s(iq)+v2_s(iq,ik,ip)*phik_s(ik)*phik2(ip)+v1_s(iq,ik,ip)*phik_s(ik)
          end do
        end do
      end do

      do ia=0,a_size
        phik_s2(ia)=(mk_s(ia)*(1.0d0-dphi_s(1,ia))-dk_s(ia))/(1.0d0+dm_s(1,ia)+3.0d0*factd1_s(ia))
      end do

      conv_phik=0.0d0
      do ia=0,a_size
        conv_phik=conv_phik+dabs(phik2(ia)-phik(ia))+dabs(phik_s2(ia)-phik_s(ia))

        phik(ia)=phik2(ia)
        phik_s(ia)=phik_s2(ia)
      end do

    end do

    !One last calculation of phi, m and their primitives
    do ia=0,a_size
      phi(k,ia)=phik(ia)
      phi_s(k,ia)=phik_s(ia)

      m(k,ia)=mk(ia)
      m_s(k,ia)=mk_s(ia)
    end do

    !Calculation of the mean square displacement
    m_msd(k)=0.0d0
    do ia=0,a_size
      m_msd(k)=m_msd(k)+w2(ia)*phi_s(k,ia)*phi(k,ia)+w1(ia)*phi_s(k,ia)
    end do

    dk_msd=m_msd(k)*ddr2(1)-(4.0d0*dr2(k-1)-dr2(k-2))/(2.0d0*d0*hk*mult)+ck_msd

    dr2(k)=(6.0d0*d0-d0*dk_msd)/(d0*dm_msd(1)+3.0d0/(2.0d0*hk*mult))

    !Write the stuff
    !=================================================
    do ia=0,a_size
      write (phi_file,'(a10,i3.3,a4)') './dyn/phi_',ia,'.dat'
      write (phi_s_file,'(a12,i3.3,a4)') './dyn/phi_s_',ia,'.dat'

      call fileman(phi_file,len(phi_file),10+ia,2)
      call fileman(phi_s_file,len(phi_s_file),310+ia,2)

      write (10+ia,*) hk*iter,phi(k,ia),dphi(k,ia)
      write (310+ia,*) hk*iter,phi_s(k,ia),dphi_s(k,ia)

      call fileman(phi_file,len(phi_file),10+ia,0)
      call fileman(phi_s_file,len(phi_s_file),310+ia,0)
    end do

    write (dr2_file,'(a13)') './dyn/dr2.dat'
    call fileman(dr2_file,len(dr2_file),610,2)
    write (610,*) hk*iter,dr2(k)
    call fileman(dr2_file,len(dr2_file),610,0)

    conv=0.0d0
    do ia=0,a_size
      conv=conv+abs(phi(k,ia)-phi(k-1,ia))+abs(phi_s(k,ia)-phi_s(k-1,ia))
    end do

    write (6,'(es24.16,a4,es24.16)') hk*iter," || ",conv

  end do

  !Calculate the integrals
  !=================================================
  do ia=0,a_size
    do k=1,t_size/2
      dphi(k,ia)=0.5d0*(dphi(2*k-1,ia)+dphi(2*k,ia))
      dm(k,ia)=0.5d0*(dm(2*k-1,ia)+dm(2*k,ia))

      dphi_s(k,ia)=0.5d0*(dphi_s(2*k-1,ia)+dphi_s(2*k,ia))
      dm_s(k,ia)=0.5d0*(dm_s(2*k-1,ia)+dm_s(2*k,ia))
    end do
  end do

  do ia=0,a_size
    do k=t_size/2+1,t_size
      dphi(k,ia)=(phi(2*k,ia)+4.0d0*phi(2*k-1,ia)+phi(2*k-2,ia))/6.0d0
      dm(k,ia)=(m(2*k,ia)+4.0d0*m(2*k-1,ia)+m(2*k-2,ia))/6.0d0

      dphi_s(k,ia)=(phi_s(2*k,ia)+4.0d0*phi_s(2*k-1,ia)+phi_s(2*k-2,ia))/6.0d0
      dm_s(k,ia)=(m_s(2*k,ia)+4.0d0*m_s(2*k-1,ia)+m_s(2*k-2,ia))/6.0d0
    end do
  end do

  !Compress the calculation space
  !=================================================
  do ia=0,a_size
    do k=1,t_size
      phi(k,ia)=phi(2*k,ia)
      m(k,ia)=m(2*k,ia)

      phi_s(k,ia)=phi_s(2*k,ia)
      m_s(k,ia)=m_s(2*k,ia)
    end do
  end do

  !The same for the msd variables
  do k=1,t_size/2
    ddr2(k)=0.5d0*(ddr2(2*k-1)+ddr2(2*k))
    dm_msd(k)=0.5d0*(dm_msd(2*k-1)+dm_msd(2*k))
  end do
  do k=t_size/2+1,t_size
    ddr2(k)=(dr2(2*k)+4.0d0*dr2(2*k-1)+dr2(2*k-2))/6.0d0
    dm_msd(k)=(m_msd(2*k)+4.0d0*m_msd(2*k-1)+m_msd(2*k-2))/6.0d0
  end do

  do k=1,t_size
    dr2(k)=dr2(2*k)
    m_msd(k)=m_msd(2*k)
  end do

  k=t_size
  mult=2.0d0*mult

end do

end subroutine
