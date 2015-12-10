subroutine dyn (density,tlimit,phi_init,p2)

implicit none

double precision::density
double precision::tlimit
character(len=4)::phi_init
integer::p2

integer::t,l

integer::ip,iq,ik,ir
double precision::pr,qr,kr

double precision::a,b,c

integer,parameter::t_size=100
integer,parameter::q_size=299

double precision::factv1,factv2,factw1,factw2,factu2,factu1
double precision,dimension(0:q_size)::factd1
double precision,dimension(0:q_size)::factd1_s
double precision,parameter::pi=4.0d0*datan(1.0d0)
double precision::conv=1.0d0 !Convergence of the calculation
double precision,parameter::prec=1.0d-12 !Precision of the calculation

double precision::ht
double precision::hq
double precision::mult

integer::tt
double precision::iter

!To restart the calculation

logical::ex=.false. !to test the existence of a file or directory
integer::clines_dyn
integer::nlines
integer::recline=0
integer::io=0

double precision,dimension(0:q_size)::q
double precision,dimension(0:q_size)::ccq
double precision,dimension(0:q_size)::hcq
double precision,dimension(0:q_size)::hdq
double precision,dimension(0:q_size)::scq
double precision,dimension(0:q_size)::sdq

double precision,parameter::tmic=160.0d0
double precision,dimension(0:q_size)::tau=0.0d0
double precision,dimension(0:q_size)::tau_s=0.0d0

double precision,dimension(1:2*t_size,0:q_size)::phi=0.0d0
double precision,dimension(1:2*t_size,0:q_size)::dphi=0.0d0

!For the convergence of phi and m in the subiteration
double precision,dimension(0:q_size)::phit=0.0d0
double precision,dimension(0:q_size)::phit2=0.0d0
double precision::conv_phit=1.0d0
double precision::conv_phit_s=1.0d0
double precision::prec_phit=1.0d-12
double precision,dimension(0:q_size)::mt=0.0d0

double precision,dimension(1:2*t_size,0:q_size)::m=0.0d0
double precision,dimension(1:2*t_size,0:q_size)::dm=0.0d0
double precision,dimension(0:q_size)::m0=0.0d0

double precision,dimension(0:q_size)::ct=0.0d0
double precision,dimension(0:q_size)::dt=0.0d0

double precision,dimension(0:q_size,0:q_size,0:q_size)::v1,v2

!For the calculation of the self dynamics
double precision,dimension(1:2*t_size,0:q_size)::phi_s=0.0d0
double precision,dimension(1:2*t_size,0:q_size)::dphi_s=0.0d0

double precision,dimension(0:q_size)::phit_s=0.0d0
double precision,dimension(0:q_size)::phit_s2=0.0d0

double precision,dimension(0:q_size)::mt_s=0.0d0

double precision,dimension(1:2*t_size,0:q_size)::m_s=0.0d0
double precision,dimension(1:2*t_size,0:q_size)::dm_s=0.0d0
double precision,dimension(0:q_size)::m_s0=0.0d0

double precision,dimension(0:q_size)::ct_s=0.0d0
double precision,dimension(0:q_size)::dt_s=0.0d0

double precision,dimension(0:q_size,0:q_size,0:q_size)::v1_s,v2_s

!For the calculation of the mean square displacement
double precision,dimension(0:q_size)::w1,w2
double precision::d0=1.0d0/tmic

double precision::ct_msd=0.0d0
double precision::dt_msd=0.0d0

double precision,dimension(1:2*t_size)::m_msd=0.0d0
double precision,dimension(1:2*t_size)::dm_msd=0.0d0
double precision::m_msd0=0.0d0

double precision,dimension(1:2*t_size)::dr2=0.0d0
double precision,dimension(1:2*t_size)::ddr2=0.0d0

double precision,dimension(1:2*t_size)::mu

!For the calculation of the mean quartic displacement, and then the non-gaussian parameter
double precision,dimension(0:q_size)::u1,u2
double precision,dimension(1:2*t_size)::dr4=0.0d0
double precision,dimension(1:2*t_size)::ddr4=0.0d0
double precision::m_mqd0=0.0d0
double precision,dimension(1:2*t_size)::m_mqd=0.0d0
double precision,dimension(1:2*t_size)::dm_mqd=0.0d0
double precision::dt1_mqd,dt2_mqd
double precision::ct1_mqd,ct2_mqd
double precision,dimension(1:2*t_size,0:q_size)::dqphi_s=0.0d0
double precision,dimension(1:2*t_size,0:q_size)::d2qphi_s=0.0d0
double precision,dimension(1:2*t_size)::alpha=0.0d0

integer::wits

character(len=20)::phi_file
character(len=20)::phi_s_file
character(len=20)::dr2_file
character(len=20)::dr4_file
character(len=20)::alpha_file
character(len=20)::mu_file

character(len=25)::phi_inpt_file
character(len=25)::phi_s_inpt_file
character(len=25)::dr2_inpt_file
character(len=25)::dr4_inpt_file

call fileman('ccq.dat',7,11,1)
call fileman('hcq.dat',7,12,1)
call fileman('hdq.dat',7,13,1)

!do iq=0,q_size
!  do ir=0,(2**(p2-11)-1)
!    read(11,*) a,b
!    read(12,*) a,b
!    read(13,*) a,b
!  end do
!  read(11,*) q(iq),ccq(iq)
!  read(12,*) q(iq),hcq(iq)
!  read(13,*) q(iq),hdq(iq)
!end do

do iq=0,q_size
  read(11,*) a,b
  read(12,*) a,b
  read(13,*) a,b
 
  read(11,*) q(iq),ccq(iq)
  read(12,*) q(iq),hcq(iq)
  read(13,*) q(iq),hdq(iq)
end do

call fileman('ccq.dat',7,11,0)
call fileman('hcq.dat',7,12,0)
call fileman('hdq.dat',7,13,0)

do iq=0,q_size
  scq(iq)=1.0d0+density*hcq(iq)
  sdq(iq)=density*hdq(iq)

  tau(iq)=(tmic*scq(iq))/(q(iq)**2)
  tau_s(iq)=tmic/(q(iq)**2)
end do

ht=1.0d-4*minval(tau,dim=1)/dble(t_size) !the pitch on the time

hq=q(2)-q(1) !the pitch on the q coord

factv2=density*(hq**3)/(32.0d0*(pi**2))
factv1=(hq**3)/(16.0d0*(pi**2))

factw2=((hq**5)*density)/(6.0d0*pi**2)
factw1=(hq**5)/(6.0d0*pi**2)

factu2=((hq**5)*density)/(20.0d0*pi**2)
factu1=(hq**5)/(20.0d0*pi**2)

do ip=0,q_size
  pr=dble(ip)+0.5d0
  do ik=0,q_size
    kr=dble(ik)+0.5d0
    do iq=0,q_size
      qr=dble(iq)+0.5d0

      if ((ip .ge. abs(iq-ik)) .and. (ip .le. (iq+ik))) then

        v2(iq,ik,ip)=factv2*scq(iq)*scq(ik)*scq(ip)*&
        (kr*pr/(qr**5))*(((kr**2)+(qr**2)-(pr**2))*ccq(ik)+((pr**2)+(qr**2)-(kr**2))*ccq(ip))**2
        v1(iq,ik,ip)=factv1*scq(iq)*scq(ik)*hdq(ip)*&
        (kr*pr/(qr**5))*(((kr**2)+(qr**2)-(pr**2))*ccq(ik)*density+((pr**2)+(qr**2)-(kr**2)))**2

        v2_s(iq,ik,ip)=2.0d0*factv2*(ccq(ip)**2)*scq(ip)*(kr*pr/(qr**5))*((pr**2)+(qr**2)-(kr**2))**2
        v1_s(iq,ik,ip)=factv1*hdq(ip)*(kr*pr/(qr**5))*((pr**2)+(qr**2)-(kr**2))**2

      else

        v2(iq,ik,ip)=0.0d0
        v1(iq,ik,ip)=0.0d0

        v2_s(iq,ik,ip)=0.0d0
        v1_s(iq,ik,ip)=0.0d0

      end if

    end do
  end do

  w2(ip)=factw2*(pr**4)*(ccq(ip)**2)*scq(ip)
  w1(ip)=factw1*(pr**4)*hdq(ip)

  u2(ip)=factu2*(pr**4)*(ccq(ip)**2)*scq(ip)
  u1(ip)=factu1*(pr**4)*hdq(ip)

end do


!First calculation of phi and its derivative
!=================================================
do iq=0,q_size
  do t=1,t_size
    phi(t,iq)=1.0d0-(ht*dble(t))/tau(iq)
    phi_s(t,iq)=1.0d0-(ht*dble(t))/tau_s(iq)
  end do
end do

!Calculation of the mean square displacement
do t=1,t_size
  dr2(t)=6.0d0*d0*ht*dble(t)
  dr4(t)=60.0d0*(d0*ht*dble(t))**2
  alpha(t)=(3.0d0/5.0d0)*(dr4(t)/(dr2(t)**2))-1.0d0
end do

mu(1)=1.0d0
do t=2,t_size
  mu(t)=(dlog(dr2(t))-dlog(dr2(t-1)))/(dlog(ht*dble(t))-dlog(ht*dble(t-1)))
end do

!Integrals
do iq=0,q_size
  dphi(1,iq)=0.5d0*(1.0d0+phi(1,iq))
  dphi_s(1,iq)=0.5d0*(1.0d0+phi_s(1,iq))
  do t=2,t_size
    dphi(t,iq)=0.5d0*(phi(t-1,iq)+phi(t,iq))
    dphi_s(t,iq)=0.5d0*(phi_s(t-1,iq)+phi_s(t,iq))
  end do
end do

ddr2(1)=0.5d0*dr2(1)
ddr4(1)=0.5d0*dr4(1)
do t=2,t_size
  ddr2(t)=0.5d0*(dr2(t-1)+dr2(t))
  ddr4(t)=0.5d0*(dr4(t-1)+dr4(t))
end do

!Dervatives along q
do t=1,t_size
  dqphi_s(t,0)=-(1.0d0/3.0d0)*q(0)*dr2(t)*dexp(-((q(0)**2)*dr2(t))/6.0d0)
  dqphi_s(t,1)=-(1.0d0/3.0d0)*q(1)*dr2(t)*dexp(-((q(1)**2)*dr2(t))/6.0d0)
end do
do iq=2,q_size-1
  do t=1,t_size
    dqphi_s(t,iq)=(phi_s(t,iq+1)-phi_s(t,iq-1))/(2.0d0*hq)
  end do
end do
do t=1,t_size
  dqphi_s(t,q_size)=-phi_s(t,q_size-1)/(2.0d0*hq)
end do

do t=1,t_size
  d2qphi_s(t,0)=((1.0d0/9.0d0)*((q(0)*dr2(t))**2)-(1.0d0/3.0d0)*dr2(t))*dexp(-((q(0)**2)*dr2(t))/6.0d0)
  d2qphi_s(t,1)=((1.0d0/9.0d0)*((q(1)*dr2(t))**2)-(1.0d0/3.0d0)*dr2(t))*dexp(-((q(1)**2)*dr2(t))/6.0d0)
end do
do iq=2,q_size-1
  do t=1,t_size
    d2qphi_s(t,iq)=(phi_s(t,iq+1)-2.0d0*phi_s(t,iq)+phi_s(t,iq-1))/(hq**2)
  end do
end do
do t=1,t_size
  d2qphi_s(t,q_size)=(-2.0d0*phi_s(t,q_size)+phi_s(t,q_size-1))/(hq**2)
end do


!First calculation of m and its derivative
!=================================================
do iq=0,q_size
  do t=1,t_size
    m(t,iq)=0.0d0
    m_s(t,iq)=0.0d0
  end do
end do

do ip=0,q_size
  do ik=0,q_size
    do iq=0,q_size
      m0(iq)=m0(iq)+v2(iq,ik,ip)+v1(iq,ik,ip)
      m_s0(iq)=m_s0(iq)+v2_s(iq,ik,ip)+v1_s(iq,ik,ip)
      do t=1,t_size
        m(t,iq)=m(t,iq)+v2(iq,ik,ip)*phi(t,ip)*phi(t,ik)+v1(iq,ik,ip)*phi(t,ik)

        m_s(t,iq)=m_s(t,iq)+v2_s(iq,ik,ip)*phi_s(t,ik)*phi(t,ip)+v1_s(iq,ik,ip)*phi_s(t,ik)
      end do
    end do
  end do
end do

m_msd0=0.0d0
m_mqd0=0.0d0
do t=1,t_size
  m_msd(t)=0.0d0
  m_mqd(t)=0.0d0
end do
do iq=0,q_size
  m_msd0=m_msd0+w2(iq)+w1(iq)
  do t=1,t_size
    m_msd(t)=m_msd(t)+w2(iq)*phi_s(t,iq)*phi(t,iq)+w1(iq)*phi_s(t,iq)
    m_mqd(t)=m_mqd(t)+(u2(iq)*phi(t,iq)+u1(iq))*((2.0d0*dqphi_s(t,iq))/(3.0d0*q(iq))+d2qphi_s(t,iq))
  end do
end do

do iq=0,q_size
  dm(1,iq)=0.5d0*(m0(iq)+m(1,iq))
  dm_s(1,iq)=0.5d0*(m_s0(iq)+m_s(1,iq))
  dm_msd(1)=0.5d0*(m_msd0+m_msd(1))
  dm_mqd(1)=0.5d0*(m_mqd0+m_mqd(1))
  do t=2,t_size
    dm(t,iq)=0.5d0*(m(t-1,iq)+m(t,iq))
    dm_s(t,iq)=0.5d0*(m_s(t-1,iq)+m_s(t,iq))
    dm_msd(t)=0.5d0*(m_msd(t-1)+m_msd(t))
    dm_mqd(t)=0.5d0*(m_mqd(t-1)+m_mqd(t))
  end do
end do

!If the calculation was already made, store the old stuff in dyn_inpt
if (phi_init .eq. 'file') then
!  inquire (directory='dyn', exist=ex)
  inquire (file='dyn/phi_000.dat', exist=ex)
  if (ex .eqv. .true.) then
    nlines=clines_dyn()
    call system ('mv dyn/ dyn_inpt/')
    do iq=0,q_size
      write (phi_inpt_file,'(a13,i3.3,a4)') 'dyn_inpt/phi_',iq,'.dat'
      write (phi_s_inpt_file,'(a15,i3.3,a4)') 'dyn_inpt/phi_s_',iq,'.dat'
      call fileman(phi_inpt_file,len(phi_inpt_file),710+iq,1)
      call fileman(phi_s_inpt_file,len(phi_s_inpt_file),1010+iq,1)
    end do

    !read the first value at time=0
    do iq=0,q_size
      read (710+iq,*) a,phi(t,iq)
      read (1010+iq,*) a,phi_s(t,iq)
    end do

    !read the rest
    do t=1,t_size
      do iq=0,q_size
        read (710+iq,*) a,phi(t,iq)
        read (1010+iq,*) a,phi_s(t,iq)
      end do

      !and count the number of lines
      recline=recline+1
    end do

  else
    write (6,*) 'No files for phi, please check your inputs.'
    stop
  end if
end if


!Create a folder to stock all the phi files
call system ('mkdir dyn/')

write (6,'(a47)') "          time           ||         convergence"

do iq=0,q_size
  write (phi_file,'(a8,i3.3,a4)') 'dyn/phi_',iq,'.dat'
  write (phi_s_file,'(a10,i3.3,a4)') 'dyn/phi_s_',iq,'.dat'

  call fileman(phi_file,len(phi_file),10+iq,1)
  call fileman(phi_s_file,len(phi_s_file),310+iq,1)

  write (10+iq,*) 0.0d0,1.0d0
  write (310+iq,*) 0.0d0,1.0d0
  do t=1,t_size
    write (10+iq,*) ht*dble(t),phi(t,iq)
    write (310+iq,*) ht*dble(t),phi_s(t,iq)
  end do
  call fileman(phi_file,len(phi_file),10+iq,0)
  call fileman(phi_s_file,len(phi_s_file),310+iq,0)
end do

write (dr2_file,'(a11)') 'dyn/dr2.dat'
write (dr4_file,'(a11)') 'dyn/dr4.dat'
write (alpha_file,'(a13)') 'dyn/alpha.dat'
write (mu_file,'(a10)') 'dyn/mu.dat'
call fileman(dr2_file,len(dr2_file),610,1)
call fileman(dr4_file,len(dr4_file),611,1)
call fileman(alpha_file,len(alpha_file),612,1)
call fileman(mu_file,len(mu_file),613,1)
write (610,*) 0.0d0,0.0d0
write (611,*) 0.0d0,0.0d0
write (612,*) 0.0d0,0.0d0
write (613,*) 0.0d0,1.0d0
do t=1,t_size
  write (610,*) ht*dble(t),dr2(t)
  write (611,*) ht*dble(t),dr4(t)
  write (612,*) ht*dble(t),alpha(t)
  write (613,*) ht*dble(t),mu(t)
end do
call fileman(dr2_file,len(dr2_file),610,0)
call fileman(dr4_file,len(dr4_file),611,0)
call fileman(alpha_file,len(alpha_file),612,0)
call fileman(mu_file,len(mu_file),613,0)



iter=dble(t_size)
mult=1.0d0

!recline=100 !for the restart, counts the line at which we read

!ITERATION ON THE CONVERGENCE OF PHI(K)
do while ((conv .gt. prec) .or. (ht*dble(iter) .lt. tlimit))

  !ITERATION INCREMENTING K
  do t=t_size+1,t_size*2

    tt=t/2
    iter=iter+mult

    if (phi_init .eq. 'none') then

      !Calculate C
      !=================================================
      do iq=0,q_size
        ct  (iq)=m  (tt,iq)*phi  (t-tt,iq)-m  (t-1,iq)*dphi  (1,iq)-phi  (t-1,iq)*dm  (1,iq)
        ct_s(iq)=m_s(tt,iq)*phi_s(t-tt,iq)-m_s(t-1,iq)*dphi_s(1,iq)-phi_s(t-1,iq)*dm_s(1,iq)
      end do
      do iq=0,q_size
        do l=2,tt
          ct  (iq)=ct  (iq)+(phi  (t-l+1,iq)-phi  (t-l,iq))*dm  (l,iq)
          ct_s(iq)=ct_s(iq)+(phi_s(t-l+1,iq)-phi_s(t-l,iq))*dm_s(l,iq)
        end do
      end do
      do iq=0,q_size
        do l=2,t-tt
          ct  (iq)=ct  (iq)+(m  (t-l+1,iq)-m  (t-l,iq))*dphi  (l,iq)
          ct_s(iq)=ct_s(iq)+(m_s(t-l+1,iq)-m_s(t-l,iq))*dphi_s(l,iq)
        end do
      end do

      !Calculate D
      !=================================================
      do iq=0,q_size
        factd1(iq)=tau(iq)/(2.0d0*ht*mult)
        factd1_s(iq)=tau_s(iq)/(2.0d0*ht*mult)

        dt(iq)=  ct(iq)  -factd1(iq)  *(4.0d0*phi(t-1,iq)  -phi(t-2,iq))
        dt_s(iq)=ct_s(iq)-factd1_s(iq)*(4.0d0*phi_s(t-1,iq)-phi_s(t-2,iq))

        !take the last value of phi for a good guess
        phit(iq)=phi(t-1,iq)
        phit_s(iq)=phi_s(t-1,iq)
      end do

      !Iterate to converge phi
      !=================================================
      conv_phit=1.0d0
      do while (conv_phit .gt. prec_phit)

        do iq=0,q_size
          mt(iq)=0.0d0
        end do
        do ip=0,q_size
          do ik=0,q_size
            do iq=0,q_size
              mt(iq)=mt(iq)+v2(iq,ik,ip)*phit(ip)*phit(ik)+v1(iq,ik,ip)*phit(ik)
            end do
          end do
        end do

        do iq=0,q_size
          phit2(iq)=(mt(iq)*(1.0d0-dphi(1,iq))-dt(iq))/(1.0d0+dm(1,iq)+3.0d0*factd1(iq))
        end do

        conv_phit=0.0d0
        do iq=0,q_size
          conv_phit=conv_phit+dabs(phit2(iq)-phit(iq))
        end do

        do iq=0,q_size
          phit(iq)=phit2(iq)
        end do

      end do

      !Iterate to converge phi_s
      !=================================================
      conv_phit_s=1.0d0
      do while (conv_phit_s .gt. prec_phit)

        do iq=0,q_size
          mt_s(iq)=0.0d0
        end do

        do ip=0,q_size
          do ik=0,q_size
            do iq=0,q_size
              mt_s(iq)=mt_s(iq)+v2_s(iq,ik,ip)*phit_s(ik)*phit(ip)+v1_s(iq,ik,ip)*phit_s(ik)
            end do
          end do
        end do

        do iq=0,q_size
          phit_s2(iq)=(mt_s(iq)*(1.0d0-dphi_s(1,iq))-dt_s(iq))/(1.0d0+dm_s(1,iq)+3.0d0*factd1_s(iq))
        end do

        conv_phit_s=0.0d0
        do iq=0,q_size
          conv_phit_s=conv_phit_s+dabs(phit_s2(iq)-phit_s(iq))
        end do

        do iq=0,q_size
          phit_s(iq)=phit_s2(iq)
        end do

      end do

      !Calculate the memory functions one last time
      do iq=0,q_size
         mt(iq)=0.0d0
         mt_s(iq)=0.0d0
      end do
      do ip=0,q_size
        do ik=0,q_size
          do iq=0,q_size
            mt(iq)=mt(iq)+v2(iq,ik,ip)*phit(ip)*phit(ik)+v1(iq,ik,ip)*phit(ik)
            mt_s(iq)=mt_s(iq)+v2_s(iq,ik,ip)*phit_s(ik)*phit(ip)+v1_s(iq,ik,ip)*phit_s(ik)
          end do
        end do
      end do
        
      !Pass the values of phit, mk, mk_s and phit_s to the arrays
      do iq=0,q_size
        phi(t,iq)=phit(iq)
        phi_s(t,iq)=phit_s(iq)

        m(t,iq)=mt(iq)
        m_s(t,iq)=mt_s(iq)
      end do
     
      !Calculation of the mean squared displacement
      m_msd(t)=0.0d0
      do iq=0,q_size
        m_msd(t)=m_msd(t)+w2(iq)*phi_s(t,iq)*phi(t,iq)+w1(iq)*phi_s(t,iq)
      end do

      ct_msd=m_msd(tt)*dr2(t-tt)-m_msd(t-1)*ddr2(1)-dr2(t-1)*dm_msd(1)
      do l=2,tt
        ct_msd=ct_msd+(dr2(t-l+1)-dr2(t-l))*dm_msd(l)
      end do
      do l=2,t-tt
        ct_msd=ct_msd+(m_msd(t-l+1)-m_msd(t-l))*ddr2(l)
      end do

      dt_msd=m_msd(t)*ddr2(1)-(4.0d0*dr2(t-1)-dr2(t-2))/(2.0d0*d0*ht*mult)+ct_msd
    
      dr2(t)=(6.0d0*d0-d0*dt_msd)/(d0*dm_msd(1)+3.0d0/(2.0d0*ht*mult))

      mu(t)=(dlog(dr2(t))-dlog(dr2(t-1)))/(dlog(ht*iter)-dlog(ht*(iter-mult)))

      !Dervatives of first and second order of phi_s along q
      dqphi_s(t,0)=-(1.0d0/3.0d0)*q(0)*dr2(t)*dexp(-((q(0)**2)*dr2(t))/6.0d0)
      dqphi_s(t,1)=-(1.0d0/3.0d0)*q(1)*dr2(t)*dexp(-((q(1)**2)*dr2(t))/6.0d0)
      dqphi_s(t,2)=(phi_s(t,3)-dexp(-((q(1)**2)*dr2(t))/6.0d0))/(2.0d0*hq)
      do iq=3,q_size-1
        dqphi_s(t,iq)=(phi_s(t,iq+1)-phi_s(t,iq-1))/(2.0d0*hq)
      end do
      dqphi_s(t,q_size)=-phi_s(t,q_size-1)/(2.0d0*hq)

      d2qphi_s(t,0)=((1.0d0/9.0d0)*((q(0)*dr2(t))**2)-(1.0d0/3.0d0)*dr2(t))*dexp(-((q(0)**2)*dr2(t))/6.0d0)
      d2qphi_s(t,1)=((1.0d0/9.0d0)*((q(1)*dr2(t))**2)-(1.0d0/3.0d0)*dr2(t))*dexp(-((q(1)**2)*dr2(t))/6.0d0)
      d2qphi_s(t,iq)=(phi_s(t,3)-2.0d0*phi_s(t,2)+dexp(-((q(1)**2)*dr2(t))/6.0d0))/(hq**2)
      do iq=3,q_size-1
        d2qphi_s(t,iq)=(phi_s(t,iq+1)-2.0d0*phi_s(t,iq)+phi_s(t,iq-1))/(hq**2)
      end do
      d2qphi_s(t,q_size)=(-2.0d0*phi_s(t,q_size)+phi_s(t,q_size-1))/(hq**2)

      !Calculation of the mean quartic displacement and of the gaussian parameter
      m_mqd(t)=0.0d0
      do iq=0,q_size
        m_mqd(t)=m_mqd(t)+(u2(iq)*phi(t,iq)+u1(iq))*(((2.0d0*dqphi_s(t,iq))/(3.0d0*q(iq)))+d2qphi_s(t,iq))
      end do

      ct1_mqd=m_msd(tt)*dr4(t-tt)-m_msd(t-1)*ddr4(1)-dr4(t-1)*dm_msd(1)
      do l=2,tt
        ct1_mqd=ct1_mqd+(dr4(t-l+1)-dr4(t-l))*dm_msd(l)
      end do

      do l=2,t-tt
        ct1_mqd=ct1_mqd+(m_msd(t-l+1)-m_msd(t-l))*ddr4(l)
      end do
      dt1_mqd=m_msd(t)*ddr4(1)-(4.0d0*dr4(t-1)-dr4(t-2))/(2.0d0*d0*ht*mult)+ct1_mqd

      ct2_mqd=m_mqd(tt)*dr2(t-tt)
      do l=1,tt
      ct2_mqd=ct2_mqd+(dr2(t-l+1)-dr2(t-l))*dm_mqd(l)
      end do
      do l=1,t-tt
      ct2_mqd=ct2_mqd+(m_mqd(t-l+1)-m_mqd(t-l))*ddr2(l)
      end do
      dt2_mqd=ct2_mqd+dr2(t)
      
      dr4(t)=(20.0d0*d0*dt2_mqd-d0*dt1_mqd)/(d0*dm_msd(1)+3.0d0/(2.0d0*ht*mult))

      alpha(t)=(3.0d0/5.0d0)*(dr4(t)/(dr2(t)**2))-1.0d0


!***IF READING PHI FROM A FILE***
!_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
    else if (phi_init .eq. 'file') then

!Read the phi and phi_s
      do iq=0,q_size
        read (710+iq,*) a,phi(t,iq)
        read (1010+iq,*) a,phi_s(t,iq)
      end do

!Calculate m, m_s, dr2, dr4 and alpha
      do iq=0,q_size
        m(t,iq)=0.0d0
        m_s(t,iq)=0.0d0
      end do

      do ip=0,q_size
        do ik=0,q_size
          do iq=0,q_size
            m(t,iq)=m(t,iq)+v2(iq,ik,ip)*phi(t,ip)*phi(t,ik)+v1(iq,ik,ip)*phi(t,ik)
            m_s(t,iq)=m_s(t,iq)+v2_s(iq,ik,ip)*phi_s(t,ik)*phi(t,ip)+v1_s(iq,ik,ip)*phi_s(t,ik)
          end do
        end do
      end do

      m_msd(t)=0.0d0
      do iq=0,q_size
        m_msd(t)=m_msd(t)+w2(iq)*phi_s(t,iq)*phi(t,iq)+w1(iq)*phi_s(t,iq)
      end do

      ct_msd=m_msd(tt)*dr2(t-tt)-m_msd(t-1)*ddr2(1)-dr2(t-1)*dm_msd(1)
      do l=2,tt
        ct_msd=ct_msd+(dr2(t-l+1)-dr2(t-l))*dm_msd(l)
      end do
      do l=2,t-tt
        ct_msd=ct_msd+(m_msd(t-l+1)-m_msd(t-l))*ddr2(l)
      end do

      dt_msd=m_msd(t)*ddr2(1)-(4.0d0*dr2(t-1)-dr2(t-2))/(2.0d0*d0*ht*mult)+ct_msd

      dr2(t)=(6.0d0*d0-d0*dt_msd)/(d0*dm_msd(1)+3.0d0/(2.0d0*ht*mult))

      mu(t)=(dlog(dr2(t))-dlog(dr2(t-1)))/(dlog(ht*iter)-dlog(ht*(iter-mult)))

      !Dervatives of first and second order of phi_s along q
      dqphi_s(t,0)=-(1.0d0/3.0d0)*q(0)*dr2(t)*dexp(-((q(0)**2)*dr2(t))/6.0d0)
      dqphi_s(t,1)=-(1.0d0/3.0d0)*q(1)*dr2(t)*dexp(-((q(1)**2)*dr2(t))/6.0d0)
      do iq=2,q_size-1
        dqphi_s(t,iq)=(phi_s(t,iq+1)-phi_s(t,iq-1))/(2.0d0*hq)
      end do
      dqphi_s(t,q_size)=-phi_s(t,q_size-1)/(2.0d0*hq)

      d2qphi_s(t,0)=((1.0d0/9.0d0)*((q(0)*dr2(t))**2)-(1.0d0/3.0d0)*dr2(t))*dexp(-((q(0)**2)*dr2(t))/6.0d0)
      d2qphi_s(t,1)=((1.0d0/9.0d0)*((q(1)*dr2(t))**2)-(1.0d0/3.0d0)*dr2(t))*dexp(-((q(1)**2)*dr2(t))/6.0d0)
      do iq=2,q_size-1
        d2qphi_s(t,iq)=(phi_s(t,iq+1)-2.0d0*phi_s(t,iq)+phi_s(t,iq-1))/(hq**2)
      end do
      d2qphi_s(t,q_size)=(-2.0d0*phi_s(t,q_size)+phi_s(t,q_size-1))/(hq**2)

      m_mqd(t)=0.0d0
      do iq=0,q_size
        m_mqd(t)=m_mqd(t)+(u2(iq)*phi(t,iq)+u1(iq))*(((2.0d0*dqphi_s(t,iq))/(3.0d0*q(iq)))+d2qphi_s(t,iq))
      end do

      ct1_mqd=m_msd(tt)*dr4(t-tt)-m_msd(t-1)*ddr4(1)-dr4(t-1)*dm_msd(1)
      do l=2,tt
        ct1_mqd=ct1_mqd+(dr4(t-l+1)-dr4(t-l))*dm_msd(l)
      end do
      do l=2,t-tt
        ct1_mqd=ct1_mqd+(m_msd(t-l+1)-m_msd(t-l))*ddr4(l)
      end do
      dt1_mqd=m_msd(t)*ddr4(1)-(4.0d0*dr4(t-1)-dr4(t-2))/(2.0d0*d0*ht*mult)+ct1_mqd
       
      dt2_mqd=dr2(t)+m_mqd(tt)*dr2(t-tt)
      do l=1,tt
      dt2_mqd=dt2_mqd+(dr2(t-l+1)-dr2(t-l))*dm_mqd(l)
      end do
      do l=1,t-tt
      dt2_mqd=dt2_mqd+(m_mqd(t-l+1)-m_mqd(t-l))*ddr2(l)
      end do
      
      dr4(t)=(20.0d0*d0*dt2_mqd-d0*dt1_mqd)/(d0*dm_msd(1)+3.0d0/(2.0d0*ht*mult))
 
      alpha(t)=(3.0d0/5.0d0)*(dr4(t)/(dr2(t)**2))-1.0d0

      recline=recline+1

!write (6,*) recline, nlines
 
      if (recline .eq. nlines-1) then
        phi_init='none'
        do iq=0,q_size
          call fileman(phi_inpt_file,len(phi_inpt_file),710+iq,0)
          call fileman(phi_s_inpt_file,len(phi_s_inpt_file),1010+iq,0)
        end do
        goto 11
      end if

    end if



!Write the stuff
!=================================================
11  do iq=0,q_size
      write (phi_file,'(a8,i3.3,a4)') 'dyn/phi_',iq,'.dat'
      write (phi_s_file,'(a10,i3.3,a4)') 'dyn/phi_s_',iq,'.dat'

      call fileman(phi_file,len(phi_file),10+iq,2)
      call fileman(phi_s_file,len(phi_s_file),310+iq,2)

      write (10+iq,*) ht*iter,phi(t,iq)
      write (310+iq,*) ht*iter,phi_s(t,iq)

      call fileman(phi_file,len(phi_file),10+iq,0)
      call fileman(phi_s_file,len(phi_s_file),310+iq,0)
    end do

    write (dr2_file,'(a11)') 'dyn/dr2.dat'
    write (dr4_file,'(a11)') 'dyn/dr4.dat'
    write (alpha_file,'(a13)') 'dyn/alpha.dat'
    call fileman(dr2_file,len(dr2_file),610,2)
    call fileman(dr4_file,len(dr4_file),611,2)
    call fileman(alpha_file,len(alpha_file),612,2)
    call fileman(mu_file,len(mu_file),613,2)
    write (610,*) ht*iter,dr2(t)
    write (611,*) ht*iter,dr4(t)
    write (612,*) ht*iter,alpha(t)
    write (613,*) ht*iter,mu(t)
    call fileman(dr2_file,len(dr2_file),610,0)
    call fileman(dr4_file,len(dr4_file),611,0)
    call fileman(alpha_file,len(alpha_file),612,0)
    call fileman(mu_file,len(mu_file),613,0)

    conv=0.0d0
    do iq=0,q_size
      conv=conv+abs(phi(t,iq)-phi(t-1,iq))+abs(phi_s(t,iq)-phi_s(t-1,iq))
    end do

    write (6,'(es24.16,a4,es24.16)') ht*mult*iter," || ",conv


  end do

  !Calculate the integrals
  !=================================================
  do iq=0,q_size
    do t=1,t_size/2
      dphi(t,iq)=0.5d0*(dphi(2*t-1,iq)+dphi(2*t,iq))
      dm(t,iq)=0.5d0*(dm(2*t-1,iq)+dm(2*t,iq))

      dphi_s(t,iq)=0.5d0*(dphi_s(2*t-1,iq)+dphi_s(2*t,iq))
      dm_s(t,iq)=0.5d0*(dm_s(2*t-1,iq)+dm_s(2*t,iq))
    end do
  end do
  do iq=0,q_size
    do t=t_size/2+1,t_size
      dphi(t,iq)=(phi(2*t,iq)+4.0d0*phi(2*t-1,iq)+phi(2*t-2,iq))/6.0d0
      dm(t,iq)=(m(2*t,iq)+4.0d0*m(2*t-1,iq)+m(2*t-2,iq))/6.0d0

      dphi_s(t,iq)=(phi_s(2*t,iq)+4.0d0*phi_s(2*t-1,iq)+phi_s(2*t-2,iq))/6.0d0
      dm_s(t,iq)=(m_s(2*t,iq)+4.0d0*m_s(2*t-1,iq)+m_s(2*t-2,iq))/6.0d0
    end do
  end do
  do iq=0,q_size
    do t=t_size+1,2*t_size
      dphi(t,iq)=0.0d0
      dm(t,iq)=0.0d0

      dphi_s(t,iq)=0.0d0
      dm_s(t,iq)=0.0d0
    end do
  end do

  !Compress the calculation space
  !=================================================
  do iq=0,q_size
    do t=1,t_size
      phi(t,iq)=phi(2*t,iq)
      m(t,iq)=m(2*t,iq)

      phi_s(t,iq)=phi_s(2*t,iq)
      m_s(t,iq)=m_s(2*t,iq)

      dqphi_s(t,iq)=dqphi_s(2*t,iq)
      d2qphi_s(t,iq)=d2qphi_s(2*t,iq)
    end do
  end do
  do iq=0,q_size
    do t=t_size+1,2*t_size
      phi(t,iq)=0.0d0
      m(t,iq)=0.0d0

      phi_s(t,iq)=0.0d0
      m_s(t,iq)=0.0d0

      dqphi_s(t,iq)=0.0d0
      d2qphi_s(t,iq)=0.0d0
    end do
  end do

  !The same for the msd and mqd variables
  do t=1,t_size/2
    ddr2(t)=0.5d0*(ddr2(2*t-1)+ddr2(2*t))
    dm_msd(t)=0.5d0*(dm_msd(2*t-1)+dm_msd(2*t))

    ddr4(t)=0.5d0*(ddr4(2*t-1)+ddr4(2*t))
    dm_mqd(t)=0.5d0*(dm_mqd(2*t-1)+dm_mqd(2*t))
  end do
  do t=t_size/2+1,t_size
    ddr2(t)=(dr2(2*t)+4.0d0*dr2(2*t-1)+dr2(2*t-2))/6.0d0
    dm_msd(t)=(m_msd(2*t)+4.0d0*m_msd(2*t-1)+m_msd(2*t-2))/6.0d0

    ddr4(t)=(dr4(2*t)+4.0d0*dr4(2*t-1)+dr4(2*t-2))/6.0d0
    dm_mqd(t)=(m_mqd(2*t)+4.0d0*m_mqd(2*t-1)+m_mqd(2*t-2))/6.0d0
  end do
  do t=t_size+1,2*t_size
    ddr2(t)=0.0d0
    dm_msd(t)=0.0d0

    ddr4(t)=0.0d0
    dm_mqd(t)=0.0d0
  end do

  do t=1,t_size
    dr2(t)=dr2(2*t)
    m_msd(t)=m_msd(2*t)

    dr4(t)=dr4(2*t)
    m_mqd(t)=m_mqd(2*t)
  end do
  do t=t_size+1,2*t_size
    dr2(t)=0.0d0
    m_msd(t)=0.0d0

    dr4(t)=0.0d0
    m_mqd(t)=0.0d0
  end do

  mult=2.0d0*mult

end do

end subroutine
