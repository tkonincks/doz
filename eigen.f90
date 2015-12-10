subroutine eigen (trans_mode,density,eigenvalue,lambda,eigenvalue_inflex&
,lambda_inflex,flag_inflex,fast,fcutup,fcutdown,fq_init,disp_iter,p2)

implicit none


!Input/Output
!=================================================
character(len=4)::trans_mode,fq_init
double precision::density,eigenvalue,lambda,eigenvalue_inflex,lambda_inflex,fcutdown
logical::flag_inflex,fast
integer::fcutup,disp_iter,p2

!Iteration variables
!=================================================
integer::ik,ip,iq,iter,ir
double precision::kr,qr,pr
double precision::convergence

!Constants
!=================================================
double precision,parameter::pi=4.0d0*datan(1.0d0)
double precision,parameter::prec=1.0d-12 !Precision of the calculation
double precision::hq

!Arrays
!=================================================
integer,parameter::q_size=299
double precision,dimension(0:q_size)::convtab=0.0d0
integer::convloc
double precision::a,b
double precision,dimension(0:q_size,0:q_size,0:q_size)::v2
double precision,dimension(0:q_size,0:q_size,0:q_size)::v1
double precision,dimension(0:q_size,0:q_size,0:q_size)::v2_s
double precision,dimension(0:q_size,0:q_size,0:q_size)::v1_s
double precision,dimension(0:q_size)::w1
double precision,dimension(0:q_size)::w2
double precision,dimension(0:q_size)::u1
double precision,dimension(0:q_size)::u2

double precision::factv1,factv2,factw2,factw1,factu2,factu1

double precision,dimension(0:q_size)::ffq
double precision,dimension(0:q_size)::ffq_s

double precision,dimension(0:q_size)::fq
double precision,dimension(0:q_size)::fq_s

double precision,dimension(0:q_size)::fq2
double precision,dimension(0:q_size)::fq_s2

double precision,dimension(0:q_size)::dqfq_s
double precision,dimension(0:q_size)::d2qfq_s

double precision,dimension(0:q_size)::fq_inflex

double precision::oodr2 !OneOverdr2
double precision::oodr4 !OneOverdr4

double precision::superr,err2

double precision,dimension(0:q_size)::q !F-space

double precision,dimension(0:q_size)::scq
double precision,dimension(0:q_size)::sdq

double precision,dimension(0:q_size)::ccq !Will contain the values of the arrays printed with the struct subroutine
double precision,dimension(0:q_size)::hcq
double precision,dimension(0:q_size)::hdq


!Where is the singularity?
!=================================================
double precision,dimension(0:q_size,0:q_size)::v1_sing
double precision,dimension(0:q_size)::er
double precision,dimension(0:q_size)::er2
double precision,dimension(0:q_size)::el
double precision,dimension(0:q_size)::el2
double precision::n1,n2 !n1 and n2 don't have ant signification, sometimes normes, sometimes junks...
double precision::coeff_er,coeff_el !Will contain the coefficients used to norm e and ec


double precision,dimension(0:q_size,0:q_size)::v1_sing_inflex
double precision,dimension(0:q_size)::er_inflex
double precision,dimension(0:q_size)::er2_inflex
double precision,dimension(0:q_size)::el_inflex
double precision,dimension(0:q_size)::el2_inflex
double precision::n1_inflex,n2_inflex
double precision::coeff_er_inflex,coeff_el_inflex !Will contain the coefficients used to norm e and ec

!File names
!=================================================
character(len=7)::ccqfile='ccq.dat'
character(len=7)::hcqfile='hcq.dat'
character(len=7)::hdqfile='hdq.dat'
character(len=6)::fqfile='fq.dat'
character(len=8)::fq_sfile='fq_s.dat'
character(len=13)::fq_inflexfile='fq_inflex.dat'

!Initialisations
!=================================================
n1=0.0d0
n2=0.0d0
n1_inflex=0.0d0
n2_inflex=0.0d0
coeff_er=0.0d0
coeff_el=0.0d0
coeff_er_inflex=0.0d0
coeff_el_inflex=0.0d0

factv1=0.0d0
factv2=0.0d0
factw1=0.0d0
factw2=0.0d0
factu1=0.0d0
factu2=0.0d0

flag_inflex=.false.

call fileman(ccqfile,len(ccqfile),11,1)
call fileman(hcqfile,len(hcqfile),12,1)
call fileman(hdqfile,len(hdqfile),13,1)
do iq=0,q_size
  
  do ir=1,(2**(p2-11)-1)
    read(11,*) a,b
    read(12,*) a,b
    read(13,*) a,b
  end do  

  read(11,*) q(iq),ccq(iq)
  read(12,*) q(iq),hcq(iq)
  read(13,*) q(iq),hdq(iq)

end do
call fileman(ccqfile,len(ccqfile),11,0)
call fileman(hcqfile,len(hcqfile),12,0)
call fileman(hdqfile,len(hdqfile),13,0)

hq=q(2)-q(1) !h, the pitch

do ip=0,q_size
  scq(ip)=1.0d0+density*hcq(ip)
  sdq(ip)=density*hdq(ip)
  ffq(ip)=0.0d0
  fq(ip)=1.0d0
  fq_s(ip)=1.0d0
  fq2(ip)=1.0d0
  fq_s2(ip)=1.0d0
  fq_inflex(ip)=0.0d0
  er(ip)=1.0d0
  er2(ip)=0.0d0
  el(ip)=1.0d0 
  el2(ip)=0.0d0
  er_inflex(ip)=1.0d0
  er2_inflex(ip)=0.0d0
  el_inflex(ip)=1.0d0 
  el2_inflex(ip)=0.0d0
  w2(ip)=0.0d0
  w1(ip)=0.0d0
  u2(ip)=0.0d0
  u1(ip)=0.0d0
  do ik=0,q_size
    v1_sing(ik,ip)=0.0d0
    v1_sing_inflex(ik,ip)=0.0d0
    do iq=0,q_size
      v2(iq,ik,ip)=0.0d0
      v1(iq,ik,ip)=0.0d0
      v2_s(iq,ik,ip)=0.0d0
      v1_s(iq,ik,ip)=0.0d0
    end do
  end do
end do
oodr2=0.0d0
oodr4=0.0d0

!Calculate v1 and v2
!=================================================
factv2=density*(hq**3)/(32.0d0*(pi**2))
factv1=(hq**3)/(16.0d0*(pi**2))

factw2=((hq**5)*density)/(6.0d0*pi**2)
factw1=(hq**5)/(6.0d0*pi**2)

factu2=((hq**5)*density)/(20.0d0*pi**2)
factu1=(hq**5)/(20.0d0*pi**2)

if (trans_mode .eq. 'loca') then

  do ip=0,q_size
    pr=dble(ip)+0.5d0
    do ik=0,q_size
      kr=dble(ik)+0.5d0
      do iq=0,q_size
        qr=dble(iq)+0.5d0
  
        if ((ip .ge. abs(iq-ik)) .and. (ip .le. (iq+ik))) then
          v1(iq,ik,ip)=factv1*hdq(ip)*(kr*pr/(qr**5))*((pr**2)+(qr**2)-(kr**2))**2
        else
          v1(iq,ik,ip)=0.0d0
        end if

      end do
    end do
  end do


else

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

end if

!Calculate fq and fq_inflex in the case it exists 
!for the discontinuous transition
!=================================================
if (trans_mode .eq. 'disc') then

  if (fq_init .eq. 'file') then
    call fileman(fqfile,len(fqfile),11,1)
    call fileman(fq_sfile,len(fq_sfile),12,1)
    do iq=0,q_size
      read (11,*) q(iq),fq(iq)
      read (12,*) q(iq),fq_s(iq)
    end do
    call fileman(fqfile,len(fqfile),11,0)
    call fileman(fq_sfile,len(fq_sfile),12,0)

  else if (fq_init .eq. 'unit') then
    do iq=0,q_size
     fq(ip)=1.0d0
     fq_s(ip)=1.0d0
    end do

  end if
  
  !Some initializations
  !=================================================
  iter=0
  flag_inflex=.false.
  err2=10.0d0

  write (6,*) ""
  write (6,'(a38)') "subiteration ||       f(q) convergence"

  !Iterate and converge fq.dat
  !=================================================
  convergence=1.0d0
  do while(convergence .gt. prec)

    iter=iter+1

    do iq=0,q_size
      ffq(iq)=0.0d0
    end do

    !Calculatie the new fq(fq)
    !=================================================
    do ip=0,q_size
      do ik=0,q_size
        do iq=0,q_size
          ffq(iq)=ffq(iq)+v2(iq,ik,ip)*fq(ip)*fq(ik)+v1(iq,ik,ip)*fq(ik)
        end do
      end do
    end do
  
    !Calculate the new fq
    !=================================================
    do iq=0,q_size
      fq(iq)=ffq(iq)/(1.0d0+ffq(iq))
    end do

    convergence=0.0d0
    do iq=0,q_size
      convtab(iq)=dabs(fq(iq)-fq2(iq)) !store the diferences in convtab
      if (convtab(iq) .gt. convergence) convergence=convtab(iq)
      fq2(iq)=fq(iq)
    end do


    !Store the fq at the inflexion point
    !=================================================
    if ((flag_inflex .eqv. .false.) .and. (iter .gt. fcutup) .and. (convergence .gt. fcutdown)) then
    !flag_inflex is false if the inflexion point is not found yet

    !fcutdown is the value below which the inflexion point is not considered anymore
    !to prevent the finidng of inflexion points at a convergence near 1.0d-12 which is false
    !fcutup is the value below it can find an inflexion point
      if (convergence .lt. err2) then
        err2=convergence
      else
        do iq=0,q_size
          fq_inflex(iq)=fq(iq)
        end do
        write (6,'(i7,a28)') iter,"      ||           INFLEXION"
        flag_inflex=.true.
      end if
    end if

    !Write the convergence
    !=================================================
    if (modulo(iter,disp_iter) .eq. 0) write (6,'(i7,a10,es24.16)') iter,"      ||  ",convergence


    if ((flag_inflex .eqv. .true.) .and. (fast .eqv. .true.)) exit
    !No need to calculate the rest if we use the 'fast' option, just stick at the fq_inflex which is a better approximation anyway

    do iq=0,q_size
      if (fq(iq) .ne. fq(iq)) then
        write (6,*) "DYNAMICS CRASHED"
        stop
      end if
    end do

  end do

  write (6,*) ""

  convergence=1.0d0
  !Iterate and converge fq_s.dat
  !=================================================
  do while(convergence .gt. prec)

    iter=iter+1

    do iq=0,q_size
      ffq_s(iq)=0.0d0
    end do

    !Calculatie the new ffq_s
    !=================================================
    do ip=0,q_size
      do ik=0,q_size
        do iq=0,q_size
          ffq_s(iq)=ffq_s(iq)+v2_s(iq,ik,ip)*fq_s(ik)*fq(ip)+v1_s(iq,ik,ip)*fq_s(ik)
        end do
      end do
    end do
  
    !Calculate the new fq
    !=================================================
    do iq=0,q_size
      fq_s(iq)=ffq_s(iq)/(1.0d0+ffq_s(iq))
    end do

    convergence=0.0d0
    do iq=0,q_size
      convtab(iq)=dabs(fq_s(iq)-fq_s2(iq)) !store the diferences in convtab
      if (convtab(iq) .gt. convergence) convergence=convtab(iq)
      fq_s2(iq)=fq_s(iq)
    end do

    if (modulo(iter,disp_iter) .eq. 0) write (6,'(i7,a10,es24.16)') iter,"      ||  ",convergence

  end do
 
  !And calculate 1/dr2
  oodr2=0.0d0
 !Here oodr2 is the memory function
  do iq=0,q_size
    oodr2=oodr2+w2(iq)*fq_s(iq)*fq(iq)+w1(iq)*fq_s(iq)
  end do
  !and here oodr2 is 1/dr2
  oodr2=oodr2/6.0d0

  dqfq_s(0)=-(1.0d0/3.0d0)*q(0)*(1.0d0/oodr2)*dexp(-((q(0)**2)*(1.0d0/oodr2)/6.0d0))
  dqfq_s(1)=-(1.0d0/3.0d0)*q(1)*(1.0d0/oodr2)*dexp(-((q(1)**2)*(1.0d0/oodr2)/6.0d0))
  dqfq_s(2)=(fq_s(3)-dexp(-((q(1)**2)*(1.0d0/oodr2))/6.0d0))/(2.0d0*hq)
  do iq=3,q_size-1
    dqfq_s(iq)=(fq_s(iq+1)-fq_s(iq-1))/(2.0d0*hq)
  end do
  dqfq_s(q_size)=-fq_s(q_size-1)/(2.0d0*hq)

  d2qfq_s(0)=((1.0d0/9.0d0)*((q(0)*(1.0d0/oodr2))**2)-(1.0d0/3.0d0)*(1.0d0/oodr2))&
  *dexp(-((q(0)**2)*(1.0d0/oodr2))/6.0d0)
  d2qfq_s(1)=((1.0d0/9.0d0)*((q(1)*(1.0d0/oodr2))**2)-(1.0d0/3.0d0)*(1.0d0/oodr2))*&
  dexp(-((q(1)**2)*(1.0d0/oodr2))/6.0d0)
  d2qfq_s(2)=(fq_s(3)-2.0d0*fq_s(2)+dexp(-((q(1)**2)*(1.0d0/oodr2))/6.0d0))/(hq**2)
  do iq=3,q_size-1
    d2qfq_s(iq)=(fq_s(iq+1)-2.0d0*fq_s(iq)+fq_s(iq-1))/(hq**2)
  end do
  d2qfq_s(q_size)=(-2.0d0*fq_s(q_size)+fq_s(q_size-1))/(hq**2)

  oodr4=0.0d0
  !oodr4 is the memory function
  do iq=0,q_size
    oodr4=oodr4+(u2(iq)*fq(iq)+u1(iq))*((2.0d0*dqfq_s(iq))/(3.0d0*q(iq))+d2qfq_s(iq))
  end do
  !oodr4 is 1/dr4
  oodr4=(3.0d0/10.0d0)*(oodr2**2)*(1.0d0/(1.0d0+oodr4))

  !Open the output files and write the result
  !=================================================
  call fileman(fqfile,len(fqfile),11,1)
  call fileman(fq_sfile,len(fq_sfile),12,1)
  do iq=0,q_size
    write(11,*) q(iq),fq(iq)
    write(12,*) q(iq),fq_s(iq)
  end do
  call fileman(fqfile,len(fqfile),11,0)
  call fileman(fq_sfile,len(fq_sfile),12,0)
 
  write (6,*) "1/dr2(inf)=",oodr2
  write (6,*) "1/dr4(inf)=",oodr4
  write (6,*) "alpha(inf)=",(3.0d0/5.0d0)*(oodr2**2)*(1.0d0/oodr4)-1.0d0
 
  if (flag_inflex .eqv. .true.) then
    call fileman(fq_inflexfile,len(fq_inflexfile),11,1)
    do iq=0,q_size
      write(11,*) q(iq),fq_inflex(iq)
    end do
    call fileman(fq_inflexfile,len(fq_inflexfile),11,0)
  end if

end if









!FOR THE CONTINUOUS TRANSITION
if (trans_mode .eq. 'cont') then
  !Calculate the v1_sing arrays
  !=================================================
  do ip=0,q_size
    do ik=0,q_size
      do iq=0,q_size
        v1_sing(iq,ik)=v1_sing(iq,ik)+v1(iq,ik,ip)
      end do
    end do
  end do
  
  !Multiply sing_line by e and ec until the accuracy is reached
  !=================================================
  convergence=1.0d0
  do while (convergence .gt. prec)

!  write (6,*) "convergence ",convergence
 
    do iq=0,q_size
      er2(iq)=0.0d0
      el2(iq)=0.0d0
    end do
 
    do ik=0,q_size
      do iq=0,q_size
        er2(iq)=er2(iq)+er(ik)*v1_sing(iq,ik)
        el2(iq)=el2(iq)+el(ik)*v1_sing(ik,iq)
      end do
    end do
  
  !Calculate a and ac to norm e and ec
  !=================================================
    n1=0.0d0
    n2=0.0d0
    do iq=0,q_size
      n1=n1+el2(iq)*er2(iq)
      n2=n2+el2(iq)*(er2(iq)**2)
    end do
    coeff_er=n1/n2
    coeff_el=n2/(n1**2)
  
  !Normalize e and ec
  !=================================================
    do iq=0,q_size
      er2(iq)=coeff_er*er2(iq)
      el2(iq)=coeff_el*el2(iq)
    end do
  
    convergence=0.0d0
    do iq=0,q_size !Calculate the convergence
      convergence=convergence+(er2(iq)-er(iq))**2+(el2(iq)-el(iq))**2
    end do
    convergence=dsqrt(convergence)
  
  !Passing the values
  !=================================================
    do iq=0,q_size
      er(iq)=er2(iq)
      el(iq)=el2(iq)
    end do
  
  end do
  
  !Calculate lambda
  !=================================================
  lambda=0.0d0
  do ip=0,q_size
    do ik=0,q_size
      do iq=0,q_size
        lambda=lambda+el2(iq)*(v2(iq,ik,ip)+v2(iq,ip,ik))*er2(ik)*er2(ip)
      end do
    end do
  end do
 
  lambda=0.5d0*lambda

  !Calculate the eigenvalueue
  !=================================================
  do iq=0,q_size
    er2(iq)=0.0d0
    el2(iq)=0.0d0
  end do

  do ik=0,q_size
    do iq=0,q_size
      er2(iq)=er2(iq)+er(ik)*v1_sing(iq,ik)
      el2(iq)=el2(iq)+el(ik)*v1_sing(ik,iq)
    end do
  end do
  
  n1=0.0d0 
  n2=0.0d0
  do iq=0,q_size
    n1=n1+er(iq)**2
    n2=n2+er2(iq)**2
  end do
  n1=dsqrt(n1)
  n2=dsqrt(n2)
  
  eigenvalue=n2/n1
  
  
  









else if (trans_mode .eq. 'loca') then

  !Calculate the v1_sing arrays
  !=================================================
  do ip=0,q_size
    do ik=0,q_size
      do iq=0,q_size
        v1_sing(iq,ik)=v1_sing(iq,ik)+v1(iq,ik,ip)
      end do
    end do
  end do
  
  !Multiply sing_line by e and ec until the accuracy is reached
  !=================================================
  convergence=1.0d0
  do while (convergence .gt. prec)
 
    do iq=0,q_size
      er2(iq)=0.0d0
      el2(iq)=0.0d0
    end do
 
    do ik=0,q_size
      do iq=0,q_size
        er2(iq)=er2(iq)+er(ik)*v1_sing(iq,ik)
        el2(iq)=el2(iq)+el(ik)*v1_sing(ik,iq)
      end do
    end do
  
  !Calculate a and ac to norm e and ec
  !=================================================
    n1=0.0d0
    n2=0.0d0
    do iq=0,q_size
      n1=n1+el2(iq)*er2(iq)
      n2=n2+el2(iq)*(er2(iq)**2)
    end do
    coeff_er=n1/n2
    coeff_el=n2/(n1**2)
  
  !Normalize e and ec
  !=================================================
    do iq=0,q_size
      er2(iq)=coeff_er*er2(iq)
      el2(iq)=coeff_el*el2(iq)
    end do
  
    convergence=0.0d0
    do iq=0,q_size !Calculate the convergence
      convergence=convergence+(er2(iq)-er(iq))**2+(el2(iq)-el(iq))**2
    end do
    convergence=dsqrt(convergence)
  
  !Passing the values
  !=================================================
    do iq=0,q_size
      er(iq)=er2(iq)
      el(iq)=el2(iq)
    end do
  
  end do
  
  !Calculate lambda
  !=================================================
  lambda=0.0d0
  
  !Calculate the eigenvalueue
  !=================================================
  do iq=0,q_size
    er2(iq)=0.0d0
    el2(iq)=0.0d0
  end do

  do ik=0,q_size
    do iq=0,q_size
      er2(iq)=er2(iq)+er(ik)*v1_sing(iq,ik)
      el2(iq)=el2(iq)+el(ik)*v1_sing(ik,iq)
    end do
  end do
  
  n1=0.0d0 
  n2=0.0d0
  do iq=0,q_size
    n1=n1+er(iq)**2
    n2=n2+er2(iq)**2
  end do
  n1=dsqrt(n1)
  n2=dsqrt(n2)
  
  eigenvalue=n2/n1



















!FOR THE DISCONTINUOUS TRANSITION
else if (trans_mode .eq. 'disc') then
  !Calculate the v1_sing arrays
  !=================================================
  do iq=0,q_size
    do ik=0,q_size
      do ip=0,q_size
        v1_sing(iq,ik)=v1_sing(iq,ik)+(v1(iq,ik,ip)+v2(iq,ik,ip)&
        *fq(ip)+v2(iq,ip,ik)*fq(ip))*(1.0d0-fq(ik))**2
      end do
    end do
  end do
  
  !Multiply v1_sing by er and el until the accuracy is reached
  !=================================================
  convergence=1.0d0
  do while (convergence .gt. prec)
  
    do iq=0,q_size
      er2(iq)=0.0d0
      el2(iq)=0.0d0
      do ik=0,q_size
        er2(iq)=er2(iq)+er(ik)*v1_sing(iq,ik)
        el2(iq)=el2(iq)+el(ik)*v1_sing(ik,iq)
      end do
    end do
  
  !Calculate n1 and n2 to normalize er and el
  !=================================================
    n1=0.0d0
    n2=0.0d0
  
    do iq=0,q_size
      n1=n1+el2(iq)*er2(iq)
      n2=n2+el2(iq)*(1.0d0-fq(iq))*(er2(iq)**2)
    end do
  
    coeff_er=n1/n2
    coeff_el=n2/(n1**2)
  
  !Normalize er and el
  !=================================================
    do iq=0,q_size
      er2(iq)=coeff_er*er2(iq)
      el2(iq)=coeff_el*el2(iq)
    end do
  
    convergence=0.0d0
  
    do iq=0,q_size !Calculate the convergence
      convergence=convergence+(er2(iq)-er(iq))**2+(el2(iq)-el(iq))**2
    end do
  
    convergence=dsqrt(convergence)
  
  !Passing the values
  !=================================================
    do iq=0,q_size
      er(iq)=er2(iq)
      el(iq)=el2(iq)
    end do
  
  end do
  
  !Calculate lambda
  !=================================================
  lambda=0.0d0
  do iq=0,q_size
    do ik=0,q_size
      do ip=0,q_size
        lambda=lambda+el2(iq)*((1.0d0-fq(ik))**2)*((1.0d0-fq(ip))**2)&
        *(v2(iq,ik,ip)+v2(iq,ip,ik))*er2(ik)*er2(ip)
      end do
    end do
  end do
  
  lambda=0.5d0*lambda
  
  !Calculate the eigenvalueue
  !=================================================
  do iq=0,q_size
    er2(iq)=0.0d0
    el2(iq)=0.0d0
    do ik=0,q_size
      er2(iq)=er2(iq)+er(ik)*v1_sing(iq,ik)
      el2(iq)=el2(iq)+el(ik)*v1_sing(ik,iq)
    end do
  end do
  
  n1=0.0d0 
  n2=0.0d0
  
  do iq=0,q_size
    n1=n1+er(iq)**2
    n2=n2+er2(iq)**2
  end do
  
  n1=dsqrt(n1)
  n2=dsqrt(n2)
  
  eigenvalue=n2/n1
  
end if  
  
  
  
  
  
  




  
!WHEN THE INFLEXION POINT EXIST
if (flag_inflex .eqv. .true.) then

  !Calculate the v1_sing arrays
  !=================================================
  do ip=0,q_size
    do ik=0,q_size
      do iq=0,q_size
        v1_sing_inflex(iq,ik)=v1_sing_inflex(iq,ik)+(v1(iq,ik,ip)+v2(iq,ik,ip)&
        *fq_inflex(ip)+v2(iq,ip,ik)*fq_inflex(ip))*(1.0d0-fq_inflex(ik))**2
      end do
    end do
  end do
  
  !Multiply v1_sing by er and el until the accuracy is reached
  !=================================================
  convergence=1.0d0

  do while (convergence .gt. prec)
  
    do iq=0,q_size
      er2_inflex(iq)=0.0d0
      el2_inflex(iq)=0.0d0
    end do

    do ik=0,q_size
      do iq=0,q_size
        er2_inflex(iq)=er2_inflex(iq)+er_inflex(ik)*v1_sing_inflex(iq,ik)
        el2_inflex(iq)=el2_inflex(iq)+el_inflex(ik)*v1_sing_inflex(ik,iq)
      end do
    end do
  
  !Calculate n1 and n2 to normalize er and el
  !=================================================
    n1_inflex=0.0d0
    n2_inflex=0.0d0
  
    do iq=0,q_size
      n1_inflex=n1_inflex+el2_inflex(iq)*er2_inflex(iq)
      n2_inflex=n2_inflex+el2_inflex(iq)*(1.0d0-fq_inflex(iq))*(er2_inflex(iq)**2)
    end do
  
    coeff_er_inflex=n1_inflex/n2_inflex
    coeff_el_inflex=n2_inflex/(n1_inflex**2)
  
  !Normalize er and el
  !=================================================
    do iq=0,q_size
      er2_inflex(iq)=coeff_er_inflex*er2_inflex(iq)
      el2_inflex(iq)=coeff_el_inflex*el2_inflex(iq)
    end do
  
    convergence=0.0d0
  
    do iq=0,q_size !Calculate the convergence
      convergence=convergence+(er2_inflex(iq)-er_inflex(iq))**2+(el2_inflex(iq)-el_inflex(iq))**2
    end do
  
    convergence=dsqrt(convergence)
  
  !Passing the values
  !=================================================
    do iq=0,q_size
      er_inflex(iq)=er2_inflex(iq)
      el_inflex(iq)=el2_inflex(iq)
    end do
  
  end do




  !Calculate lambda
  !=================================================
  lambda_inflex=0.0d0
  do ip=0,q_size
    do ik=0,q_size
      do iq=0,q_size
        lambda_inflex=lambda_inflex+el2_inflex(iq)*((1.0d0-fq_inflex(ik))**2)*((1.0d0-fq_inflex(ip))**2)&
        *(v2(iq,ik,ip)+v2(iq,ip,ik))*er2_inflex(ik)*er2_inflex(ip)
      end do
    end do
  end do
  
  lambda_inflex=0.5d0*lambda_inflex


  !Applicating a last time the er and el vectors
  !=================================================
  do iq=0,q_size
    er2_inflex(iq)=0.0d0
    el2_inflex(iq)=0.0d0
  end do

  do ik=0,q_size
    do iq=0,q_size
      er2_inflex(iq)=er2_inflex(iq)+er_inflex(ik)*v1_sing_inflex(iq,ik)
      el2_inflex(iq)=el2_inflex(iq)+el_inflex(ik)*v1_sing_inflex(ik,iq)
    end do
  end do
  
  !Calculate the eigenvalueue
  !=================================================
  n1_inflex=0.0d0 
  n2_inflex=0.0d0
  
  do iq=0,q_size
    n1_inflex=n1_inflex+er_inflex(iq)**2
    n2_inflex=n2_inflex+er2_inflex(iq)**2
  end do
  
  n1_inflex=dsqrt(n1_inflex)
  n2_inflex=dsqrt(n2_inflex)
  
  eigenvalue_inflex=n2_inflex/n1_inflex
 
end if  
  
  
end subroutine
