subroutine eigen (trans_mode,density,eigenvalue,lambda,eigenvalue_inflex,lambda_inflex,flag_inflex,fast,init_fq)

implicit none


!Input/Output
!=================================================
character(len=4)::trans_mode,init_fq
double precision::density,eigenvalue,lambda,eigenvalue_inflex,lambda_inflex
logical::flag_inflex,fast

!Iteration variables
!=================================================
integer::ik,ip,iq,iter
double precision::kr,qr,pr
double precision::convergence
double precision::conv1,conv2

!Constants
!=================================================
double precision,parameter::pi=4.0d0*datan(1.0d0)
double precision,parameter::prec=1.0d-12 !Precision of the calculation
double precision::h
integer,parameter::p2=12

!Functions
!=================================================
integer::witb

!Arrays
!=================================================
integer,parameter::a_size=299
double precision::a,b
double precision,dimension(0:a_size,0:a_size,0:a_size)::v2 !Intermediate coefficients here,
double precision,dimension(0:a_size,0:a_size,0:a_size)::v1 !and here
double precision::factv1,factv2

double precision,dimension(0:a_size)::ffq
double precision,dimension(0:a_size)::fq
double precision,dimension(0:a_size)::fq_old
double precision,dimension(0:a_size)::fq_inflex
double precision::superr,err1,err2

double precision,dimension(0:a_size)::q !F-space

double precision,dimension(0:a_size)::scq
double precision,dimension(0:a_size)::sdq

double precision,dimension(0:a_size)::ccq !Will contain the values of the arrays printed with the struct subroutine
double precision,dimension(0:a_size)::hcq
double precision,dimension(0:a_size)::hdq


!Where is the singularity?
!=================================================
double precision,dimension(0:a_size,0:a_size)::v1_sing
double precision,dimension(0:a_size)::er
double precision,dimension(0:a_size)::er2
double precision,dimension(0:a_size)::el
double precision,dimension(0:a_size)::el2
double precision::n1,n2 !n1 and n2 don't have ant signification, sometimes normes, sometimes junks...
double precision::coeff_er,coeff_el !Will contain the coefficients used to norm e and ec


double precision,dimension(0:a_size,0:a_size)::v1_sing_inflex
double precision,dimension(0:a_size)::er_inflex
double precision,dimension(0:a_size)::er2_inflex
double precision,dimension(0:a_size)::el_inflex
double precision,dimension(0:a_size)::el2_inflex
double precision::n1_inflex,n2_inflex
double precision::coeff_er_inflex,coeff_el_inflex !Will contain the coefficients used to norm e and ec

!File names
!=================================================
character(len=7)::ccqfile='ccq.dat'
character(len=7)::hcqfile='hcq.dat'
character(len=7)::hdqfile='hdq.dat'
character(len=6)::fqfile='fq.dat'
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

flag_inflex=.false.

call fileman(ccqfile,len(ccqfile),11,1)
call fileman(hcqfile,len(hcqfile),12,1)
call fileman(hdqfile,len(hdqfile),13,1)
do iq=0,a_size
  read(11,*) a,b
  read(11,*) q(iq),ccq(iq)
  read(12,*) a,b
  read(12,*) q(iq),hcq(iq)
  read(13,*) a,b
  read(13,*) q(iq),hdq(iq)
end do
call fileman(ccqfile,len(ccqfile),11,0)
call fileman(hcqfile,len(hcqfile),12,0)
call fileman(hdqfile,len(hdqfile),13,0)

h=q(2)-q(1) !h, the pitch

do ip=0,a_size
  scq(ip)=1.0d0+density*hcq(ip)
  sdq(ip)=density*hdq(ip)
  ffq(ip)=0.0d0
  fq(ip)=1.0d0
  fq_old(ip)=1.0d0
  fq_inflex(ip)=0.0d0
  er(ip)=1.0d0
  er2(ip)=0.0d0
  el(ip)=1.0d0 
  el2(ip)=0.0d0
  er_inflex(ip)=1.0d0
  er2_inflex(ip)=0.0d0
  el_inflex(ip)=1.0d0 
  el2_inflex(ip)=0.0d0
  do ik=0,a_size
    v1_sing(ik,ip)=0.0d0
    v1_sing_inflex(ik,ip)=0.0d0
    do iq=0,a_size
      v2(iq,ik,ip)=0.0d0
      v1(iq,ik,ip)=0.0d0
    end do
  end do
end do











!Calculate v1 and v2
!=================================================
factv2=density*(h**3)/(32.0d0*(pi**2))
factv1=(h**3)/(16.0d0*(pi**2))

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















!Calculate fq and fq_inflex in the case it exists 
!for the discontinuous transition
!=================================================
if (trans_mode .eq. 'disc') then

  if (init_fq .eq. 'file') then
    call fileman(fqfile,len(fqfile),11,1)
    do iq=0,a_size
      read (11,*) q(iq),fq(iq)
    end do
    call fileman(fqfile,len(fqfile),11,0)

  else if (init_fq .eq. 'unit') then
    do iq=0,a_size
     fq(ip)=1.0d0
     fq_old(ip)=1.0d0
    end do

  end if
  
  !Some initializations
  !=================================================
  convergence=1.0d0
  iter=0
  flag_inflex=.false.
  err1=0.0d0
  err2=0.0d0
  superr=1.0d0

  write (6,*) ""
  write (6,'(a38)') "subiteration ||       f(q) convergence"


  !Iterate and converge fq.dat
  !=================================================
  do while(convergence .gt. prec)

    iter=iter+1

    do iq=0,a_size
      ffq(iq)=0.0d0
    end do

    !Calculate conv1
    !=================================================
    conv1=fq(witb(fq,a_size))


    !Calculatie the new fq(fq)
    !=================================================
    do ip=0,a_size
      do ik=0,a_size
        do iq=0,a_size
          ffq(iq)=ffq(iq)+v2(iq,ik,ip)*fq(ip)*fq(ik)
          ffq(iq)=ffq(iq)+v1(iq,ik,ip)*fq(ik)
        end do
      end do
    end do
  
    !Calculate the new fq
    !=================================================
    do iq=0,a_size
      fq(iq)=ffq(iq)/(1.0d0+ffq(iq))
    end do
  
    !Store the fq at the inflexion point
    !=================================================
    err2=0.0d0
    err1=0.0d0
    do iq=0,a_size
      err1=abs(fq_old(iq)-fq(iq))
      if (err1 .gt. err2) then
        err2=err1
      end if
    end do
  
    if (flag_inflex .eqv. .false.) then
      if (err2 .gt. superr) then
        do iq=0,a_size
          fq_inflex(iq)=fq(iq)
        end do
        write (6,'(i7,a28)') iter,"      ||           INFLEXION"
        flag_inflex=.true.
      else
       superr=err2 
      end if
    end if
  
    do iq=0,a_size
      fq_old(iq)=fq(iq)
    end do
  
    !Calculate conv2 and the convergence
    !=================================================
    conv2=fq(witb(fq,a_size))
    convergence=dabs(conv2-conv1) !calculate the convergence
  
    !Write the convergence
    !=================================================
    if (modulo(iter,10) .eq. 0) write (6,'(i7,a10,es24.16)') iter,"      ||  ",convergence

    if ((flag_inflex .eqv. .true.) .and. (fast .eqv. .true.)) exit
    !No need to calculate the rest if we use the 'fdic' option, just stick at the fq_inflex which is a better approximation anyway

  end do
  
  write (6,'(i7,a10,es24.16)') iter,"      ||  ",convergence

  !Open the output files and write the result
  !=================================================
  call fileman(fqfile,len(fqfile),11,1)
  do iq=0,a_size
    write(11,*) q(iq),fq(iq)
  end do
  call fileman(fqfile,len(fqfile),11,0)
  
  if (flag_inflex .eqv. .true.) then
    call fileman(fq_inflexfile,len(fq_inflexfile),11,1)
    do iq=0,a_size
      write(11,*) q(iq),fq_inflex(iq)
    end do
    call fileman(fq_inflexfile,len(fq_inflexfile),11,0)
  end if

end if








  
  









!FOR THE CONTINUOUS TRANSITION
if (trans_mode .eq. 'cont') then
  !Calculate the v1_sing arrays
  !=================================================
  do ip=0,a_size
    do ik=0,a_size
      do iq=0,a_size
        v1_sing(iq,ik)=v1_sing(iq,ik)+v1(iq,ik,ip)
      end do
    end do
  end do
  
  !Multiply sing_line by e and ec until the accuracy is reached
  !=================================================
  convergence=1.0d0
  do while (convergence .gt. prec)
 
    do iq=0,a_size
      er2(iq)=0.0d0
      el2(iq)=0.0d0
    end do
 
    do ik=0,a_size
      do iq=0,a_size
        er2(iq)=er2(iq)+er(ik)*v1_sing(iq,ik)
        el2(iq)=el2(iq)+el(ik)*v1_sing(ik,iq)
      end do
    end do
  
  !Calculate a and ac to norm e and ec
  !=================================================
    n1=0.0d0
    n2=0.0d0
    do iq=0,a_size
      n1=n1+el2(iq)*er2(iq)
      n2=n2+el2(iq)*(er2(iq)**2)
    end do
    coeff_er=n1/n2
    coeff_el=n2/(n1**2)
  
  !Normalize e and ec
  !=================================================
    do iq=0,a_size
      er2(iq)=coeff_er*er2(iq)
      el2(iq)=coeff_el*el2(iq)
    end do
  
    convergence=0.0d0
    do iq=0,a_size !Calculate the convergence
      convergence=convergence+(er2(iq)-er(iq))**2+(el2(iq)-el(iq))**2
    end do
    convergence=dsqrt(convergence)
  
  !Passing the values
  !=================================================
    do iq=0,a_size
      er(iq)=er2(iq)
      el(iq)=el2(iq)
    end do
  
  end do
  
  !Calculate lambda
  !=================================================
  lambda=0.0d0
  do ip=0,a_size
    do ik=0,a_size
      do iq=0,a_size
        lambda=lambda+el2(iq)*(v2(iq,ik,ip)+v2(iq,ip,ik))*er2(ik)*er2(ip)
      end do
    end do
  end do
 
  lambda=0.5d0*lambda
  
  !Calculate the eigenvalueue
  !=================================================
  do iq=0,a_size
    er2(iq)=0.0d0
    el2(iq)=0.0d0
  end do

  do ik=0,a_size
    do iq=0,a_size
      er2(iq)=er2(iq)+er(ik)*v1_sing(iq,ik)
      el2(iq)=el2(iq)+el(ik)*v1_sing(ik,iq)
    end do
  end do
  
  n1=0.0d0 
  n2=0.0d0
  do iq=0,a_size
    n1=n1+er(iq)**2
    n2=n2+er2(iq)**2
  end do
  n1=dsqrt(n1)
  n2=dsqrt(n2)
  
  eigenvalue=n2/n1
  
  
  









!MAKE THE CHANGES FOR THE LOOPS HERE
else if (trans_mode .eq. 'loca') then

  !Calculate the v1_sing arrays
  !=================================================
  do ip=0,a_size
    do ik=0,a_size
      do iq=0,a_size
        v1_sing(iq,ik)=v1_sing(iq,ik)+v1(iq,ik,ip)
      end do
    end do
  end do
  
  !Multiply sing_line by e and ec until the accuracy is reached
  !=================================================
  convergence=1.0d0
  do while (convergence .gt. prec)
 
    do iq=0,a_size
      er2(iq)=0.0d0
      el2(iq)=0.0d0
    end do
 
    do ik=0,a_size
      do iq=0,a_size
        er2(iq)=er2(iq)+er(ik)*v1_sing(iq,ik)
        el2(iq)=el2(iq)+el(ik)*v1_sing(ik,iq)
      end do
    end do
  
  !Calculate a and ac to norm e and ec
  !=================================================
    n1=0.0d0
    n2=0.0d0
    do iq=0,a_size
      n1=n1+el2(iq)*er2(iq)
      n2=n2+el2(iq)*(er2(iq)**2)
    end do
    coeff_er=n1/n2
    coeff_el=n2/(n1**2)
  
  !Normalize e and ec
  !=================================================
    do iq=0,a_size
      er2(iq)=coeff_er*er2(iq)
      el2(iq)=coeff_el*el2(iq)
    end do
  
    convergence=0.0d0
    do iq=0,a_size !Calculate the convergence
      convergence=convergence+(er2(iq)-er(iq))**2+(el2(iq)-el(iq))**2
    end do
    convergence=dsqrt(convergence)
  
  !Passing the values
  !=================================================
    do iq=0,a_size
      er(iq)=er2(iq)
      el(iq)=el2(iq)
    end do
  
  end do
  
  !Calculate lambda
  !=================================================
  lambda=0.0d0
  
  !Calculate the eigenvalueue
  !=================================================
  do iq=0,a_size
    er2(iq)=0.0d0
    el2(iq)=0.0d0
  end do

  do ik=0,a_size
    do iq=0,a_size
      er2(iq)=er2(iq)+er(ik)*v1_sing(iq,ik)
      el2(iq)=el2(iq)+el(ik)*v1_sing(ik,iq)
    end do
  end do
  
  n1=0.0d0 
  n2=0.0d0
  do iq=0,a_size
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
  do iq=0,a_size
    do ik=0,a_size
      do ip=0,a_size
        v1_sing(iq,ik)=v1_sing(iq,ik)+(v1(iq,ik,ip)+v2(iq,ik,ip)&
        *fq(ip)+v2(iq,ip,ik)*fq(ip))*(1.0d0-fq(ik))**2
      end do
    end do
  end do
  
  !Multiply v1_sing by er and el until the accuracy is reached
  !=================================================
  convergence=1.0d0
  do while (convergence .gt. prec)
  
    do iq=0,a_size
      er2(iq)=0.0d0
      el2(iq)=0.0d0
      do ik=0,a_size
        er2(iq)=er2(iq)+er(ik)*v1_sing(iq,ik)
        el2(iq)=el2(iq)+el(ik)*v1_sing(ik,iq)
      end do
    end do
  
  !Calculate n1 and n2 to normalize er and el
  !=================================================
    n1=0.0d0
    n2=0.0d0
  
    do iq=0,a_size
      n1=n1+el2(iq)*er2(iq)
      n2=n2+el2(iq)*(1.0d0-fq(iq))*(er2(iq)**2)
    end do
  
    coeff_er=n1/n2
    coeff_el=n2/(n1**2)
  
  !Normalize er and el
  !=================================================
    do iq=0,a_size
      er2(iq)=coeff_er*er2(iq)
      el2(iq)=coeff_el*el2(iq)
    end do
  
    convergence=0.0d0
  
    do iq=0,a_size !Calculate the convergence
      convergence=convergence+(er2(iq)-er(iq))**2+(el2(iq)-el(iq))**2
    end do
  
    convergence=dsqrt(convergence)
  
  !Passing the values
  !=================================================
    do iq=0,a_size
      er(iq)=er2(iq)
      el(iq)=el2(iq)
    end do
  
  end do
  
  !Calculate lambda
  !=================================================
  lambda=0.0d0
  do iq=0,a_size
    do ik=0,a_size
      do ip=0,a_size
        lambda=lambda+el2(iq)*((1.0d0-fq(ik))**2)*((1.0d0-fq(ip))**2)&
        *(v2(iq,ik,ip)+v2(iq,ip,ik))*er2(ik)*er2(ip)
      end do
    end do
  end do
  
  lambda=0.5d0*lambda
  
  !Calculate the eigenvalueue
  !=================================================
  do iq=0,a_size
    er2(iq)=0.0d0
    el2(iq)=0.0d0
    do ik=0,a_size
      er2(iq)=er2(iq)+er(ik)*v1_sing(iq,ik)
      el2(iq)=el2(iq)+el(ik)*v1_sing(ik,iq)
    end do
  end do
  
  n1=0.0d0 
  n2=0.0d0
  
  do iq=0,a_size
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
  do ip=0,a_size
    do ik=0,a_size
      do iq=0,a_size
        v1_sing_inflex(iq,ik)=v1_sing_inflex(iq,ik)+(v1(iq,ik,ip)+v2(iq,ik,ip)&
        *fq_inflex(ip)+v2(iq,ip,ik)*fq_inflex(ip))*(1.0d0-fq_inflex(ik))**2
      end do
    end do
  end do
  
  !Multiply v1_sing by er and el until the accuracy is reached
  !=================================================
  convergence=1.0d0

  do while (convergence .gt. prec)
  
    do iq=0,a_size
      er2_inflex(iq)=0.0d0
      el2_inflex(iq)=0.0d0
    end do

    do ik=0,a_size
      do iq=0,a_size
        er2_inflex(iq)=er2_inflex(iq)+er_inflex(ik)*v1_sing_inflex(iq,ik)
        el2_inflex(iq)=el2_inflex(iq)+el_inflex(ik)*v1_sing_inflex(ik,iq)
      end do
    end do
  
  !Calculate n1 and n2 to normalize er and el
  !=================================================
    n1_inflex=0.0d0
    n2_inflex=0.0d0
  
    do iq=0,a_size
      n1_inflex=n1_inflex+el2_inflex(iq)*er2_inflex(iq)
      n2_inflex=n2_inflex+el2_inflex(iq)*(1.0d0-fq_inflex(iq))*(er2_inflex(iq)**2)
    end do
  
    coeff_er_inflex=n1_inflex/n2_inflex
    coeff_el_inflex=n2_inflex/(n1_inflex**2)
  
  !Normalize er and el
  !=================================================
    do iq=0,a_size
      er2_inflex(iq)=coeff_er_inflex*er2_inflex(iq)
      el2_inflex(iq)=coeff_el_inflex*el2_inflex(iq)
    end do
  
    convergence=0.0d0
  
    do iq=0,a_size !Calculate the convergence
      convergence=convergence+(er2_inflex(iq)-er_inflex(iq))**2+(el2_inflex(iq)-el_inflex(iq))**2
    end do
  
    convergence=dsqrt(convergence)
  
  !Passing the values
  !=================================================
    do iq=0,a_size
      er_inflex(iq)=er2_inflex(iq)
      el_inflex(iq)=el2_inflex(iq)
    end do
  
  end do




  !Calculate lambda
  !=================================================
  lambda_inflex=0.0d0
  do ip=0,a_size
    do ik=0,a_size
      do iq=0,a_size
        lambda_inflex=lambda_inflex+el2_inflex(iq)*((1.0d0-fq_inflex(ik))**2)*((1.0d0-fq_inflex(ip))**2)&
        *(v2(iq,ik,ip)+v2(iq,ip,ik))*er2_inflex(ik)*er2_inflex(ip)
      end do
    end do
  end do
  
  lambda_inflex=0.5d0*lambda_inflex


  !Applicating a last time the er and el vectors
  !=================================================
  do iq=0,a_size
    er2_inflex(iq)=0.0d0
    el2_inflex(iq)=0.0d0
  end do

  do ik=0,a_size
    do iq=0,a_size
      er2_inflex(iq)=er2_inflex(iq)+er_inflex(ik)*v1_sing_inflex(iq,ik)
      el2_inflex(iq)=el2_inflex(iq)+el_inflex(ik)*v1_sing_inflex(ik,iq)
    end do
  end do
  
  !Calculate the eigenvalueue
  !=================================================
  n1_inflex=0.0d0 
  n2_inflex=0.0d0
  
  do iq=0,a_size
    n1_inflex=n1_inflex+er_inflex(iq)**2
    n2_inflex=n2_inflex+er2_inflex(iq)**2
  end do
  
  n1_inflex=dsqrt(n1_inflex)
  n2_inflex=dsqrt(n2_inflex)
  
  eigenvalue_inflex=n2_inflex/n1_inflex
  
  if (fast .eqv. .true.) then
    eigenvalue=eigenvalue_inflex
    lambda=lambda_inflex
  end if

end if  
  
  
  end subroutine
