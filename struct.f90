subroutine struct (density,delta,sigma,closure,mix_param,cr_init,correl,p2)

implicit none

double precision::density
double precision::delta
double precision::sigma
character(len=4)::closure,cr_init,correl
double precision::mix_param

!Comment faire en sorte qu'il l'accepte ? Mettre une valeur par defaut ?
integer::p2

integer::disc
integer::io

!Iteration variables
!==================================================
integer::i
integer::nbiter

!Constants
!==================================================
double precision::prec
double precision::pi

!prec=1.0d-12
!pi=4.0d0*datan(1.0d0)

!For the initialization by existing files
!==================================================
logical::fileex

!And these for the gaussian potential of the disorder
!==================================================
double precision,dimension(0:2**p2)::k

!For the Percus-Yevick initialization
!==================================================
double precision::eta
double precision::lamb1
double precision::lamb2

!For the Waisman parametrization
!==================================================
double precision::cs=0.0d0
double precision::xred=0.0d0
double precision::fw=0.0d0
double precision::qw=0.0d0
double precision::v0=0.0d0
double precision::z1=0.0d0
double precision::sigma1=0.0d0
double precision::tau1=0.0d0
double precision::alpha1=0.0d0
double precision::v1v0=0.0d0
double precision::k1=0.0d0

!Arrays
!==================================================
double precision,dimension(0:2**p2)::r
double precision,dimension(0:2**p2)::q
double precision::dr=0.01d0
!double precision::dq=pi/((2.0d0**(p2)*dr))
!double precision::rmax=dr*(2**(p2))
double precision::dq
double precision::rmax

double precision,dimension(0:2**p2)::cr !Direct interaction
double precision,dimension(0:2**p2)::cr2 !Array used to stock the value of cr(r)
double precision,dimension(0:2**p2)::cq !The FT of cr
double precision,dimension(0:2**p2)::cdr !Direct interaction
double precision,dimension(0:2**p2)::cdr2 !Array used to stock the value of cdr(r)
double precision,dimension(0:2**p2)::cdq !The FT of cdr
double precision,dimension(0:2**p2)::ccr !Connected direct interaction
double precision,dimension(0:2**p2)::ccq !Same, in the F-space

double precision,dimension(0:2**p2)::gamr !gam = h(r) (total interaction) - c(r) (direct interaction)
double precision,dimension(0:2**p2)::gamq !The FT of gam
double precision,dimension(0:2**p2)::gamdr !gam = hd(r) (total interaction) - cd(r) (direct interaction)
double precision,dimension(0:2**p2)::gamdq !The FT of gam

double precision,dimension(0:2**p2)::hr !Total interaction
double precision,dimension(0:2**p2)::hq !Same, in the F-space
double precision,dimension(0:2**p2)::hdr !Total disorder interaction
double precision,dimension(0:2**p2)::hdq !Same, in the F-space
double precision,dimension(0:2**p2)::hcr !Connected total interaction
double precision,dimension(0:2**p2)::hcq !Same, in the F-space

double precision,dimension(0:2**p2)::scq
double precision,dimension(0:2**p2)::sdq

!These are for the EXPC and OEXP closure
double precision,dimension(0:2**p2)::hrref
double precision,dimension(0:2**p2)::crref
double precision,dimension(0:2**p2)::crref2
double precision,dimension(0:2**p2)::ccrref
double precision,dimension(0:2**p2)::cqref
double precision,dimension(0:2**p2)::gamqref
double precision,dimension(0:2**p2)::gamrref

!These are for the N3 closure
double precision,dimension(0:2**p2)::hrorpa
double precision,dimension(0:2**p2)::hdrorpa
double precision,dimension(0:2**p2)::ssr
double precision,dimension(0:2**p2)::ssdr
double precision,dimension(0:2**p2)::hhr
double precision,dimension(0:2**p2)::hhdr

double precision,dimension(0:2**p2)::ssq
double precision,dimension(0:2**p2)::ssdq
double precision,dimension(0:2**p2)::hhq
double precision,dimension(0:2**p2)::hhdq

double precision::convergence !Used to calculate the convergence of the calcultaion

!Functions
!==================================================
double precision::gaussian
double precision::exponential
double precision::ltzsquare
double precision::cosecant

!File names
!==================================================
character(len=6)::hrfile='hr.dat'
character(len=7)::hdrfile='hdr.dat'
character(len=7)::hcqfile='hcq.dat'
character(len=7)::hdqfile='hdq.dat'
character(len=7)::ccqfile='ccq.dat'
character(len=6)::crfile='cr.dat'
character(len=7)::cdrfile='cdr.dat'
character(len=7)::scqfile='scq.dat'
character(len=7)::sdqfile='sdq.dat'
character(len=7)::cdqfile='cdq.dat'

!Initializations
!==================================================
io=0
nbiter=0
convergence=1.0d0
fileex=.false.

prec=1.0d-12
pi=4.0d0*datan(1.0d0)

dq=pi/((2.0d0**(p2)*dr))
rmax=dr*(2**(p2))

do i=0,2**p2
  k(i)=0.0d0
  r(i)=0.0d0
  q(i)=0.0d0
  cr(i)=0.0d0
  cr2(i)=0.0d0
  cq(i)=0.0d0
  cdr(i)=0.0d0
  cdr2(i)=0.0d0
  cdq(i)=0.0d0
  ccr(i)=0.0d0
  ccq(i)=0.0d0
  gamr(i)=0.0d0
  gamq(i)=0.0d0
  gamdr(i)=0.0d0
  gamdq(i)=0.0d0
  hr(i)=0.0d0
  hq(i)=0.0d0
  hdr(i)=0.0d0
  hdq(i)=0.0d0
  hcr(i)=0.0d0
  hcq(i)=0.0d0
  scq(i)=0.0d0
  sdq(i)=0.0d0
  hrref(i)=0.0d0
  crref(i)=0.0d0
  crref2(i)=0.0d0
  ccrref(i)=0.0d0
  cqref(i)=0.0d0
  gamqref(i)=0.0d0
  gamrref(i)=0.0d0

  hrorpa(i)=0.0d0
  hdrorpa(i)=0.0d0
  ssr(i)=0.0d0
  ssdr(i)=0.0d0
  hhr(i)=0.0d0
  hhdr(i)=0.0d0
  ssq(i)=0.0d0
  ssdq(i)=0.0d0
  hhq(i)=0.0d0
  hhdq(i)=0.0d0
end do

eta=density*(pi/6.0d0)
lamb1=((1.0d0+2.0d0*eta)**2)/((1.0d0-eta)**4)
lamb2=(-(1.0d0+0.5d0*eta)**2)/((1.0d0-eta)**4)

do i=0,2**p2 !Distances
  r(i)=i*dr
  if (dabs(r(i)-1.0d0) .lt. (0.01d0*dr)) disc=i
  q(i)=i*dq
end do

!Calculation of the disorder correl function according to the user's choice
if (correl .eq. 'gaus') then
  do i=0,2**p2 !K(r)
  k(i)=delta*gaussian(sigma,r(i))
  end do

else if (correl .eq. 'expl') then
  do i=0,2**p2 !K(r)
  k(i)=delta*exponential(sigma,r(i))
  end do

else if (correl .eq. 'lzsq') then
  do i=0,2**p2 !K(r)
  k(i)=delta*ltzsquare(sigma,r(i))
  end do

else if (correl .eq. 'exga') then
  do i=0,2**p2 !K(r)
  k(i)=delta*(dexp(gaussian(sigma,r(i)))-1d0)
  end do

else if (correl .eq. 'cosc') then
  do i=0,2**p2 !K(r)
  k(i)=delta*cosecant(sigma,r(i))
  end do
end if


!Test the existence of a file named cr.dat
if (cr_init .eq. 'file') then
  inquire(file='cr.dat',exist=fileex)
  if (fileex .eqv. .true.) then
    call fileman(crfile,len(crfile),11,1)
    do i=0,2**p2
      read (11,*) r(i),cr(i)
    end do
    call fileman(crfile,len(crfile),11,0)
  else
    do i=0,disc-1 !cr
      cr(i)=-lamb1-6.0d0*eta*lamb2*r(i)-0.5d0*eta*lamb1*r(i)**3 !py
    end do
    cr(disc)=(-lamb1-6.0d0*eta*lamb2*r(disc)-0.5d0*eta*lamb1*r(disc)**3)/2.0d0 !py
    do i=disc+1,2**p2
      cr(i)=0.0d0 !py
    end do
  end if

else if (cr_init .eq. 'pyev') then
  do i=0,disc-1 !cr
    cr(i)=-lamb1-6.0d0*eta*lamb2*r(i)-0.5d0*eta*lamb1*r(i)**3 !py
  end do
  cr(disc)=(-lamb1-6.0d0*eta*lamb2*r(disc)-0.5d0*eta*lamb1*r(disc)**3)/2.0d0 !py
  do i=disc+1,2**p2
    cr(i)=0.0d0 !py
  end do

else if (cr_init .eq. 'hncc') then
  do i=0,disc-1 !cr
    cr(i)=-1.0d0 !hnc
  end do
  cr(disc)=(exp(k(disc))-2.0d0)/2.0d0 !hnc
  do i=disc+1,2**p2
    cr(i)=exp(k(i))-1.0d0 !hnc
  end do

else if (cr_init .eq. 'zero') then
  do i=0,2**p2
    cr(i)=0.0d0
  end do

end if

if (cr_init .eq. 'file') then
  inquire(file='cdr.dat',exist=fileex)
  if (fileex .eqv. .true.) then
    call fileman(cdrfile,len(cdrfile),11,1)
    do i=0,2**p2
      read (11,*) r(i),cdr(i)
    end do
    call fileman(cdrfile,len(cdrfile),11,0)
  else
    do i=0,2**p2
      cdr(i)=0.0d0 !py
    end do
  end if

else if (cr_init .eq. 'hncc') then
  do i=0,2**p2
    cdr(i)=exp(k(i))-1.0d0 !hnc
  end do

else if ((cr_init .eq. 'pyev') .or. (cr_init .eq. 'zero')) then
  do i=0,2**p2 !cdr
    cdr(i)=0.0d0 !py
  end do

end if


!Calculation of the Waisman parametrization for the ORPApproximation
!==================================================
if ((closure .eq. 'orpa') .or. (closure .eq. 'oexp')) then
  xred=1.0d0/((8.0d0*eta-2.0d0*eta**2)/(1.0d0-eta)**4+1.0d0)
  cs=(1+eta+eta**2-eta**3)/((1-eta)**3)
  fw=(1.0d0-eta)*dsqrt(1.0d0/xred)
  qw=(1.0d0+2.0d0*eta)**2/(1.0d0-eta)**2
  v0=(6.0d0/4.0d0)*(cs-1.0d0)-fw**2+1.0d0
  z1=(2.0d0/(qw-fw**2))*((v0+fw**2-qw)*fw+dsqrt((v0+fw**2-qw)*v0*qw))
  sigma1=(1.0d0/(2.0d0*z1))*((z1-2.0d0)/(z1+2.0d0)+dexp(-z1))
  tau1=(1.0d0/(2.0d0*z1))*((z1**2+2.0d0*z1-4.0d0)/(4.0d0+2.0d0*z1-z1**2)+dexp(-z1))
  alpha1=((4.0d0+2.0d0*z1-z1**2)*tau1)/(2.0d0*(2.0d0+z1)*sigma1)
  v1v0=2.0d0-dsqrt(qw)-(1.0d0/(2.0d0*v0*dsqrt(qw)))*((v0+fw**2-qw)*(v0+fw**2)+0.25d0*(z1**2)*(qw-fw**2))
  k1=((2.0d0*((z1+2.0d0)**2)*(sigma1**2))/(3.0d0*eta*(z1**2)))*v0*(v1v0-alpha1)**2
end if

!Iterate!!!!
!==================================================
nbiter=0
convergence=1.0d0

do while (convergence .gt. prec)

!  write (6,*) "convergence ",convergence

  nbiter=nbiter+1

!Fourier transform of cr and cdr
!==================================================
  cq(0)=0.0d0
  cdq(0)=0.0d0

  do i=0,2**p2
    cq(0)=cq(0)+4.0d0*pi*(r(i)**2)*cr(i)*dr
    cdq(0)=cdq(0)+4.0d0*pi*(r(i)**2)*cdr(i)*dr
  end do

  call fft3s(0.0,cr(0),p2,rmax,cq(0),1,0)
  call fft3s(0.0,cdr(0),p2,rmax,cdq(0),1,0)


!Calculating the gamq and gamdq arrays
!==================================================
  do i=0,2**p2
    gamq(i)=(cq(i)-cdq(i))/(1.0d0-density*(cq(i)-cdq(i)))&
    +cdq(i)/((1.0d0-density*(cq(i)-cdq(i)))**2)-cq(i)

    gamdq(i)=cdq(i)/((1.0d0-density*(cq(i)-cdq(i)))**2)-cdq(i)
  end do

!Inverse FT of gamq and gamdq
!==================================================
  gamr(0)=0.0d0
  gamdr(0)=0.0d0

  do i=0,2**p2
    gamr(0)=gamr(0)+(1.0d0/(2.0d0*pi**2))*(q(i)**2)*gamq(i)*dq
    gamdr(0)=gamdr(0)+(1.0d0/(2.0d0*pi**2))*(q(i)**2)*gamdq(i)*dq
  end do


  call fft3s(0.0,gamq(0),p2,rmax,gamr(0),-1,0)
  call fft3s(0.0,gamdq(0),p2,rmax,gamdr(0),-1,0)

!Calculating the cr2 and the cdr2 arrays
!==================================================

  if (closure .eq. 'hncc') then

    do i=0,disc-1
      cr2(i)=-1.0d0-gamr(i)
      cdr2(i)=dexp(k(i)+gamdr(i))-1.0d0-gamdr(i)
    end do    
    cr2(disc)=0.5d0*dexp(k(disc)+gamr(disc))-1.0d0-gamr(disc)
    cdr2(disc)=dexp(k(disc)+gamdr(disc))-1.0d0-gamdr(disc)
    do i=disc+1,2**p2
      cr2(i)=dexp(k(i)+gamr(i))-1.0d0-gamr(i)
      cdr2(i)=dexp(k(i)+gamdr(i))-1.0d0-gamdr(i)
    end do

  else if ((closure .eq. 'msac') .or. (closure .eq. 'expc') .or. (closure .eq. 'mn3c')) then

    do i=0,disc-1
      cr2(i)=-1.0d0-gamr(i)
      cdr2(i)=k(i)
    end do    
    cr2(disc)=(k(disc)-1.0d0-gamr(disc))/2.0d0
    cdr2(disc)=k(disc)
    do i=disc+1,2**p2
      cr2(i)=k(i)
      cdr2(i)=k(i)
    end do

  else if (closure .eq. 'pyev') then
    do i=0,disc-1
      cr2(i)=-1.0d0-gamr(i)
      cdr2(i)=(dexp(k(i))-1.0d0)*(1.0d0+gamdr(i))
    end do    
    cr2(disc)=((dexp(k(disc))-1.0d0)*(1.0d0+gamr(disc))-1.0d0-gamr(disc))/2.0d0
    cdr2(disc)=(dexp(k(disc))-1.0d0)*(1.0d0+gamdr(disc))
    do i=disc+1,2**p2
      cr2(i)=(dexp(k(i))-1.0d0)*(1.0d0+gamr(i))
      cdr2(i)=(dexp(k(i))-1.0d0)*(1.0d0+gamdr(i))
    end do

  else if ((closure .eq. 'orpa') .or. (closure .eq. 'oexp')) then
    do i=0,disc-1
      cr2(i)=-1.0d0-gamr(i)
      cdr2(i)=k(i)
    end do
    cr2(disc)=(k(disc)-1.0d0-gamr(disc)+k1*(dexp(-z1*(r(disc)-1.0d0))/r(disc)))/2.0d0
    cdr2(disc)=k(disc)
    do i=disc+1,2**p2
      cr2(i)=k(i)+k1*(dexp(-z1*(r(i)-1.0d0))/r(i))
      cdr2(i)=k(i)
    end do

  end if


!Take 100*mix_param% of the old and mix_param it with 100*(1-mix_param)%
!of the new and put it in the cr2 -the stock-
!==================================================
  do i=0,2**p2 !cr
    cr2(i)=mix_param*cr(i)+(1.0d0-mix_param)*cr2(i)
    cdr2(i)=mix_param*cdr(i)+(1.0d0-mix_param)*cdr2(i)
  end do

!Calculate the convergence
!==================================================
  convergence=0.0d0
  do i=0,2**p2
    convergence=convergence+dabs(cr(i)-cr2(i))+dabs(cdr(i)-cdr2(i))
  end do

!Passing the values and testing of the program did not crash
!==================================================
  do i=0,2**p2
    cr(i)=cr2(i)
    cdr(i)=cdr2(i)

    if ((cr(i) .ne. cr(i)) .or. (cdr(i) .ne. cdr(i))) then
      write (6,*) "STRUCTURE CRASHED"
      stop
    end if

  end do

!Calculate hr and hdr, and some other stuff
!==================================================
  do i=0,2**p2
    hr(i)=gamr(i)+cr(i) !hr
    hdr(i)=gamdr(i)+cdr(i) !hdr
    !Since this value as well as the exp and the ref value will be needed we have to stock it somewhere
    if (closure .eq. 'mn3c') then
      hrorpa(i)=hr(i)
      hdrorpa(i)=hdr(i)
    end if
    hq(i)=gamq(i)+cq(i) !hq
    hdq(i)=gamdq(i)+cdq(i) !hdq
    hcq(i)=hq(i)-hdq(i) !hcq
    ccq(i)=cq(i)-cdq(i) !ccq
    scq(i)=1.0d0+density*hcq(i)
    sdq(i)=density*hdq(i)
  end do


end do

!Better write it here in case you want to restart a calculation
if ((closure .eq. 'expc') .or. (closure .eq. 'oexp') .or. (closure .eq. 'mn3c')) then
  call fileman(crfile,len(crfile),16,1)
  call fileman(cdrfile,len(cdrfile),17,1)
  do i=0,2**p2
    write(16,*) r(i),cr(i)
    write(17,*) r(i),cdr(i)
  end do
  call fileman(crfile,len(crfile),16,0)
  call fileman(cdrfile,len(cdrfile),17,0)
end if




!Calculation of the reference system for EXPC
!==================================================
if ((closure .eq. 'expc') .or. (closure .eq. 'mn3c')) then

  convergence=1.0d0

  do i=0,disc-1 !cr
    crref(i)=-lamb1-6.0d0*eta*lamb2*r(i)-0.5d0*eta*lamb1*r(i)**3 !py
  end do
  crref(disc)=(-lamb1-6.0d0*eta*lamb2*r(disc)-0.5d0*eta*lamb1*r(disc)**3)/2.0d0 !py
  do i=disc+1,2**p2
    crref(i)=0.0d0 !py
  end do
 

  do while (convergence .gt. prec)

!    write (6,*) "convergence expc",convergence

    cqref(0)=0.0d0
    do i=0,2**p2
      cqref(0)=cqref(0)+4.0d0*pi*(r(i)**2)*crref(i)*dr
    end do
 
    call fft3s(0.0,crref(0),p2,rmax,cqref(0),1,0)
 
    do i=0,2**p2
      gamqref(i)=cqref(i)/(1.0d0-density*cqref(i))-cqref(i)
    end do

    gamrref(0)=0.0d0
    do i=0,2**p2
      gamrref(0)=gamrref(0)+(1.0d0/(2.0d0*pi**2))*(q(i)**2)*gamqref(i)*dq
    end do
 
    call fft3s(0.0,gamqref(0),p2,rmax,gamrref(0),-1,0)

    do i=0,disc-1
      crref2(i)=-1.0d0-gamrref(i)
    end do    
    crref2(disc)=(-1.0d0-gamrref(disc))/2.0d0
    do i=disc+1,2**p2
      crref2(i)=0.0d0
    end do

    !Take 100*mix_param% of the old and mix_param it with 100*(1-mix_param)%
    !of the new and put it in the cr2 -the stock-
    !==================================================
    do i=0,2**p2 !cr
      crref2(i)=mix_param*crref(i)+(1.0d0-mix_param)*crref2(i)
    end do

    !Calculate the convergence
    !==================================================
    convergence=0.0d0
    do i=0,2**p2
      convergence=convergence+dabs(crref(i)-crref2(i))
    end do

    !Passing the values
    !==================================================
    do i=0,2**p2
      crref(i)=crref2(i)
    end do

  end do

  do i=0,2**p2
    hrref(i)=gamrref(i)+crref(i)
  end do

  call fileman('hrref.dat',9,11,1)
  do i=0,2**p2
    write(11,*) r(i),hrref(i)
  end do
  call fileman('hrref.dat',9,11,0)

  do i=0,disc-1
    hdr(i)=dexp(hdr(i))-1.0d0 !for the disconnected part, the renormalized potential is equal to hdr(i)
    hr(i)=-1.0d0
  end do
  hdr(disc)=dexp(hdr(disc))-1.0d0
  hr(disc)=(hrref(disc)+1.0d0)*dexp(2.0d0*(hr(disc)-hrref(disc)))-1.0d0 !the limit coming from the right of the value at the discontinuity
  do i=disc+1,2**p2
    hdr(i)=dexp(hdr(i))-1.0d0
    hr(i)=(hrref(i)+1.0d0)*dexp(hr(i)-hrref(i))-1.0d0
  end do

  hq(0)=0.0d0
  hdq(0)=0.0d0
  do i=0,2**p2
    hq(0)=hq(0)+4.0d0*pi*(r(i)**2)*hr(i)*dr
    hdq(0)=hdq(0)+4.0d0*pi*(r(i)**2)*hdr(i)*dr
  end do

  call fft3s(0.0,hr(0),p2,rmax,hq(0),1,0)
  call fft3s(0.0,hdr(0),p2,rmax,hdq(0),1,0)

  do i=0,2**p2
    hcq(i)=hq(i)-hdq(i)
    cq(i)=hq(i)/(1.0d0+density*hq(i))
    ccq(i)=hcq(i)/(1.0d0+density*hcq(i))
    cdq(i)=cq(i)-ccq(i)
    scq(i)=1.0d0+density*hcq(i)
    sdq(i)=density*hdq(i)
  end do

  cr(0)=0.0d0
  do i=0,2**p2
    cr(0)=cr(0)+(1.0d0/(2.0d0*pi**2))*(q(i)**2)*cq(i)*dq
  end do

  call fft3s(0.0,cq(0),p2,rmax,cr(0),-1,0)
  call fft3s(0.0,cdq(0),p2,rmax,cdr(0),-1,0)

end if











!Calculation of the reference system for EXPC
!==================================================
if (closure .eq. 'oexp') then

  convergence=1.0d0

  do i=0,disc-1 !cr
    crref(i)=0.0d0
  end do
  crref(disc)=(k1*(dexp(-z1*(r(disc)-1.0d0))/r(disc)))/2.0d0
  do i=disc+1,2**p2
    crref(i)=k1*(dexp(-z1*(r(i)-1.0d0))/r(i))
  end do
 

  do while (convergence .gt. prec)

!    write (6,*) "convergence expc",convergence

    cqref(0)=0.0d0
    do i=0,2**p2
      cqref(0)=cqref(0)+4.0d0*pi*(r(i)**2)*crref(i)*dr
    end do
 
    call fft3s(0.0,crref(0),p2,rmax,cqref(0),1,0)
 
    do i=0,2**p2
      gamqref(i)=cqref(i)/(1.0d0-density*cqref(i))-cqref(i)
    end do

    gamrref(0)=0.0d0
    do i=0,2**p2
      gamrref(0)=gamrref(0)+(1.0d0/(2.0d0*pi**2))*(q(i)**2)*gamqref(i)*dq
    end do
 
    call fft3s(0.0,gamqref(0),p2,rmax,gamrref(0),-1,0)

    do i=0,disc-1
      crref2(i)=-1.0d0-gamrref(i)
    end do    
    crref2(disc)=(-1.0d0-gamrref(disc)+k1*(dexp(-z1*(r(disc)-1.0d0))/r(disc)))/2.0d0
    do i=disc+1,2**p2
      crref2(i)=k1*(dexp(-z1*(r(i)-1.0d0))/r(i))
    end do

    !Take 100*mix_param% of the old and mix_param it with 100*(1-mix_param)%
    !of the new and put it in the cr2 -the stock-
    !==================================================
    do i=0,2**p2 !cr
      crref2(i)=mix_param*crref(i)+(1.0d0-mix_param)*crref2(i)
    end do

    !Calculate the convergence
    !==================================================
    convergence=0.0d0
    do i=0,2**p2
      convergence=convergence+dabs(crref(i)-crref2(i))
    end do

    !Passing the values
    !==================================================
    do i=0,2**p2
      crref(i)=crref2(i)
    end do

  end do

  do i=0,2**p2
    hrref(i)=gamrref(i)+crref(i)
  end do

  call fileman('hrref.dat',9,11,1)
  do i=0,2**p2
    write(11,*) r(i),hrref(i)
  end do
  call fileman('hrref.dat',9,11,0)

  do i=0,disc-1
    hdr(i)=dexp(hdr(i))-1.0d0 !for the disconnected part, the renormalized potential is equal to hdr(i)
    hr(i)=-1.0d0
  end do
  hdr(disc)=dexp(hdr(disc))-1.0d0
  hr(disc)=(hrref(disc)+1.0d0)*dexp(2.0d0*(hr(disc)-hrref(disc)))-1.0d0 !the limit coming from the right of the value at the discontinuity
  do i=disc+1,2**p2
    hdr(i)=dexp(hdr(i))-1.0d0
    hr(i)=(hrref(i)+1.0d0)*dexp(hr(i)-hrref(i))-1.0d0
  end do

  hq(0)=0.0d0
  hdq(0)=0.0d0
  do i=0,2**p2
    hq(0)=hq(0)+4.0d0*pi*(r(i)**2)*hr(i)*dr
    hdq(0)=hdq(0)+4.0d0*pi*(r(i)**2)*hdr(i)*dr
  end do

  call fft3s(0.0,hr(0),p2,rmax,hq(0),1,0)
  call fft3s(0.0,hdr(0),p2,rmax,hdq(0),1,0)

  do i=0,2**p2
    hcq(i)=hq(i)-hdq(i)
    cq(i)=hq(i)/(1.0d0+density*hq(i))
    ccq(i)=hcq(i)/(1.0d0+density*hcq(i))
    cdq(i)=cq(i)-ccq(i)
    scq(i)=1.0d0+density*hcq(i)
    sdq(i)=density*hdq(i)
  end do

  cr(0)=0.0d0
  do i=0,2**p2
    cr(0)=cr(0)+(1.0d0/(2.0d0*pi**2))*(q(i)**2)*cq(i)*dq
  end do

  call fft3s(0.0,cq(0),p2,rmax,cr(0),-1,0)
  call fft3s(0.0,cdq(0),p2,rmax,cdr(0),-1,0)

end if







if (closure .eq. 'mn3c') then

  do i=0,2**p2
    ssr(i)=(hrref(i)+1.0d0)*dexp(hrorpa(i)-hrref(i))-1.0d0-hrorpa(i)+hrref(i)
    ssdr(i)=dexp(hdrorpa(i))-1.0d0-hdrorpa(i)

    hhr(i)=hrorpa(i)
    hhdr(i)=hdrorpa(i)
  end do

  call fft3s(0.0,ssq(0),p2,rmax,ssq(0),1,0)
  call fft3s(0.0,ssdq(0),p2,rmax,ssdq(0),1,0)

  call fft3s(0.0,hhq(0),p2,rmax,hhq(0),1,0)
  call fft3s(0.0,hhdq(0),p2,rmax,hhdq(0),1,0)

  do i=0,2**p2
    hq(i)=(hq(i)+1.0d0)*(1.0d0+2.0d0*density*ssq(i)*hhq(i)+density*ssq(i)**2)-1.0d0
    hdq(i)=(hdq(i)+1.0d0)*(1.0d0+2.0d0*density*ssdq(i)*hhdq(i)+density*ssdq(i)**2)-1.0d0
  end do

  call fft3s(0.0,hq(0),p2,rmax,hr(0),-1,0)
  call fft3s(0.0,hdq(0),p2,rmax,hdr(0),-1,0)

  do i=0,2**p2
    hcq(i)=hq(i)-hdq(i)
    cq(i)=hq(i)/(1.0d0+density*hq(i))
    ccq(i)=hcq(i)/(1.0d0+density*hcq(i))
    cdq(i)=cq(i)-ccq(i)
    scq(i)=1.0d0+density*hcq(i)
    sdq(i)=density*hdq(i)
  end do

  cr(0)=0.0d0
  do i=0,2**p2
    cr(0)=cr(0)+(1.0d0/(2.0d0*pi**2))*(q(i)**2)*cq(i)*dq
  end do

  call fft3s(0.0,cq(0),p2,rmax,cr(0),-1,0)
  call fft3s(0.0,cdq(0),p2,rmax,cdr(0),-1,0)

end if




!Print the stuff
!==================================================
call fileman(hrfile,len(hrfile),11,1)
call fileman(hdrfile,len(hdrfile),12,1)
call fileman(hcqfile,len(hcqfile),13,1)
call fileman(hdqfile,len(hdqfile),14,1)
call fileman(ccqfile,len(ccqfile),15,1)
if (closure .eq. 'expc') then
  call fileman('cr_expc.dat',11,16,1)
  call fileman('cdr_expc.dat',12,17,1)
else if (closure .eq. 'oexp') then
  call fileman('cr_oexp.dat',11,16,1)
  call fileman('cdr_oexp.dat',12,17,1)
else if (closure .eq. 'mn3c') then
  call fileman('cr_mn3c.dat',11,16,1)
  call fileman('cdr_mn3c.dat',12,17,1)
else
  call fileman(crfile,len(crfile),16,1)
  call fileman(cdrfile,len(cdrfile),17,1)
end if
call fileman(scqfile,len(scqfile),18,1)
call fileman(sdqfile,len(sdqfile),19,1)
call fileman(cdqfile,len(cdqfile),20,1)

do i=0,2**p2
  write(11,*) r(i),hr(i)
  write(12,*) r(i),hdr(i)
  write(13,*) q(i),hcq(i)
  write(14,*) q(i),hdq(i)
  write(15,*) q(i),ccq(i)
  write(16,*) r(i),cr(i)
  write(17,*) r(i),cdr(i)
  write(18,*) q(i),scq(i)
  write(19,*) q(i),sdq(i)
  write(20,*) q(i),cdq(i)
end do

call fileman(hrfile,len(hrfile),11,0)
call fileman(hdrfile,len(hdrfile),12,0)
call fileman(hcqfile,len(hcqfile),13,0)
call fileman(hdqfile,len(hdqfile),14,0)
call fileman(ccqfile,len(ccqfile),15,0)
if (closure .eq. 'expc') then
  call fileman('cr_expc.dat',11,16,0)
  call fileman('cdr_expc.dat',12,17,0)
else if (closure .eq. 'oexp') then
  call fileman('cr_oexp.dat',11,16,0)
  call fileman('cdr_oexp.dat',12,17,0)
else if (closure .eq. 'mn3c') then
  call fileman('cr_mn3c.dat',11,16,0)
  call fileman('cdr_mn3c.dat',12,17,0)
else
  call fileman(crfile,len(crfile),16,0)
  call fileman(cdrfile,len(cdrfile),17,0)
end if
call fileman(scqfile,len(scqfile),18,0)
call fileman(sdqfile,len(sdqfile),19,0)
call fileman(cdqfile,len(cdqfile),20,0)

end subroutine
