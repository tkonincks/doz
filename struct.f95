subroutine struct (density,delta,sigma,closure,mix_param,cr_init)

implicit none

double precision::density
double precision::delta
double precision::sigma
character(len=4)::closure,cr_init
double precision::mix_param


integer::disc
integer::io

!Iteration variables
!==================================================
integer::i
integer::nbiter

!Constants
!==================================================
double precision,parameter::prec=1.0d-12
double precision,parameter::pi=4.0d0*datan(1.0d0)
integer,parameter::p2=12

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

!Arrays
!==================================================
double precision,dimension(0:2**p2)::r
double precision,dimension(0:2**p2)::q
double precision,parameter::dr=0.01d0
double precision,parameter::dq=pi/((2.0d0**p2)*dr)

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

!These are for the EXPC closure
double precision,dimension(0:2**p2)::hrref
double precision,dimension(0:2**p2)::hdrref
double precision,dimension(0:2**p2)::crref
double precision,dimension(0:2**p2)::cdrref
double precision,dimension(0:2**p2)::ccrref
double precision,dimension(0:2**p2)::cqref
double precision,dimension(0:2**p2)::cdqref
double precision,dimension(0:2**p2)::gamqref
double precision,dimension(0:2**p2)::gamdqref
double precision,dimension(0:2**p2)::gamrref
double precision,dimension(0:2**p2)::gamdrref
double precision,dimension(0:2**p2)::krenorm !the renormalized potential


double precision,dimension(0:2**p2)::conv1 !Array that contains the discrete convergence
double precision::convergence !Used to calculate the convergence of the calcultaion

!Functions
!==================================================
double precision::gaussian

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

!Initializations
!==================================================
io=0
nbiter=0
eta=0.0d0
lamb1=0.0d0
lamb2=0.0d0
convergence=1.0d0
fileex=.false.

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
  hdrref(i)=0.0d0
  crref(i)=0.0d0
  cdrref(i)=0.0d0
  ccrref(i)=0.0d0
  cqref(i)=0.0d0
  cdqref(i)=0.0d0
  gamqref(i)=0.0d0
  gamdqref(i)=0.0d0
  gamrref(i)=0.0d0
  gamdrref(i)=0.0d0
  krenorm(i)=0.0d0 
  conv1(i)=0.0d0
end do

eta=density*(4.0d0/3.0d0)*pi
lamb1=((1.0d0+2.0d0*eta)**2)/((1.0d0-eta)**4)
lamb2=(-(1.0d0+0.5d0*eta)**2)/((1.0d0-eta)**4)

do i=0,2**p2 !Distances
  r(i)=i*dr
  if (dabs(r(i)-1.0d0) .lt. (0.01d0*dr)) disc=i
  q(i)=i*dq
end do

do i=0,2**p2 !K(r)
  k(i)=delta*gaussian(sigma,r(i))
end do

!Test the existence of a file named cr.dat
if (cr_init .eq. 'file') then
  call fileman(crfile,len(crfile),11,1)
  do i=0,2**p2
    read (11,*) r(i),cr(i)
  end do
  call fileman(crfile,len(crfile),11,0)

else if (cr_init .eq. 'pyev') then
  do i=0,disc-1 !cr
    cr(i)=-lamb1-6.0d0*eta*lamb2*r(i)-0.5d0*eta*lamb1*r(i)**3 !py
  end do
  cr(disc)=(-lamb1-6.0d0*eta*lamb2*r(i)-0.5d0*eta*lamb1*r(i)**3)/2.0d0 !py
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
  call fileman(cdrfile,len(cdrfile),11,1)
  do i=0,2**p2
    read (11,*) r(i),cdr(i)
  end do
  call fileman(cdrfile,len(cdrfile),11,0)

else if (cr_init .eq. 'hncc') then
  do i=0,2**p2
    cdr(i)=exp(k(i))-1.0d0 !hnc
  end do

else if ((cr_init .eq. 'pyev') .or. (cr_init .eq. 'zero')) then
  do i=0,2**p2 !cdr
    cdr(i)=0.0d0 !py
  end do

end if

!Iterate!!!!
!==================================================
nbiter=0
convergence=1.0d0

do while (convergence .gt. prec)

  nbiter=nbiter+1

!Fourrier transform of cr and cdr
!==================================================
  cq(0)=0.0d0
  cdq(0)=0.0d0

  do i=0,2**p2
    cq(0)=cq(0)+4.0d0*pi*(r(i)**2)*cr(i)*dr
    cdq(0)=cdq(0)+4.0d0*pi*(r(i)**2)*cdr(i)*dr
  end do

  call fft3s(0.0,cr(0),p2,r(2**p2),cq(0),1,0)
  call fft3s(0.0,cdr(0),p2,r(2**p2),cdq(0),1,0)


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

  do i=0,2**p2 !i=0 or i=1?
    gamr(0)=gamr(0)+(1.0d0/(2.0d0*pi**2))*(q(i)**2)*gamq(i)*dq
    gamdr(0)=gamdr(0)+(1.0d0/(2.0d0*pi**2))*(q(i)**2)*gamdq(i)*dq
  end do


  call fft3s(0.0,gamq(0),p2,r(2**p2),gamr(0),-1,0)
  call fft3s(0.0,gamdq(0),p2,r(2**p2),gamdr(0),-1,0)

!Calculating the cr2 and the cdr2 arrays
!==================================================

  if (closure .eq. 'hncc') then

    do i=0,disc-1
      cr2(i)=-1.0d0-gamr(i)
      cdr2(i)=dexp(k(i)+gamdr(i))-1.0d0-gamdr(i)
    end do    
    cr2(disc)=(dexp(k(disc)+gamr(disc))-2.0d0-2.0d0*gamr(disc))/2.0d0
    cdr2(disc)=dexp(k(disc)+gamdr(disc))-1.0d0-gamdr(disc)
    do i=disc+1,2**p2
      cr2(i)=dexp(k(i)+gamr(i))-1.0d0-gamr(i)
      cdr2(i)=dexp(k(i)+gamdr(i))-1.0d0-gamdr(i)
    end do

  else if ((closure .eq. 'msac') .or. (closure .eq. 'expc')) then

    do i=0,disc-1
      cr2(i)=-1.0d0-gamr(i)
      cdr2(i)=k(i)
    end do    
    cr2(disc)=(k(disc+1)-1.0d0-gamr(disc+1))/2.0d0
    cdr2(disc)=k(disc)
    do i=disc+1,2**p2
      cr2(i)=k(i)
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

!Calculate hr and hdr
!==================================================

  if (closure .eq. 'hncc') then

    do i=0,2**p2
      hr(i)=gamr(i)+cr(i) !hr
      hdr(i)=gamdr(i)+cdr(i) !hdr
      hq(i)=gamq(i)+cq(i) !hq
      hdq(i)=gamdq(i)+cdq(i) !hdq
      hcq(i)=hq(i)-hdq(i) !hcq
      ccq(i)=cq(i)-cdq(i) !ccq
      scq(i)=1.0d0+density*hcq(i)
      sdq(i)=density*hdq(i)
    end do

  else if ((closure .eq. 'msac') .or. (closure .eq. 'expc')) then

    do i=0,disc-1
      hr(i)=-1.0d0
    end do
    hr(disc)=(gamr(disc+1)+cr(disc+1)-1.0d0)/2.0d0
    do i=disc+1,2**p2
      hr(i)=gamr(i)+cr(i) !hr
    end do
 
    do i=0,2**p2
      hdr(i)=gamdr(i)+cdr(i) !hdr
      hq(i)=gamq(i)+cq(i) !hq
      hdq(i)=gamdq(i)+cdq(i) !hdq
      hcq(i)=hq(i)-hdq(i) !hcq
      ccq(i)=cq(i)-cdq(i) !ccq
      scq(i)=1.0d0+density*hcq(i)
      sdq(i)=density*hdq(i)
    end do

  end if


end do

!If the density is zero, take the analytical solution of hdr
!==================================================
if (density .eq. 0.0d0) then

  if (closure .eq. 'hncc') then

    do i=0,2**p2
      hdr(i)=exp(k(i))-1.0d0
    end do

  else if ((closure .eq. 'msac') .or. (closure .eq. 'expc')) then

    do i=0,disc-1
      hdr(i)=k(i)
      hr(i)=-1.0d0
    end do
      hdr(disc)=k(disc)
      hr(disc)=(k(disc+1)-1.0d0)/2.0d0
    do i=disc+1,2**p2
      hdr(i)=k(i)
      hr(i)=k(i)
    end do

  end if

end if

call fft3s(0.0,hdr(0),p2,r(2**p2),hdq(0),1,0)


!Calculation of the reference system for EXPC
!==================================================
if (closure .eq. 'expc') then

  do i=0,disc-1 !cr
    crref(i)=-lamb1-6.0d0*eta*lamb2*r(i)-0.5d0*eta*lamb1*r(i)**3 !py
  end do
  crref(disc)=(-lamb1-6.0d0*eta*lamb2*r(disc+1)-0.5d0*eta*lamb1*r(disc+1)**3)/2.0d0 !py
  do i=disc+1,2**p2
    crref(i)=0.0d0 !py
  end do

  do i=0,2**p2
    cqref(0)=cqref(0)+4.0d0*pi*(r(i)**2)*crref(i)*dr
    cdqref(0)=cdqref(0)+4.0d0*pi*(r(i)**2)*cdrref(i)*dr
  end do

  call fft3s(0.0,crref(0),p2,r(2**p2),cqref(0),1,0)
  call fft3s(0.0,cdrref(0),p2,r(2**p2),cdqref(0),1,0)

  do i=0,2**p2
    gamqref(i)=(cqref(i)-cdqref(i))/(1.0d0-density*(cqref(i)-cdqref(i)))&
    +cdqref(i)/((1.0d0-density*(cqref(i)-cdqref(i)))**2)-cqref(i)
    gamdqref(i)=cdqref(i)/((1.0d0-density*(cqref(i)-cdqref(i)))**2)-cdqref(i)
  end do

  do i=0,2**p2
    gamrref(0)=gamrref(0)+(1.0d0/(2.0d0*pi**2))*(q(i)**2)*gamqref(i)*dq
    gamdrref(0)=gamdrref(0)+(1.0d0/(2.0d0*pi**2))*(q(i)**2)*gamdqref(i)*dq
  end do

  call fft3s(0.0,gamqref(0),p2,r(2**p2),gamrref(0),-1,0)
  call fft3s(0.0,gamdqref(0),p2,r(2**p2),gamdrref(0),-1,0)

  do i=0,disc-1
    hrref(i)=-1.0d0
    hdrref(i)=gamdrref(i)+cdrref(i)
    krenorm(i)=hr(i)-hrref(i) !calculation of the renomalized potential
  end do
  hrref(disc)=(gamrref(disc+1)+crref(disc+1)-1.0d0)/2.0d0
  hdrref(disc)=gamdrref(disc)+cdrref(disc)
  krenorm(disc)=hr(disc)-hrref(disc) !calculation of the renomalized potential
  do i=disc+1,2**p2
    hrref(i)=gamrref(i)+crref(i) !hr
    hdrref(i)=gamdrref(i)+cdrref(i) !hdr
    krenorm(i)=hr(i)-hrref(i) !calculation of the renomalized potential
  end do

  do i=0,disc-1
    hdr(i)=dexp(hdr(i))-1.0d0 !for the disconnected part, the renormalized potential is equal to hdr(i)
    hr(i)=(hrref(i)+1.0d0)*dexp(krenorm(i))-1.0d0
  end do
  hdr(disc)=dexp(hdr(disc))-1.0d0
  hr(disc)=0.5d0*(2.0d0*hrref(disc)+1.0d0)*dexp(2.0d0*(hr(disc)-hrref(disc)))-1.0d0
  do i=disc+1,2**p2
    hdr(i)=dexp(hdr(i))-1.0d0
    hr(i)=(hrref(i)+1.0d0)*dexp(krenorm(i))-1.0d0
  end do

  hq(0)=0.0d0
  hdq(0)=0.0d0

  do i=0,2**p2
    hq(0)=hq(0)+4.0d0*pi*(r(i)**2)*hr(i)*dr
    hdq(0)=hdq(0)+4.0d0*pi*(r(i)**2)*hdr(i)*dr
  end do

  call fft3s(0.0,hr(0),p2,r(2**p2),hq(0),1,0)
  call fft3s(0.0,hdr(0),p2,r(2**p2),hdq(0),1,0)

  do i=0,2**p2
    hcq(i)=hq(i)-hdq(i)
    scq(i)=1.0d0+density*hcq(i)
    sdq(i)=density*hdq(i)
  end do



end if

!Print the stuff
!==================================================
call fileman(hrfile,len(hrfile),11,1)
call fileman(hdrfile,len(hdrfile),12,1)
call fileman(hcqfile,len(hcqfile),13,1)
call fileman(hdqfile,len(hdqfile),14,1)
call fileman(ccqfile,len(ccqfile),15,1)
call fileman(crfile,len(crfile),16,1)
call fileman(cdrfile,len(cdrfile),17,1)
call fileman(scqfile,len(scqfile),18,1)
call fileman(sdqfile,len(sdqfile),19,1)

call fileman('hq.dat',6,20,1)

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

  write(20,*) q(i),hq(i)
end do

call fileman(hrfile,len(hrfile),11,0)
call fileman(hdrfile,len(hdrfile),12,0)
call fileman(hcqfile,len(hcqfile),13,0)
call fileman(hdqfile,len(hdqfile),14,0)
call fileman(ccqfile,len(ccqfile),15,0)
call fileman(crfile,len(crfile),16,0)
call fileman(cdrfile,len(cdrfile),17,0)
call fileman(scqfile,len(scqfile),18,0)
call fileman(sdqfile,len(sdqfile),19,0)

call fileman('hq.dat',6,20,0)


end subroutine
