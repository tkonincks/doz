program main

implicit none

integer::iter

!Inputs of the input.doz file
!==================================================
character(len=4)::trans_mode,closure,calc_mode,var_param,cr_init,fq_init
double precision::density,delta,sigma
double precision::var_prec

double precision::var_liq
double precision::var_glas

double precision::dens_lo
double precision::delt_li_lo
double precision::delt_gl_lo
double precision::delt_lo
double precision::lamb_lo
          
double precision::dens_hi
double precision::delt_li_hi
double precision::delt_gl_hi
double precision::delt_hi
double precision::lamb_hi

logical::fast

double precision::mix_param

logical::inflex

!Ouputs of the run subroutine
!==================================================
double precision::eigenvalue,lambda
double precision::eigenvalue_inflex,lambda_inflex

!For the dichotomy, calculation of the precision and stuff
!==================================================
double precision::prec_eigen_cont
double precision::prec_eigen_disc
double precision::prec_fq
double precision::conv_fq
double precision::delt_liq
double precision::delt_glas
double precision::dens_liq
double precision::dens_glas

!For the dichotomy around lambda
!==================================================

double precision::prec_lamb

double precision::var_incr !Incrementation of the var_param variable to find the transition in the dichotomy mode

!Array and stuff to calculate the convergence in the case of the discontinuous transition
!==================================================
integer::i=0
integer,parameter::a_size=299
double precision,dimension(0:a_size)::q
double precision,dimension(0:a_size)::fq

!To recover the time and date
!==================================================
character(len=8)::d
character(len=10)::t

!Functions
!==================================================
integer::witb

!For the calculation time
!==================================================
double precision::calc_end,iter_start,iter_end










!Initializations
!==================================================
trans_mode='xxxx'
calc_mode='xxxx'
var_param='xxxx'
closure='xxxx'

density=0.0d0
delta=0.0d0
sigma=0.0

eigenvalue=0.0d0
lambda=0.0d0
eigenvalue_inflex=0.0d0
lambda_inflex=0.0d0
inflex=.false.

var_liq=0.0d0
var_glas=0.0d0

conv_fq=0.0d0
delt_liq=0.0d0
delt_glas=0.0d0
dens_liq=0.0d0
dens_glas=0.0d0

dens_lo=0.0d0
delt_li_lo=0.0d0
delt_gl_lo=0.0d0
delt_lo=0.0d0
lamb_lo=0.0d0
                 
dens_hi=0.0d0
delt_li_hi=0.0d0
delt_gl_hi=0.0d0
delt_hi=0.0d0
lamb_hi=0.0d0
                 
var_incr=1.05d0 !Has to be higher than 1.0d0
var_prec=0.0d0
prec_eigen_cont=1.0d-6
prec_eigen_disc=1.0d-4
prec_fq=1.0d-6
prec_lamb=1.0d-6

fast=.false.

iter=0

do i=0,a_size
  q(i)=0.0d0
  fq(i)=0.0d0
end do











!Start the calculation
!==================================================
call date_and_time(DATE=d,TIME=t)
write (6,*) '       ..                              '
write (6,*) '     dF                                '
write (6,*) '    ''88bu.             u.        ..    '
write (6,*) '    ''*88888bu    ...ue888b     .@88i   '
write (6,*) '      ^"*8888N   888R Y888r   ""%888>  '
write (6,*) '     beWE "888L  888R I888>     ''88%   '
write (6,*) '     888E  888E  888R I888>   ..dILr~` '
write (6,*) '     888E  888E  888R I888>  ''".-%88b  '
write (6,*) '     888E  888F u8888cJ888    @  ''888k '
write (6,*) '    .888N..888   "*888*P"    8F   8888 '
write (6,*) '     `"888*""      ''Y"      ''8    8888 '
write (6,*) '        ""                  ''8    888F '
write (6,*) '                             %k  <88F  '
write (6,*) '                              "+:*%`   '
write (6,*) '         VERSION 05012015         '
write (6,*) ''
write (6,'(a24,a8,a,a10)') "CALCULATION LAUNCHED ON ",d," ",t
                                   









!Call the read subroutine to collect all the data
call read (trans_mode,closure,calc_mode,var_param,density,delta,sigma,var_incr&
,var_prec,var_liq,var_glas,dens_lo,delt_li_lo,delt_gl_lo,dens_hi,delt_li_hi&
,delt_gl_hi,fast,mix_param,cr_init,fq_init)

prec_eigen_cont=var_prec
prec_eigen_disc=var_prec

if ((calc_mode .eq. 'rest') .and. (var_param .eq. 'delt')) then
  delt_liq=var_liq
  delt_glas=var_glas
else if ((calc_mode .eq. 'rest') .and. (var_param .eq. 'dens')) then
  dens_liq=var_liq
  dens_glas=var_glas
end if

!if (trans_mode .ne. 'disc') .and. ((calc_mode .ne. 'dich') .or. (calc_mode .ne. 'rest')) fast=.false.

800 format (a53)
900 format (a13,a4)
901 format (a13,f25.16)
902 format (a13,es25.16)
903 format (a13,l8)








write (6,*) ""
write (6,800) "           PARAMETERS OF THE CALCULATION             "
write (6,800) "-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-"
write (6,900) 'trans_mode = ',trans_mode
write (6,900) 'closure    = ',closure
write (6,901) 'mix_param  = ',mix_param
write (6,900) 'calc_mode  = ',calc_mode
if (fast .eqv. .true.) write (6,'(a19)') "Fast Dichotomy mode"

!Print the parameters
if (var_param .eq. 'lamb') then
  write (6,900) 'var_param  = ',var_param
  write (6,901) 'var_prec   = ',var_prec
  write (6,901) 'dens_lo    = ',dens_lo
  write (6,901) 'delt_li_lo = ',delt_li_lo
  write (6,901) 'delt_gl_lo = ',delt_gl_lo
  write (6,901) 'dens_hi    = ',dens_hi
  write (6,901) 'delt_li_hi = ',delt_li_hi
  write (6,901) 'delt_gl_hi = ',delt_gl_hi
  write (6,901) 'sigma      = ',sigma

else if (calc_mode .eq. 'sing') then
  write (6,901) 'density    = ',density
  write (6,901) 'delta      = ',delta
  write (6,901) 'sigma      = ',sigma

else
  write (6,900) 'var_param  = ',var_param
  write (6,901) 'var_incr   = ',var_incr
  write (6,901) 'var_prec   = ',var_prec
  write (6,901) 'density    = ',density
  write (6,901) 'delta      = ',delta
  write (6,901) 'sigma      = ',sigma

end if

write (6,*) ""
write (6,900) 'cr_init    = ',cr_init
write (6,900) 'fq_init    = ',fq_init






!SINGLE POINT
!=================================================
if (calc_mode .eq. 'sing') then
    write (6,*) ""
    write (6,800) "                     RESULTS                         "
    write (6,800) "-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-"
    call struct (density,delta,sigma,closure,mix_param,cr_init)
    call eigen (trans_mode,density,eigenvalue,lambda,eigenvalue_inflex,lambda_inflex,inflex,fast,fq_init)
    write (6,*) ""
    write (6,901) "eigenvalue = ",eigenvalue
    write (6,901) "lambda     = ",lambda
    if (trans_mode .eq. 'disc') then
      write (6,903) "inflexion  = ",inflex
      if (inflex .eqv. .true.) then
        write (6,901) "eigen_infl = ",eigenvalue_inflex
        write (6,901) "lamb_infl  = ",lambda_inflex
      end if
      write (6,901) "conv_fq    = ",fq(witb(fq,a_size))
    end if
  call fileman('final_res',9,11,1)
  write (11,'(f22.16,f22.16,f22.16,f22.16)') density,density,delta,delta,lambda,eigenvalue
  call fileman('final_res',9,11,0)











!DICH DELT CONT/LOCA
!=================================================
else if ((calc_mode .eq. 'dich') .and. (var_param .eq. 'delt') .and. ((trans_mode .eq. 'cont') .or. (trans_mode .eq. 'loca'))) then
  write (6,*) ""
  write (6,800) "                  FIRST RESULTS                      "
  write (6,800) "-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-"
  call struct (density,delta,sigma,closure,mix_param,cr_init)
  call eigen (trans_mode,density,eigenvalue,lambda,eigenvalue_inflex,lambda_inflex,inflex,fast,fq_init)
  write (6,*) ""
  write (6,901) "density    = ",density
  write (6,901) "delta      = ",delta
  write (6,901) "eigenvalue = ",eigenvalue
  write (6,901) "lambda     = ",lambda
  if (eigenvalue .lt. 1.0d0) then
    do while (eigenvalue .lt. 1.0d0)
      delt_liq=delta
      delta=delta*var_incr
      write (6,*) ""
      write (6,800) "          Looking for the transition point           "
      write (6,800) "-----------------------------------------------------"
      call struct (density,delta,sigma,closure,mix_param,cr_init)
      call eigen (trans_mode,density,eigenvalue,lambda,eigenvalue_inflex,lambda_inflex,inflex,fast,fq_init)
      write (6,*) ""
      write (6,901) "delt_liq   = ",delt_liq
      write (6,901) "delt_glas  = ",delt_glas
      write (6,*) ""
      write (6,901) "density    = ",density
      write (6,901) "eigenvalue = ",eigenvalue
      write (6,901) "lambda     = ",lambda
    end do
    delt_glas=delta
  else if (eigenvalue .gt. 1.0d0) then
    do while (eigenvalue .gt. 1.0d0)
      delt_glas=delta
      delta=delta/var_incr
      write (6,*) ""
      write (6,800) "          Looking for the transition point           "
      write (6,800) "-----------------------------------------------------"
      call struct (density,delta,sigma,closure,mix_param,cr_init)
      call eigen (trans_mode,density,eigenvalue,lambda,eigenvalue_inflex,lambda_inflex,inflex,fast,fq_init)
      write (6,*) ""
      write (6,901) "delt_liq   = ",delt_liq
      write (6,901) "delt_glas  = ",delt_glas
      write (6,901) ""
      write (6,901) "density    = ",density
      write (6,901) "eigenvalue = ",eigenvalue
      write (6,901) "lambda     = ",lambda
    end do
    delt_liq=delta
  end if
  write (6,*) ""
  write (6,800) "-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-"
  write (6,800) "||           ENTERING THE DICHOTOMY                ||"
  write (6,800) "-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-"
  write (6,*) ""
  write (6,800) "STARTING PARAMETERS                                  "
  write (6,*) ""
  write (6,901) "delt_liq   = ",delt_liq
  write (6,901) "delt_glas  = ",delt_glas
  write (6,901) ""
  write (6,901) "density    = ",density
  write (6,901) "eigenvalue = ",eigenvalue
  write (6,901) "lambda     = ",lambda
  do 
    call cpu_time(iter_start)
    iter=iter+1
    delta=(delt_liq+delt_glas)/2.0d0
    write (6,*) ""
    write (6,*) "                 ITERATION ",iter
    write (6,800) "-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-"
    call struct (density,delta,sigma,closure,mix_param,cr_init)
    call eigen (trans_mode,density,eigenvalue,lambda,eigenvalue_inflex,lambda_inflex,inflex,fast,fq_init)
    if (eigenvalue .lt. 1.0d0) then
      delt_liq=delta
    else if (eigenvalue .gt. 1.0d0) then
      delt_glas=delta
    end if
    call cpu_time(iter_end)
    write (6,*) ""
    write (6,901) "delt_liq   = ",delt_liq
    write (6,901) "delt_glas  = ",delt_glas
    write (6,*) ""
    write (6,901) "density    = ",density
    write (6,901) "eigenvalue = ",eigenvalue
    write (6,901) "lambda     = ",lambda
    write (6,901) ""
    write (6,*) "ITERATION TIME ", iter_end-iter_start,"s"
    if (dabs(eigenvalue-1.0d0) .le. prec_eigen_cont) then
      exit
    end if
  end do
  write (6,*) ""
  write (6,800) "               CONVERGENCE REACHED                   "
  write (6,800) "-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-"
  write (6,901) "delt_liq   =  ",delt_liq
  write (6,901) "delt_glas  = ",delt_glas
  write (6,*) ""
  write (6,901) "density    = ",density
  write (6,901) "eigenvalue = ",eigenvalue
  write (6,901) "lambda     = ",lambda
  call fileman('final_res',9,11,1)
  write (11,'(f22.16,f22.16,f22.16,f22.16)') density,density,delt_liq,delt_glas,lambda,eigenvalue
  call fileman('final_res',9,11,0)















!DICH DELT DISC
!=================================================
else if ((calc_mode .eq. 'dich') .and. (var_param .eq. 'delt') .and. (trans_mode .eq. 'disc') .and. (fast .eqv. .false.)) then
  write (6,*) ""
  write (6,800) "                  FIRST RESULTS                      "
  write (6,800) "-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-"
  call struct (density,delta,sigma,closure,mix_param,cr_init)
  call eigen (trans_mode,density,eigenvalue,lambda,eigenvalue_inflex,lambda_inflex,inflex,fast,fq_init)
  call fileman('fq.dat',6,11,1)
  do i=0,a_size
    read (11,*) q(i),fq(i)
  end do
  call fileman('fq.dat',6,11,0)
  conv_fq=fq(witb(fq,a_size))
  write (6,*) ""
  write (6,901) "density    = ",density
  write (6,901) "delta      = ",delta
  write (6,901) "eigenvalue = ",eigenvalue
  write (6,901) "lambda     = ",lambda
  write (6,903) "inflexion  = ",inflex
  if (inflex .eqv. .true.) then
    write (6,901) "eigen_infl = ",eigenvalue_inflex
    write (6,901) "lamb_infl  = ",lambda_inflex
  end if
  write (6,901) "conv_fq    = ",conv_fq
  if (conv_fq .lt. prec_fq) then
    do while (conv_fq .lt. prec_fq)
      delt_liq=delta
      delta=delta*var_incr
      write (6,*) ""
      write (6,800) "          Looking for the transition point           "
      write (6,800) "-----------------------------------------------------"
      call struct (density,delta,sigma,closure,mix_param,cr_init)
      call eigen (trans_mode,density,eigenvalue,lambda,eigenvalue_inflex,lambda_inflex,inflex,fast,fq_init)
      call fileman('fq.dat',6,11,1)
      do i=0,a_size
        read (11,*) q(i),fq(i)
      end do
      call fileman('fq.dat',6,11,0)
      conv_fq=fq(witb(fq,a_size))
      write (6,*) ""
      write (6,901) "delt_liq   = ",delt_liq
      write (6,901) "delt_glas  = ",delt_glas
      write (6,*) ""
      write (6,901) "density    = ",density
      write (6,901) "eigenvalue = ",eigenvalue
      write (6,901) "lambda     = ",lambda
      write (6,903) "inflexion  = ",inflex
      if (inflex .eqv. .true.) then
        write (6,901) "eigen_infl = ",eigenvalue_inflex
        write (6,901) "lamb_infl  = ",lambda_inflex
      end if
      write (6,901) "conv_fq    = ",conv_fq
    end do
    delt_glas=delta
  else if (conv_fq .gt. prec_fq) then
    do while (conv_fq .gt. prec_fq)
      delt_glas=delta
      delta=delta/var_incr
      write (6,*) ""
      write (6,800) "          Looking for the transition point           "
      write (6,800) "-----------------------------------------------------"
      call struct (density,delta,sigma,closure,mix_param,cr_init)
      call eigen (trans_mode,density,eigenvalue,lambda,eigenvalue_inflex,lambda_inflex,inflex,fast,fq_init)
      call fileman('fq.dat',6,11,1)
      do i=0,a_size
        read (11,*) q(i),fq(i)
      end do
      call fileman('fq.dat',6,11,0)
      conv_fq=fq(witb(fq,a_size))
      write (6,*) ""
      write (6,901) "delt_liq   = ",delt_liq
      write (6,901) "delt_glas  = ",delt_glas
      write (6,901) ""
      write (6,901) "density    = ",density
      write (6,901) "eigenvalue = ",eigenvalue
      write (6,901) "lambda     = ",lambda
      write (6,903) "inflexion  = ",inflex
      if (inflex .eqv..true.) then
        write (6,901) "eigen_infl = ",eigenvalue_inflex
        write (6,901) "lamb_infl  = ",lambda_inflex
      end if
      write (6,901) "conv_fq    = ",conv_fq
    end do
    delt_liq=delta
  end if
  write (6,*) ""
  write (6,800) "-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-"
  write (6,800) "||           ENTERING THE DICHOTOMY                ||"
  write (6,800) "-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-"
  write (6,*) ""
  write (6,800) "STARTING PARAMETERS                                  "
  write (6,*) ""
  write (6,901) "delt_liq   = ",delt_liq
  write (6,901) "delt_glas  = ",delt_glas
  write (6,901) ""
  write (6,901) "density    = ",density
  write (6,901) "eigenvalue = ",eigenvalue
  write (6,901) "lambda     = ",lambda
  write (6,903) "inflexion  = ",inflex
  if (inflex .eqv. .true.) then
    write (6,901) "eigen_infl = ",eigenvalue_inflex
    write (6,901) "lamb_infl  = ",lambda_inflex
  end if
  write (6,901) "conv_fq    = ",conv_fq
  do 
    call cpu_time(iter_start)
    iter=iter+1
    delta=(delt_liq+delt_glas)/2.0d0
    write (6,*) ""
    write (6,*) "                 ITERATION ",iter
    write (6,800) "-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-"
    call struct (density,delta,sigma,closure,mix_param,cr_init)
    call eigen (trans_mode,density,eigenvalue,lambda,eigenvalue_inflex,lambda_inflex,inflex,fast,fq_init)
    call fileman('fq.dat',6,11,1)
    do i=0,a_size
      read (11,*) q(i),fq(i)
    end do
    call fileman('fq.dat',6,11,0)
    conv_fq=fq(witb(fq,a_size))
    if (conv_fq .lt. prec_fq) then
      delt_liq=delta
    else if (conv_fq .gt. prec_fq) then
      delt_glas=delta
    end if
    call cpu_time(iter_end)
    write (6,*) ""
    write (6,901) "delt_liq   = ",delt_liq
    write (6,901) "delt_glas  = ",delt_glas
    write (6,*) ""
    write (6,901) "density    = ",density
    write (6,901) "eigenvalue = ",eigenvalue
    write (6,901) "lambda     = ",lambda
    write (6,903) "inflexion  = ",inflex
    if (inflex .eqv. .true.) then
      write (6,901) "eigen_infl = ",eigenvalue_inflex
      write (6,901) "lamb_infl  = ",lambda_inflex
    end if
    write (6,901) "conv_fq    = ",conv_fq
    write (6,901) ""
    write (6,*) "ITERATION TIME ", iter_end-iter_start,"s"
    if (dabs(eigenvalue-1.0d0) .le. prec_eigen_disc) then
      exit
    end if
  end do
  write (6,*) ""
  write (6,800) "               CONVERGENCE REACHED                   "
  write (6,800) "-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-"
  write (6,901) "delt_liq   = ",delt_liq
  write (6,901) "delt_glas  = ",delt_glas
  write (6,*) ""
  write (6,901) "density    = ",density
  write (6,901) "eigenvalue = ",eigenvalue
  write (6,901) "lambda     = ",lambda
  write (6,903) "inflexion  = ",inflex
  if (inflex .eqv. .true.) then
    write (6,901) "eigen_infl = ",eigenvalue_inflex
    write (6,901) "lamb_infl  = ",lambda_inflex
  end if
  write (6,901) "conv_fq    = ",conv_fq
  call fileman('final_res',9,11,1)
  write (11,'(f22.16,f22.16,f22.16,f22.16)') density,density,delt_liq,delt_glas,lambda,eigenvalue
  call fileman('final_res',9,11,0)







!FAST DICH DELT DISC
else if ((calc_mode .eq. 'dich') .and. (fast .eqv. .true.) .and. (var_param .eq. 'delt') .and. (trans_mode .eq. 'disc')) then
  fast=.false.
  write (6,*) ""
  write (6,800) "                  FIRST RESULTS                      "
  write (6,800) "-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-"
  call struct (density,delta,sigma,closure,mix_param,cr_init)
  call eigen (trans_mode,density,eigenvalue,lambda,eigenvalue_inflex,lambda_inflex,inflex,fast,fq_init)
  write (6,*) ""
  write (6,901) "density    = ",density
  write (6,901) "delta      = ",delta
  write (6,901) "eigenvalue = ",eigenvalue
  write (6,901) "lambda     = ",lambda
  write (6,903) "inflexion  = ",inflex
  call fileman('fq.dat',6,11,1)
  do i=0,a_size
    read (11,*) q(i),fq(i)
  end do
  call fileman('fq.dat',6,11,0)
  conv_fq=fq(witb(fq,a_size))
  write (6,901) "conv_fq    = ",conv_fq
  if (conv_fq .lt. prec_fq) then
    do while (conv_fq .lt. prec_fq)
      delt_liq=delta
      delta=delta*var_incr
      write (6,*) ""
      write (6,800) "          Looking for the transition point           "
      write (6,800) "-----------------------------------------------------"
      call struct (density,delta,sigma,closure,mix_param,cr_init)
      call eigen (trans_mode,density,eigenvalue,lambda,eigenvalue_inflex,lambda_inflex,inflex,fast,fq_init)
      write (6,*) ""
      write (6,901) "delt_liq   = ",delt_liq
      write (6,901) "delt_glas  = ",delt_glas
      write (6,*) ""
      write (6,901) "density    = ",density
      write (6,901) "eigenvalue = ",eigenvalue
      write (6,901) "lambda     = ",lambda
      write (6,903) "inflexion  = ",inflex
      call fileman('fq.dat',6,11,1)
      do i=0,a_size
        read (11,*) q(i),fq(i)
      end do
      call fileman('fq.dat',6,11,0)
      conv_fq=fq(witb(fq,a_size))
      write (6,901) "conv_fq    = ",conv_fq
    end do
    delt_glas=delta
  else if (conv_fq .gt. prec_fq) then
    do while (conv_fq .gt. prec_fq)
      delt_glas=delta
      delta=delta/var_incr
      write (6,*) ""
      write (6,800) "          Looking for the transition point           "
      write (6,800) "-----------------------------------------------------"
      call struct (density,delta,sigma,closure,mix_param,cr_init)
      call eigen (trans_mode,density,eigenvalue,lambda,eigenvalue_inflex,lambda_inflex,inflex,fast,fq_init)
      write (6,*) ""
      write (6,901) "delt_liq   = ",delt_liq
      write (6,901) "delt_glas  = ",delt_glas
      write (6,901) ""
      write (6,901) "density    = ",density
      write (6,901) "eigenvalue = ",eigenvalue
      write (6,901) "lambda     = ",lambda
      write (6,903) "inflexion  = ",inflex
      call fileman('fq.dat',6,11,1)
      do i=0,a_size
        read (11,*) q(i),fq(i)
      end do
      call fileman('fq.dat',6,11,0)
      conv_fq=fq(witb(fq,a_size))
      write (6,901) "conv_fq    = ",conv_fq
    end do
    delt_liq=delta
  end if
  fast=.true.
  write (6,*) ""
  write (6,800) "-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-"
  write (6,800) "||           ENTERING THE DICHOTOMY                ||"
  write (6,800) "-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-"
  write (6,*) ""
  write (6,800) "STARTING PARAMETERS                                  "
  write (6,*) ""
  write (6,901) "delt_liq   = ",delt_liq
  write (6,901) "delt_glas  = ",delt_glas
  write (6,901) ""
  write (6,901) "density    = ",density
  write (6,901) "eigenvalue = ",eigenvalue
  write (6,901) "lambda     = ",lambda
  write (6,903) "inflexion  = ",inflex
  if (inflex .eqv. .false.) then
    call fileman('fq.dat',6,11,1)
    do i=0,a_size
      read (11,*) q(i),fq(i)
    end do
    call fileman('fq.dat',6,11,0)
    conv_fq=fq(witb(fq,a_size))
    write (6,901) "conv_fq    = ",conv_fq
  end if
  do 
    call cpu_time(iter_start)
    iter=iter+1
    delta=(delt_liq+delt_glas)/2.0d0
    write (6,*) ""
    write (6,*) "                 ITERATION ",iter
    write (6,800) "-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-"
    call struct (density,delta,sigma,closure,mix_param,cr_init)
    call eigen (trans_mode,density,eigenvalue,lambda,eigenvalue_inflex,lambda_inflex,inflex,fast,fq_init)
    if (inflex .eqv. .true.) then
      delt_liq=delta
    else if (inflex .eqv. .false.) then
      delt_glas=delta
    end if
    call cpu_time(iter_end)
    write (6,*) ""
    write (6,901) "delt_liq   = ",delt_liq
    write (6,901) "delt_glas  = ",delt_glas
    write (6,*) ""
    write (6,901) "density    = ",density
    write (6,901) "eigenvalue = ",eigenvalue
    write (6,901) "lambda     = ",lambda
    write (6,903) "inflexion  = ",inflex
    if (inflex .eqv. .false.) then
      call fileman('fq.dat',6,11,1)
      do i=0,a_size
        read (11,*) q(i),fq(i)
      end do
      call fileman('fq.dat',6,11,0)
      conv_fq=fq(witb(fq,a_size))
      write (6,901) "conv_fq    = ",conv_fq
    end if
    write (6,901) ""
    write (6,*) "ITERATION TIME ", iter_end-iter_start,"s"
    if (dabs(eigenvalue-1.0d0) .le. prec_eigen_disc) then
      exit
    end if
  end do
  write (6,*) ""
  write (6,800) "               CONVERGENCE REACHED                   "
  write (6,800) "-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-"
  write (6,901) "delt_liq   = ",delt_liq
  write (6,901) "delt_glas  = ",delt_glas
  write (6,*) ""
  write (6,901) "density    = ",density
  write (6,901) "eigenvalue = ",eigenvalue
  write (6,901) "lambda     = ",lambda
  write (6,903) "inflexion  = ",inflex
  if (inflex .eqv. .false.) then
    call fileman('fq.dat',6,11,1)
    do i=0,a_size
      read (11,*) q(i),fq(i)
    end do
    call fileman('fq.dat',6,11,0)
    conv_fq=fq(witb(fq,a_size))
    write (6,901) "conv_fq    = ",conv_fq
  end if
  call fileman('final_res',9,11,1)
  write (11,'(f22.16,f22.16,f22.16,f22.16)') density,density,delt_liq,delt_glas,lambda,eigenvalue
  call fileman('final_res',9,11,0)


















!DICH DENS CONT/LOCA
!=================================================
else if ((calc_mode .eq. 'dich') .and. (var_param .eq. 'dens') .and. ((trans_mode .eq. 'cont') .or. (trans_mode .eq. 'loca'))) then
  write (6,*) ""
  write (6,800) "                  FIRST RESULTS                      "
  write (6,800) "-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-"
  call struct (density,delta,sigma,closure,mix_param,cr_init)
  call eigen (trans_mode,density,eigenvalue,lambda,eigenvalue_inflex,lambda_inflex,inflex,fast,fq_init)
  write (6,901) "density    = ",density
  write (6,901) "delta      = ",delta
  write (6,901) "eigenvalue = ",eigenvalue
  write (6,901) "lambda     = ",lambda
  if (eigenvalue .gt. 1.0d0) then
    do while (eigenvalue .gt. 1.0d0)
      dens_liq=density
      density=density*var_incr
      write (6,*) ""
      write (6,800) "          Looking for the transition point           "
      write (6,800) "-----------------------------------------------------"
      call struct (density,delta,sigma,closure,mix_param,cr_init)
      call eigen (trans_mode,density,eigenvalue,lambda,eigenvalue_inflex,lambda_inflex,inflex,fast,fq_init)
      write (6,901) "dens_liq   = ",dens_liq
      write (6,901) "dens_glas  = ",dens_glas
      write(6,*) ""
      write (6,901) "delta      = ",delta
      write (6,901) "eigenvalue = ",eigenvalue
      write (6,901) "lambda     = ",lambda
    end do
    dens_glas=density
  else if (eigenvalue .lt. 1.0d0) then
    do while (eigenvalue .lt. 1.0d0)
      dens_glas=density
      density=density/var_incr
      write (6,*) ""
      write (6,800) "          Looking for the transition point           "
      write (6,800) "-----------------------------------------------------"
      call struct (density,delta,sigma,closure,mix_param,cr_init)
      call eigen (trans_mode,density,eigenvalue,lambda,eigenvalue_inflex,lambda_inflex,inflex,fast,fq_init)
      write (6,901) "dens_liq   = ",dens_liq
      write (6,901) "dens_glas  = ",dens_glas
      write (6,901) ""
      write (6,901) "delta      = ",delta
      write (6,901) "eigenvalue = ",eigenvalue
      write (6,901) "lambda     = ",lambda
    end do
    dens_liq=density
  end if
  write (6,*) ""
  write (6,800) "-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-"
  write (6,800) "||           ENTERING THE DICHOTOMY                ||"
  write (6,800) "-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-"
  write (6,*) ""
  write (6,800) "STARTING PARAMETERS                                  "
  write (6,*) ""
  write (6,901) "dens_liq   = ",dens_liq
  write (6,901) "dens_glas  = ",dens_glas
  write (6,901) ""
  write (6,901) "delta      = ",delta
  write (6,901) "eigenvalue = ",eigenvalue
  write (6,901) "lambda     = ",lambda
  do 
    call cpu_time(iter_start)
    iter=iter+1
    density=(dens_liq+dens_glas)/2.0d0
    write (6,*) ""
    write (6,*) "                 ITERATION ",iter
    write (6,800) "-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-"
    call struct (density,delta,sigma,closure,mix_param,cr_init)
    call eigen (trans_mode,density,eigenvalue,lambda,eigenvalue_inflex,lambda_inflex,inflex,fast,fq_init)
    if (eigenvalue .gt. 1.0d0) then
      dens_liq=density
    else if (eigenvalue .lt. 1.0d0) then
      dens_glas=density
    end if
    call cpu_time(iter_end)
    write (6,901) "dens_liq   = ",dens_liq
    write (6,901) "dens_glas  = ",dens_glas
    write (6,901) ""
    write (6,901) "delta      = ",delta
    write (6,901) "eigenvalue = ",eigenvalue
    write (6,901) "lambda     = ",lambda
    write (6,901) ""
    write (6,*) "ITERATION TIME ", iter_end-iter_start,"s"
    if (dabs(eigenvalue-1.0d0) .le. prec_eigen_cont) then
      exit
    end if
  end do
  write (6,*) ""
  write (6,800) "               CONVERGENCE REACHED                   "
  write (6,800) "-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-"
  write (6,901) "dens_liq =   ",dens_liq
  write (6,901) "dens_glas =  ",dens_glas
  write (6,901) ""
  write (6,901) "delta      = ",delta
  write (6,901) "eigenvalue = ",eigenvalue
  write (6,901) "lambda     = ",lambda
  call fileman('final_res',9,11,1)
  write (11,'(f22.16,f22.16,f22.16,f22.16)') dens_liq,dens_glas,delta,delta,lambda,eigenvalue
  call fileman('final_res',9,11,0)











!DICH DENS DISC
!=================================================
else if ((calc_mode .eq. 'dich') .and. (var_param .eq. 'dens') .and. (trans_mode .eq. 'disc') .and. (fast .eqv. .false.)) then
  eigenvalue=0.0d0
  lambda=0.0d0
  eigenvalue_inflex=0.0d0
  lambda_inflex=0.0d0
  write (6,*) ""
  write (6,800) "                  FIRST RESULTS                      "
  write (6,800) "-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-"
  call struct (density,delta,sigma,closure,mix_param,cr_init)
  call eigen (trans_mode,density,eigenvalue,lambda,eigenvalue_inflex,lambda_inflex,inflex,fast,fq_init)
  call fileman('fq.dat',6,11,1)
  do i=0,a_size
    read (11,*) q(i),fq(i)
  end do
  call fileman('fq.dat',6,11,0)
  conv_fq=fq(witb(fq,a_size))
  write (6,*) ""
  write (6,901) "density    = ",density
  write (6,901) "delta      = ",delta
  write (6,901) "eigenvalue = ",eigenvalue
  write (6,901) "lambda     = ",lambda
  write (6,903) "inflexion  = ",inflex
  if (inflex .eqv. .true.) then
    write (6,901) "eigen_infl = ",eigenvalue_inflex
    write (6,901) "lamb_infl  = ",lambda_inflex
    inflex=.false.
  end if
  write (6,902) "conv_fq    = ",conv_fq
  if (conv_fq .lt. prec_fq) then
    do while (conv_fq .lt. prec_fq)
      dens_liq=density
      density=density*var_incr
      write (6,*) ""
      write (6,800) "          Looking for the transition point           "
      write (6,800) "-----------------------------------------------------"
      call struct (density,delta,sigma,closure,mix_param,cr_init)
      call eigen (trans_mode,density,eigenvalue,lambda,eigenvalue_inflex,lambda_inflex,inflex,fast,fq_init)
      call fileman('fq.dat',6,11,1)
      do i=0,a_size
        read (11,*) q(i),fq(i)
      end do
      call fileman('fq.dat',6,11,0)
      conv_fq=fq(witb(fq,a_size))
      write (6,*) ""
      write (6,901) "dens_liq   = ",dens_liq
      write (6,901) "dens_glas  = ",dens_glas
      write(6,*) ""
      write (6,901) "delta      = ",delta
      write (6,901) "eigenvalue = ",eigenvalue
      write (6,901) "lambda     = ",lambda
      write (6,903) "inflexion  = ",inflex
      if (inflex .eqv. .true.) then
        write (6,901) "eigen_infl = ",eigenvalue_inflex
        write (6,901) "lamb_infl  = ",lambda_inflex
      end if
      write (6,902) "conv_fq    = ",conv_fq
    end do
    dens_glas=density
  else if (conv_fq .gt. prec_fq) then
    do while (conv_fq .gt. prec_fq)
      dens_glas=density
      density=density/var_incr
      write (6,*) ""
      write (6,800) "          Looking for the transition point           "
      write (6,800) "-----------------------------------------------------"
      call struct (density,delta,sigma,closure,mix_param,cr_init)
      call eigen (trans_mode,density,eigenvalue,lambda,eigenvalue_inflex,lambda_inflex,inflex,fast,fq_init)
      call fileman('fq.dat',6,11,1)
      do i=0,a_size
        read (11,*) q(i),fq(i)
      end do
      call fileman('fq.dat',6,11,0)
      conv_fq=fq(witb(fq,a_size))
      write (6,*) ""
      write (6,901) "dens_liq   = ",dens_liq
      write (6,901) "dens_glas  = ",dens_glas
      write (6,901) ""
      write (6,901) "delta      = ",delta
      write (6,901) "eigenvalue = ",eigenvalue
      write (6,901) "lambda     = ",lambda
      write (6,903) "inflexion  = ",inflex
      if (inflex .eqv. .true.) then
        write (6,901) "eigen_infl = ",eigenvalue_inflex
        write (6,901) "lamb_infl  = ",lambda_inflex
      end if
      write (6,902) "conv_fq    = ",conv_fq
    end do
    dens_liq=density
  end if
  write (6,*) ""
  write (6,800) "-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-"
  write (6,800) "||           ENTERING THE DICHOTOMY                ||"
  write (6,800) "-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-"
  write (6,*) ""
  write (6,800) "STARTING PARAMETERS                                  "
  write (6,*) ""
  write (6,901) "dens_liq   = ",dens_liq
  write (6,901) "dens_glas  = ",dens_glas
  write (6,901) ""
  write (6,901) "delta      = ",delta
  write (6,901) "eigenvalue = ",eigenvalue
  write (6,901) "lambda     = ",lambda
  write (6,903) "inflexion  = ",inflex
  if (inflex .eqv. .true.) then
    write (6,901) "eigen_infl = ",eigenvalue_inflex
    write (6,901) "lamb_infl  = ",lambda_inflex
  end if
  write (6,902) "conv_fq    = ",conv_fq
  do 
    call cpu_time(iter_start)
    iter=iter+1
    density=(dens_liq+dens_glas)/2.0d0
    write (6,*) ""
    write (6,*) "                 ITERATION ",iter
    write (6,800) "-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-"
    call struct (density,delta,sigma,closure,mix_param,cr_init)
    call eigen (trans_mode,density,eigenvalue,lambda,eigenvalue_inflex,lambda_inflex,inflex,fast,fq_init)
    call fileman('fq.dat',6,11,1)
    do i=0,a_size
      read (11,*) q(i),fq(i)
    end do
    call fileman('fq.dat',6,11,0)
    conv_fq=fq(witb(fq,a_size))
    if (conv_fq .lt. prec_fq) then
      dens_liq=density
    else if (conv_fq .gt. prec_fq) then
      dens_glas=density
    end if
    call cpu_time(iter_end)
    write (6,*) ""
    write (6,901) "dens_liq   = ",dens_liq
    write (6,901) "dens_glas  = ",dens_glas
    write (6,901) ""
    write (6,901) "delta      = ",delta
    write (6,901) "eigenvalue = ",eigenvalue
    write (6,901) "lambda     = ",lambda
    write (6,903) "inflexion  = ",inflex
    if (inflex .eqv. .true.) then
      write (6,901) "eigen_infl = ",eigenvalue_inflex
      write (6,901) "lamb_infl  = ",lambda_inflex
    end if
    write (6,902) "conv_fq    = ",conv_fq
    write (6,901) ""
    write (6,*) "ITERATION TIME ", iter_end-iter_start,"s"
    if (dabs(eigenvalue-1.0d0) .le. prec_eigen_disc) then
      exit
    end if
  end do
  write (6,*) ""
  write (6,800) "               CONVERGENCE REACHED                   "
  write (6,800) "-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-"
  write (6,901) "dens_liq   = ",dens_liq
  write (6,901) "dens_glas  = ",dens_glas
  write (6,901) ""
  write (6,901) "delta      = ",delta
  write (6,901) "eigenvalue = ",eigenvalue
  write (6,901) "lambda     = ",lambda
  write (6,903) "inflexion  = ",inflex
 if (inflex .eqv. .true.) then
    write (6,901) "eigen_infl = ",eigenvalue_inflex
    write (6,901) "lamb_infl  = ",lambda_inflex
  end if
  write (6,902) "conv_fq    = ",conv_fq
  call fileman('final_res',9,11,1)
  write (11,'(f22.16,f22.16,f22.16,f22.16)') dens_liq,dens_glas,delta,delta,lambda,eigenvalue
  call fileman('final_res',9,11,0)













!FAST DICH DENS DISC
!=================================================
else if ((calc_mode .eq. 'dich') .and. (fast .eqv. .true.) .and. (var_param .eq. 'dens') .and. (trans_mode .eq. 'disc')) then
  fast=.false.
  write (6,*) ""
  write (6,800) "                  FIRST RESULTS                      "
  write (6,800) "-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-"
  call struct (density,delta,sigma,closure,mix_param,cr_init)
  call eigen (trans_mode,density,eigenvalue,lambda,eigenvalue_inflex,lambda_inflex,inflex,fast,fq_init)
  write (6,*) ""
  write (6,901) "density    = ",density
  write (6,901) "delta      = ",delta
  write (6,901) "eigenvalue = ",eigenvalue
  write (6,901) "lambda     = ",lambda
  write (6,903) "inflexion  = ",inflex
  call fileman('fq.dat',6,11,1)
  do i=0,a_size
    read (11,*) q(i),fq(i)
  end do
  call fileman('fq.dat',6,11,0)
  conv_fq=fq(witb(fq,a_size))
  write (6,901) "conv_fq    = ",conv_fq
  if (conv_fq .lt. prec_fq) then
    do while (conv_fq .lt. prec_fq)
      dens_liq=density
      density=density*var_incr
      write (6,*) ""
      write (6,800) "          Looking for the transition point           "
      write (6,800) "-----------------------------------------------------"
      call struct (density,delta,sigma,closure,mix_param,cr_init)
      call eigen (trans_mode,density,eigenvalue,lambda,eigenvalue_inflex,lambda_inflex,inflex,fast,fq_init)
      write (6,*) ""
      write (6,901) "dens_liq   = ",dens_liq
      write (6,901) "dens_glas  = ",dens_glas
      write(6,*) ""
      write (6,901) "delta      = ",delta
      write (6,901) "eigenvalue = ",eigenvalue
      write (6,901) "lambda     = ",lambda
      write (6,903) "inflexion  = ",inflex
      call fileman('fq.dat',6,11,1)
      do i=0,a_size
        read (11,*) q(i),fq(i)
      end do
      call fileman('fq.dat',6,11,0)
      conv_fq=fq(witb(fq,a_size))
      write (6,901) "conv_fq    = ",conv_fq
    end do
    dens_glas=density
  else if (conv_fq .gt. prec_fq) then
    do while (conv_fq .gt. prec_fq)
      dens_glas=density
      density=density/var_incr
      write (6,*) ""
      write (6,800) "          Looking for the transition point           "
      write (6,800) "-----------------------------------------------------"
      call struct (density,delta,sigma,closure,mix_param,cr_init)
      call eigen (trans_mode,density,eigenvalue,lambda,eigenvalue_inflex,lambda_inflex,inflex,fast,fq_init)
      write (6,*) ""
      write (6,901) "dens_liq   = ",dens_liq
      write (6,901) "dens_glas  = ",dens_glas
      write (6,901) ""
      write (6,901) "delta      = ",delta
      write (6,901) "eigenvalue = ",eigenvalue
      write (6,901) "lambda     = ",lambda
      write (6,903) "inflexion  = ",inflex
      call fileman('fq.dat',6,11,1)
      do i=0,a_size
        read (11,*) q(i),fq(i)
      end do
      call fileman('fq.dat',6,11,0)
      conv_fq=fq(witb(fq,a_size))
      write (6,901) "conv_fq    = ",conv_fq
    end do
    dens_liq=density
  end if
  fast=.true.
  write (6,*) ""
  write (6,800) "-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-"
  write (6,800) "||           ENTERING THE DICHOTOMY                ||"
  write (6,800) "-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-"
  write (6,*) ""
  write (6,800) "STARTING PARAMETERS                                  "
  write (6,*) ""
  write (6,901) "dens_liq   = ",dens_liq
  write (6,901) "dens_glas  = ",dens_glas
  write (6,901) ""
  write (6,901) "delta      = ",delta
  write (6,901) "eigenvalue = ",eigenvalue
  write (6,901) "lambda     = ",lambda
  write (6,903) "inflexion  = ",inflex
  if (inflex .eqv. .false.) then
    call fileman('fq.dat',6,11,1)
    do i=0,a_size
      read (11,*) q(i),fq(i)
    end do
    call fileman('fq.dat',6,11,0)
    conv_fq=fq(witb(fq,a_size))
    write (6,901) "conv_fq    = ",conv_fq
  end if
  do 
    call cpu_time(iter_start)
    iter=iter+1
    density=(dens_liq+dens_glas)/2.0d0
    write (6,*) ""
    write (6,*) "                 ITERATION ",iter
    write (6,800) "-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-"
    call struct (density,delta,sigma,closure,mix_param,cr_init)
    call eigen (trans_mode,density,eigenvalue,lambda,eigenvalue_inflex,lambda_inflex,inflex,fast,fq_init)
    if (inflex .eqv. .true.) then
      dens_liq=density
    else if (inflex .eqv. .false.) then
      dens_glas=density
    end if
    call cpu_time(iter_end)
    write (6,*) ""
    write (6,901) "dens_liq   = ",dens_liq
    write (6,901) "dens_glas  = ",dens_glas
    write (6,901) ""
    write (6,901) "delta      = ",delta
    write (6,901) "eigenvalue = ",eigenvalue
    write (6,901) "lambda     = ",lambda
    write (6,903) "inflexion  = ",inflex
    if (inflex .eqv. .false.) then
      call fileman('fq.dat',6,11,1)
      do i=0,a_size
        read (11,*) q(i),fq(i)
      end do
      call fileman('fq.dat',6,11,0)
      conv_fq=fq(witb(fq,a_size))
      write (6,901) "conv_fq    = ",conv_fq
    end if
    write (6,901) ""
    write (6,*) "ITERATION TIME ", iter_end-iter_start,"s"
    if (dabs(eigenvalue-1.0d0) .le. prec_eigen_disc) then
      exit
    end if
  end do
  write (6,*) ""
  write (6,800) "               CONVERGENCE REACHED                   "
  write (6,800) "-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-"
  write (6,901) "dens_liq   = ",dens_liq
  write (6,901) "dens_glas  = ",dens_glas
  write (6,901) ""
  write (6,901) "delta      = ",delta
  write (6,901) "eigenvalue = ",eigenvalue
  write (6,901) "lambda     = ",lambda
  write (6,903) "inflexion  = ",inflex
  if (inflex .eqv. .false.) then
    call fileman('fq.dat',6,11,1)
    do i=0,a_size
      read (11,*) q(i),fq(i)
    end do
    call fileman('fq.dat',6,11,0)
    conv_fq=fq(witb(fq,a_size))
    write (6,901) "conv_fq    = ",conv_fq
  end if
  call fileman('final_res',9,11,1)
  write (11,'(f22.16,f22.16,f22.16,f22.16)') dens_liq,dens_glas,delta,delta,lambda,eigenvalue
  call fileman('final_res',9,11,0)




















!REST DELT DISC
!=================================================
else if ((trans_mode .eq. 'disc') .and. (calc_mode .eq. 'rest') .and. (var_param .eq. 'delt') .and. (fast .eqv. .false.)) then
 conv_fq=1.0d0
  write (6,*) ""
  write (6,800) "-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-"
  write (6,800) "||           ENTERING THE DICHOTOMY                ||"
  write (6,800) "-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-"
  write (6,*) ""
  write (6,800) "STARTING PARAMETERS                                  "
  write (6,*) ""
  write (6,901) "delt_liq   = ",delt_liq
  write (6,901) "delt_glas  = ",delt_glas
  write (6,901) ""
  write (6,901) "density    = ",density
  write (6,901) "eigenvalue = ",eigenvalue
  write (6,901) "lambda     = ",lambda
  write (6,903) "inflexion  = ",inflex
  if (inflex .eqv. .true.) then
    write (6,901) "eigen_infl = ",eigenvalue_inflex
    write (6,901) "lamb_infl  = ",lambda_inflex
  end if
  write (6,901) "conv_fq    = ",conv_fq
  do 
    call cpu_time(iter_start)
    iter=iter+1
    delta=(delt_liq+delt_glas)/2.0d0
    write (6,*) ""
    write (6,*) "                 ITERATION ",iter
    write (6,800) "-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-"
    call struct (density,delta,sigma,closure,mix_param,cr_init)
    call eigen (trans_mode,density,eigenvalue,lambda,eigenvalue_inflex,lambda_inflex,inflex,fast,fq_init)
    call fileman('fq.dat',6,11,1)
    do i=0,a_size
      read (11,*) q(i),fq(i)
    end do
    call fileman('fq.dat',6,11,0)
    conv_fq=fq(witb(fq,a_size))
    if (conv_fq .lt. prec_fq) then
      delt_liq=delta
    else if (conv_fq .gt. prec_fq) then
      delt_glas=delta
    end if
    call cpu_time(iter_end)
    write (6,*) ""
    write (6,901) "delt_liq   = ",delt_liq
    write (6,901) "delt_glas  = ",delt_glas
    write (6,*) ""
    write (6,901) "density    = ",density
    write (6,901) "eigenvalue = ",eigenvalue
    write (6,901) "lambda     = ",lambda
    write (6,903) "inflexion  = ",inflex
    if (inflex .eqv. .true.) then
      write (6,901) "eigen_infl = ",eigenvalue_inflex
      write (6,901) "lamb_infl  = ",lambda_inflex
    end if
    write (6,901) "conv_fq    = ",conv_fq
    write (6,901) ""
    write (6,*) "ITERATION TIME ", iter_end-iter_start,"s"
    if (dabs(eigenvalue-1.0d0) .le. prec_eigen_disc) then
      exit
    end if
  end do
  write (6,*) ""
  write (6,800) "               CONVERGENCE REACHED                   "
  write (6,800) "-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-"
  write (6,901) "delt_liq   = ",delt_liq
  write (6,901) "delt_glas  = ",delt_glas
  write (6,*) ""
  write (6,901) "density    = ",density
  write (6,901) "eigenvalue = ",eigenvalue
  write (6,901) "lambda     = ",lambda
  write (6,903) "inflexion  = ",inflex
  if (inflex .eqv. .true.) then
    write (6,901) "eigen_infl = ",eigenvalue_inflex
    write (6,901) "lamb_infl  = ",lambda_inflex
  end if
  write (6,901) "conv_fq    = ",conv_fq
  call fileman('final_res',9,11,1)
  write (11,'(f22.16,f22.16,f22.16,f22.16)') density,density,delt_liq,delt_glas,lambda,eigenvalue
  call fileman('final_res',9,11,0)











!FAST REST DELT DISC
!=================================================
else if ((trans_mode .eq. 'disc') .and. (calc_mode .eq. 'rest') .and. (var_param .eq. 'delt') .and. (fast .eqv. .true.)) then
 conv_fq=1.0d0
  write (6,*) ""
  write (6,800) "-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-"
  write (6,800) "||           ENTERING THE DICHOTOMY                ||"
  write (6,800) "-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-"
  write (6,*) ""
  write (6,800) "STARTING PARAMETERS                                  "
  write (6,*) ""
  write (6,901) "delt_liq   = ",delt_liq
  write (6,901) "delt_glas  = ",delt_glas
  write (6,901) ""
  write (6,901) "density    = ",density
  write (6,901) "eigenvalue = ",eigenvalue
  write (6,901) "lambda     = ",lambda
  write (6,903) "inflexion  = ",inflex
  if (inflex .eqv. .false.) then
    call fileman('fq.dat',6,11,1)
    do i=0,a_size
      read (11,*) q(i),fq(i)
    end do
    call fileman('fq.dat',6,11,0)
    conv_fq=fq(witb(fq,a_size))
    write (6,901) "conv_fq    = ",conv_fq
  end if
  do 
    call cpu_time(iter_start)
    iter=iter+1
    delta=(delt_liq+delt_glas)/2.0d0
    write (6,*) ""
    write (6,*) "                 ITERATION ",iter
    write (6,800) "-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-"
    call struct (density,delta,sigma,closure,mix_param,cr_init)
    call eigen (trans_mode,density,eigenvalue,lambda,eigenvalue_inflex,lambda_inflex,inflex,fast,fq_init)
    call fileman('fq.dat',6,11,1)
    do i=0,a_size
      read (11,*) q(i),fq(i)
    end do
    call fileman('fq.dat',6,11,0)
    conv_fq=fq(witb(fq,a_size))
    if (inflex .eqv. .true.) then
      delt_liq=delta
    else if (inflex .eqv. .false.) then
      delt_glas=delta
    end if
    call cpu_time(iter_end)
    write (6,*) ""
    write (6,901) "delt_liq   = ",delt_liq
    write (6,901) "delt_glas  = ",delt_glas
    write (6,*) ""
    write (6,901) "density    = ",density
    write (6,901) "eigenvalue = ",eigenvalue
    write (6,901) "lambda     = ",lambda
    write (6,903) "inflexion  = ",inflex
    if (inflex .eqv. .false.) then
      call fileman('fq.dat',6,11,1)
      do i=0,a_size
        read (11,*) q(i),fq(i)
      end do
      call fileman('fq.dat',6,11,0)
      conv_fq=fq(witb(fq,a_size))
      write (6,901) "conv_fq    = ",conv_fq
    end if
    write (6,901) ""
    write (6,*) "ITERATION TIME ", iter_end-iter_start,"s"
    if (dabs(eigenvalue-1.0d0) .le. prec_eigen_disc) exit
  end do
  write (6,*) ""
  write (6,800) "               CONVERGENCE REACHED                   "
  write (6,800) "-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-"
  write (6,901) "delt_liq   = ",delt_liq
  write (6,901) "delt_glas  = ",delt_glas
  write (6,*) ""
  write (6,901) "density    = ",density
  write (6,901) "eigenvalue = ",eigenvalue
  write (6,901) "lambda     = ",lambda
  write (6,903) "inflexion  = ",inflex
  if (inflex .eqv. .false.) then
    call fileman('fq.dat',6,11,1)
    do i=0,a_size
      read (11,*) q(i),fq(i)
    end do
    call fileman('fq.dat',6,11,0)
    conv_fq=fq(witb(fq,a_size))
    write (6,901) "conv_fq    = ",conv_fq
  end if
  call fileman('final_res',9,11,1)
  write (11,'(f22.16,f22.16,f22.16,f22.16)') density,density,delt_liq,delt_glas,lambda,eigenvalue
  call fileman('final_res',9,11,0)





















!DISC REST DENS
else if ((trans_mode .eq. 'disc') .and. (calc_mode .eq. 'rest') .and. (var_param .eq. 'dens') .and. (fast .eqv. .false.)) then
  conv_fq=1.0d0
  write (6,*) ""
  write (6,800) "-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-"
  write (6,800) "||           ENTERING THE DICHOTOMY                ||"
  write (6,800) "-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-"
  write (6,*) ""
  write (6,800) "STARTING PARAMETERS                                  "
  write (6,*) ""
  write (6,901) "dens_liq   = ",dens_liq
  write (6,901) "dens_glas  = ",dens_glas
  write (6,901) ""
  write (6,901) "delta      = ",delta
  write (6,901) "eigenvalue = ",eigenvalue
  write (6,901) "lambda     = ",lambda
  write (6,903) "inflexion  = ",inflex
  if (inflex .eqv. .true.) then
    write (6,901) "eigen_infl = ",eigenvalue_inflex
    write (6,901) "lamb_infl  = ",lambda_inflex
  end if
  write (6,902) "conv_fq    = ",conv_fq
  do 
    call cpu_time(iter_start)
    iter=iter+1
    density=(dens_liq+dens_glas)/2.0d0
    write (6,*) ""
    write (6,*) "                 ITERATION ",iter
    write (6,800) "-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-"
    call struct (density,delta,sigma,closure,mix_param,cr_init)
    call eigen (trans_mode,density,eigenvalue,lambda,eigenvalue_inflex,lambda_inflex,inflex,fast,fq_init)
    call fileman('fq.dat',6,11,1)
    do i=0,a_size
      read (11,*) q(i),fq(i)
    end do
    call fileman('fq.dat',6,11,0)
    conv_fq=fq(witb(fq,a_size))
    if (conv_fq .lt. prec_fq) then
      dens_liq=density
    else if (conv_fq .gt. prec_fq) then
      dens_glas=density
    end if
    call cpu_time(iter_end)
    write (6,*) ""
    write (6,901) "dens_liq   = ",dens_liq
    write (6,901) "dens_glas  = ",dens_glas
    write (6,901) ""
    write (6,901) "delta      = ",delta
    write (6,901) "eigenvalue = ",eigenvalue
    write (6,901) "lambda     = ",lambda
    write (6,903) "inflexion  = ",inflex
    if (inflex .eqv. .true.) then
      write (6,901) "eigen_infl = ",eigenvalue_inflex
      write (6,901) "lamb_infl  = ",lambda_inflex
    end if
    write (6,902) "conv_fq    = ",conv_fq
    write (6,901) ""
    write (6,*) "ITERATION TIME ", iter_end-iter_start,"s"
    if (dabs(eigenvalue-1.0d0) .le. prec_eigen_disc) then
      exit
    end if
  end do
  write (6,*) ""
  write (6,800) "               CONVERGENCE REACHED                   "
  write (6,800) "-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-"
  write (6,901) "dens_liq   = ",dens_liq
  write (6,901) "dens_glas  = ",dens_glas
  write (6,901) ""
  write (6,901) "delta      = ",delta
  write (6,901) "eigenvalue = ",eigenvalue
  write (6,901) "lambda     = ",lambda
  write (6,903) "inflexion  = ",inflex
  if (inflex .eqv. .true.) then
    write (6,901) "eigen_infl = ",eigenvalue_inflex
    write (6,901) "lamb_infl  = ",lambda_inflex
  end if
  write (6,902) "conv_fq    = ",conv_fq
  call fileman('final_res',9,11,1)
  write (11,'(f22.16,f22.16,f22.16,f22.16)') dens_liq,dens_glas,delta,delta,lambda,eigenvalue
  call fileman('final_res',9,11,0)













!FAST DISC REST DENS
else if ((trans_mode .eq. 'disc') .and. (calc_mode .eq. 'rest') .and. (var_param .eq. 'dens') .and. (fast .eqv. .true.)) then
  conv_fq=1.0d0
  write (6,*) ""
  write (6,800) "-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-"
  write (6,800) "||           ENTERING THE DICHOTOMY                ||"
  write (6,800) "-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-"
  write (6,*) ""
  write (6,800) "STARTING PARAMETERS                                  "
  write (6,*) ""
  write (6,901) "dens_liq   = ",dens_liq
  write (6,901) "dens_glas  = ",dens_glas
  write (6,901) ""
  write (6,901) "delta      = ",delta
  write (6,901) "eigenvalue = ",eigenvalue
  write (6,901) "lambda     = ",lambda
  write (6,903) "inflexion  = ",inflex
  if (inflex .eqv. .false.) then
    call fileman('fq.dat',6,11,1)
    do i=0,a_size
      read (11,*) q(i),fq(i)
    end do
    call fileman('fq.dat',6,11,0)
    conv_fq=fq(witb(fq,a_size))
    write (6,901) "conv_fq    = ",conv_fq
  end if
  do 
    call cpu_time(iter_start)
    iter=iter+1
    density=(dens_liq+dens_glas)/2.0d0
    write (6,*) ""
    write (6,*) "                 ITERATION ",iter
    write (6,800) "-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-"
    call struct (density,delta,sigma,closure,mix_param,cr_init)
    call eigen (trans_mode,density,eigenvalue,lambda,eigenvalue_inflex,lambda_inflex,inflex,fast,fq_init)
    call fileman('fq.dat',6,11,1)
    do i=0,a_size
      read (11,*) q(i),fq(i)
    end do
    call fileman('fq.dat',6,11,0)
    conv_fq=fq(witb(fq,a_size))
    if (inflex .eqv. .true.) then
      dens_liq=density
    else if (fast .eqv. .false.) then
      dens_glas=density
    end if
    call cpu_time(iter_end)
    write (6,*) ""
    write (6,901) "dens_liq   = ",dens_liq
    write (6,901) "dens_glas  = ",dens_glas
    write (6,901) ""
    write (6,901) "delta      = ",delta
    write (6,901) "eigenvalue = ",eigenvalue
    write (6,901) "lambda     = ",lambda
    write (6,903) "inflexion  = ",inflex
    if (inflex .eqv. .false.) then
      call fileman('fq.dat',6,11,1)
      do i=0,a_size
        read (11,*) q(i),fq(i)
      end do
      call fileman('fq.dat',6,11,0)
      conv_fq=fq(witb(fq,a_size))
      write (6,901) "conv_fq    = ",conv_fq
    end if
    write (6,901) ""
    write (6,*) "ITERATION TIME ", iter_end-iter_start,"s"
    if (dabs(eigenvalue-1.0d0) .le. prec_eigen_disc) then
      exit
    end if
  end do
  write (6,*) ""
  write (6,800) "               CONVERGENCE REACHED                   "
  write (6,800) "-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-"
  write (6,901) "dens_liq   = ",dens_liq
  write (6,901) "dens_glas  = ",dens_glas
  write (6,901) ""
  write (6,901) "delta      = ",delta
  write (6,901) "eigenvalue = ",eigenvalue
  write (6,901) "lambda     = ",lambda
  write (6,903) "inflexion  = ",inflex
  if (inflex .eqv. .false.) then
    call fileman('fq.dat',6,11,1)
    do i=0,a_size
      read (11,*) q(i),fq(i)
    end do
    call fileman('fq.dat',6,11,0)
    conv_fq=fq(witb(fq,a_size))
    write (6,901) "conv_fq    = ",conv_fq
  end if
  call fileman('final_res',9,11,1)
  write (11,'(f22.16,f22.16,f22.16,f22.16)') dens_liq,dens_glas,delta,delta,lambda,eigenvalue
  call fileman('final_res',9,11,0)



























!DICH LAMB CONT
!=================================================
else if ((trans_mode .eq. 'cont') .and. (calc_mode .eq. 'dich') .and. (var_param .eq. 'lamb')) then
  write (6,*) ""
  write (6,800) "-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-"
  write (6,800) "||           ENTERING THE DICHOTOMY                ||"
  write (6,800) "-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-"
  write (6,*) ""
  write (6,800) "STARTING PARAMETERS                                  "
  write (6,*) ""
  write (6,901) "dens_lo    = ",dens_lo
  write (6,901) "delt_li_lo = ",delt_li_lo
  write (6,901) "delt_gl_lo = ",delt_gl_lo
  write (6,901) "lamb_lo    = ",lamb_lo
  write (6,*) ""                          
  write (6,901) "dens_hi    = ",dens_hi
  write (6,901) "delt_li_hi = ",delt_li_hi
  write (6,901) "delt_gl_hi = ",delt_gl_hi
  write (6,901) "lamb_hi    = ",lamb_hi
  do
    call cpu_time(iter_start)
    iter=iter+1
    density=(dens_lo+dens_hi)/2.0d0
    delt_lo=(delt_gl_lo+delt_li_lo)/2.0d0
    delt_hi=(delt_gl_hi+delt_li_hi)/2.0d0
    delta=(delt_lo+delt_hi)/2.0d0
!Run a dichotomy on the delta parameter
    call struct (density,delta,sigma,closure,mix_param,cr_init)
    call eigen (trans_mode,density,eigenvalue,lambda,eigenvalue_inflex,lambda_inflex,inflex,fast,fq_init)
    if (eigenvalue .lt. 1.0d0) then
      do while (eigenvalue .lt. 1.0d0)
        delt_liq=delta
        delta=delta*var_incr
    call struct (density,delta,sigma,closure,mix_param,cr_init)
    call eigen (trans_mode,density,eigenvalue,lambda,eigenvalue_inflex,lambda_inflex,inflex,fast,fq_init)
      end do
      delt_glas=delta
    else if (eigenvalue .gt. 1.0d0) then
      do while (eigenvalue .gt. 1.0d0)
        delt_glas=delta
        delta=delta/var_incr
        call struct (density,delta,sigma,closure,mix_param,cr_init)
        call eigen (trans_mode,density,eigenvalue,lambda,eigenvalue_inflex,lambda_inflex,inflex,fast,fq_init)
      end do
      delt_liq=delta
    end if
    do 
      delta=(delt_liq+delt_glas)/2.0d0
      call struct (density,delta,sigma,closure,mix_param,cr_init)
      call eigen (trans_mode,density,eigenvalue,lambda,eigenvalue_inflex,lambda_inflex,inflex,fast,fq_init)
      if (eigenvalue .lt. 1.0d0) then
        delt_liq=delta
      else if (eigenvalue .gt. 1.0d0) then
        delt_glas=delta
      end if
      if (dabs(eigenvalue-1.0d0) .le. prec_eigen_cont) exit
    end do
    if (lambda .lt. 1.0d0) then
      delt_li_lo=delt_liq
      delt_gl_lo=delt_glas
      dens_lo=density
      lamb_lo=lambda
    else if (lambda .gt. 1.0d0) then
      delt_li_hi=delt_liq
      delt_gl_hi=delt_glas
      dens_hi=density
      lamb_hi=lambda
    end if
    if (dabs(lambda-1.0d0) .lt. prec_lamb) exit
    call cpu_time(iter_end)
    write (6,*) ""
    write (6,*) "                 ITERATION ",iter
    write (6,800) "-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-"
    write (6,901) "dens_lo    = ",dens_lo
    write (6,901) "delt_li_lo = ",delt_li_lo
    write (6,901) "delt_gl_lo = ",delt_gl_lo
    write (6,901) "lamb_lo    = ",lamb_lo
    write (6,*) ""                          
    write (6,901) "dens_hi    = ",dens_hi
    write (6,901) "delt_li_hi = ",delt_li_hi
    write (6,901) "delt_gl_hi = ",delt_gl_hi
    write (6,901) "lamb_hi    = ",lamb_hi
    write (6,901) ""
    write (6,*) "ITERATION TIME ", iter_end-iter_start,"s"
  end do
  write (6,*) ""
  write (6,800) "               CONVERGENCE REACHED                   "
  write (6,800) "-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-"
  write (6,901) "dens_lo    = ",dens_lo
  write (6,901) "delt_li_lo = ",delt_li_lo
  write (6,901) "delt_gl_lo = ",delt_gl_lo
  write (6,901) "lamb_lo    = ",lamb_lo
  write (6,*) ""                          
  write (6,901) "dens_hi    = ",dens_hi
  write (6,901) "delt_li_hi = ",delt_li_hi
  write (6,901) "delt_gl_hi = ",delt_gl_hi
  write (6,901) "lamb_hi    = ",lamb_hi
end if

write (6,*) ""

call cpu_time(calc_end)
write (6,'(a17,f11.2,a)') "CALCULATION TIME ",calc_end,"s"

call date_and_time(DATE=d,TIME=t)
write (6,'(a21,a8,a,a10)') "CALCULATION ENDED ON ",d," ",t


end program
