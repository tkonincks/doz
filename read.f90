subroutine read (trans_mode,closure,calc_mode,var_param,density,delta,sigma,&
var_incr,var_prec,var_liq,var_glas,dens_lo,delt_li_lo,delt_gl_lo,dens_hi,&
delt_li_hi,delt_gl_hi,fast,mix_param,cr_init,fq_init,fcutoff,dynamics,tlimit,&
phi_init,correl)

implicit none

character(len=10)::varname
character(len=4)::trans_mode,closure,calc_mode,var_param,cr_init,fq_init,&
phi_init,correl
double precision::density,delta,sigma,var_incr,var_prec,var_liq,var_glas,&
dens_lo,delt_li_lo,delt_gl_lo,dens_hi,delt_li_hi,delt_gl_hi,mix_param,&
fcutoff,tlimit

logical::fast,dynamics

integer::io=0

!Initializations
!==================================================
trans_mode='xxxx'
calc_mode='xxxx'
var_param='xxxx'
closure='xxxx'
density=0.0d0
delta=0.0d0
sigma=0.0
var_liq=0.0
var_glas=0.0
dens_lo=0.0
delt_li_lo=0.0
delt_gl_lo=0.0
dens_hi=0.0
delt_li_hi=0.0
delt_gl_hi=0.0

!Some default values
!==================================================
fast=.false.
dynamics=.false.
var_incr=1.05
mix_param=0.95
cr_init='pyev'
fq_init='unit'
fcutoff=1.0d-9
tlimit=1.0d12
phi_init='none'
correl='gaus'

open (11,file='input.doz')

do

  read (11,'(a10)',advance='no',iostat=io) varname

  select case (varname)
    case ('trans_mode')
      read (11,*,iostat=io) trans_mode
    case ('closure')
      read (11,*,iostat=io) closure
    case ('calc_mode')
      read (11,*,iostat=io) calc_mode
    case ('var_param')
      read (11,*,iostat=io) var_param
    case ('var_incr')
      read (11,*,iostat=io) var_incr
    case ('var_prec')
      read (11,*,iostat=io) var_prec
    case ('density')
      read (11,*,iostat=io) density
    case ('delta')
      read (11,*,iostat=io) delta
    case ('sigma')
      read (11,*,iostat=io) sigma
    case ('var_liq')
      read (11,*,iostat=io) var_liq
    case ('var_glas')
      read (11,*,iostat=io) var_glas
    case ('dens_lo')
      read (11,*,iostat=io) dens_lo
    case ('delt_li_lo')
      read (11,*,iostat=io) delt_li_lo
    case ('delt_gl_lo')
      read (11,*,iostat=io) delt_gl_lo
    case ('dens_hi')
      read (11,*,iostat=io) dens_hi
    case ('delt_li_hi')
      read (11,*,iostat=io) delt_li_hi
    case ('delt_gl_hi')
      read (11,*,iostat=io) delt_gl_hi
    case ('fast')
      read (11,*,iostat=io) fast
    case ('mix_param')
      read (11,*,iostat=io) mix_param
    case ('cr_init')
      read (11,*,iostat=io) cr_init
    case ('fq_init')
      read (11,*,iostat=io) fq_init
    case ('fcutoff')
      read (11,*,iostat=io) fcutoff
    case ('dynamics')
      read (11,*,iostat=io) dynamics
    case ('tlimit')
      read (11,*,iostat=io) tlimit
    case ('phi_init')
      read (11,*,iostat=io) phi_init
    case ('correl')
      read (11,*,iostat=io) correl
  end select

  if (io .ne. 0) exit

end do

close (11)

end subroutine
