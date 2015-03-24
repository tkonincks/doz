double precision function conv (x,f1,f2)
!Calculate the convergence of the system
!=================================================
double precision::x,f1,f2
double precision::pi=4.0d0*datan(1.0d0)
conv=4.0d0*pi*(x**2)*abs(f2-f1)
end function


double precision function trapez (a,b,fa,fb)
!Integrate a function between two points using the
!trapez method
!=================================================
double precision::a,b,fa,fb
trapez=(b-a)*(fa+fb)/2
end function


double precision function gaussian (sigma,x)
!Gaussian function for the potential of the particles
!=================================================
double precision::sigma,x
gaussian=dexp((-x**2)/(sigma**2))
end function

integer function witb (tab,d)
!Who is the biggest in this array?
!=================================================
  integer::d,i,bigi=0
  double precision,dimension(d)::tab
  double precision::big=-1.0d100
  do i=0,d
    if (tab(i) .gt. big) then
      big=tab(i)
      bigi=i
    end if
  end do
  witb=bigi
end function

integer function wits (tab,d)
!Who is the biggest in this array?
!=================================================
  integer::d,i,smalli
  double precision,dimension(d)::tab
  double precision::small=1.0d100
!  write (6,*) "d=",d
!  write (6,*) "small=",small
  do i=0,d
    if (tab(i) .lt. small) then
      small=tab(i)
      smalli=i
    end if
write (6,*) "small=",small
  end do
  wits=smalli
end function


logical function liq_glas (trans_mode,fast,eigenvalue,conv_fq,flag_inflex)
!returns .true. for the liquid state and .false. for the glas state
!=================================================
  character(len=4)::trans_mode
  double precision::eigenvalue,conv_fq
  double precision::prec_fq=1.0d-6
  logical::fast,flag_inflex

  if ((trans_mode .eq. 'cont') .or. (trans_mode .eq. 'loca')) then
    if (eigenvalue .lt. 1.0d0) then
      liq_glas=.true.
    else
      liq_glas=.false.
    end if

  else if ((trans_mode .eq. 'disc') .and. (fast .eqv. .true.)) then
    if (flag_inflex .eqv. .true.) then
      liq_glas=.true.
    else
      liq_glas=.false.
   end if

  else if (trans_mode .eq. 'disc') then
    if (conv_fq .lt. prec_fq) then
      liq_glas=.true.
    else
      liq_glas=.false.
    end if

  end if

end function
