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
