       subroutine linint(xtab,ytab,nb,xin,yout,extrap)
        use nrtype
      implicit none

      integer,intent(in)   :: nb  ,extrap
      real(dp),intent(in) :: xtab(nb),ytab(nb),xin

      real(dp), intent(out) :: yout
      real(dp) :: den,alfa,bbeta
      integer :: x1,x2


      call locate(xtab,nb,xin,x1)

      if (x1>0 .and. x1<nb) then
      x2=x1+1
      den=xtab(x2)-xtab(x1)
      alfa=(xtab(x2)-xin)/den
      bbeta= (xin-xtab(x1))/den
      yout=alfa*ytab(x1)+bbeta*ytab(x2)

      !  EXTRAPOLATION

      elseif (extrap==0 .and. x1>=NB) then
      yout=ytab(nb)

      elseif (x1==0 .and. extrap==0 ) then
      yout=ytab(1)

      elseif (extrap==1 .and. x1>=NB) then
      yout=ytab(nb-1)+(xin-xtab(nb-1))*(ytab(nb)-ytab(nb-1))/(xtab(nb)-xtab(nb-1))

      elseif (extrap==1  .and. x1==0 ) then
      yout=ytab(1)+(xin-xtab(1))*(ytab(2)-ytab(1))/(xtab(2)-xtab(1))

      endif
      
      
             end subroutine
         ! numerical recipes routine
  subroutine locate(xx,nb,x,j)
  use nrtype
  implicit none


  ! given an array xx of length nb, and a value of x, this routine returns
  ! a value j such that x is between xx(j) and xx(j+1). the array xx must be
  ! monotonic. j=0 or j=nb indicates that x is out of range. bisection is used
  ! to find the entry

  integer   ::        nb,j
  real(dp) :: xx(nb),x


  integer           jl,ju,jm

  jl = 0
  ju = nb+1

  10    if (ju-jl .gt. 1) then
  jm = (ju+jl)/2
  if ( (xx(nb) .ge. xx(1)) .eqv. (x .ge. xx(jm)) ) then
  jl = jm
  else
  ju = jm
  end if
  goto 10
  end if
  if (x .eq. xx(1))then
  j = 1
  else if(x .eq. xx(nb))then
  j = nb - 1
  else
  j = jl
  end if

  if (j>1) then
  if (xx(j-1)>=x) then
  j=j-1
  endif
  endif

  return
  end  SUBROUTINE
