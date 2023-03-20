  !------------------------------------------------------------------------------------------------------
  !-----------------searches for a change in sign in nonlinear equation (based on numerical recipes)---------------
  !------------------------------------------------------------------------------------------------------
    
  subroutine zbrac(func,x1,x2,succes,nb,nss,is,h,b,cph,planner,shadow)
  use nrtype;
  use common_par
  implicit none
  real(sp), intent(inout) :: x1,x2
  integer ,intent(in) :: is ,h ,nb ,planner    ,nss
  real(dp), intent(in) :: b(nb),cph(nb,nss) ,shadow(nb,nss)
  logical(lgt), intent(out) :: succes

  integer(i4b), parameter :: ntry2=50
  real(sp), parameter :: factor=0.1_sp
  integer(i4b) :: j
  real(sp) :: f1,f2

  interface
  function  func(nb,nss,is,h,b,cph,cguess,planner,shadow)
  use nrtype
  use common_par
  implicit none
  integer ,intent(in) :: is   ,h ,nb   ,planner    ,nss
  real(dp), intent(in) :: b(nb),cph(nb,nss),cguess  ,shadow(nb,nss)
  real(sp) :: func
  end function func
  end interface



  f1=func(nb,nss,is,h,b,cph,x1,planner,shadow)
  f2=func(nb,nss,is,h,b,cph,x2,planner,shadow)

      succes=.true.

   do j=1,ntry
		if ((f1 > 0.0 .and. f2 < 0.0) .or. 	(f1 < 0.0 .and. f2 > 0.0))  return

	if (abs(f1) < abs(f2)) then
			x1=x1-0.02
			f1=func(nb,nss,is,h,b,cph,x1,planner,shadow)
    else
			x2=x2+0.02
			f2=func(nb,nss,is,h,b,cph,x2,planner,shadow)
		end if
  end do


  succes=.false.

  return

  end subroutine zbrac

