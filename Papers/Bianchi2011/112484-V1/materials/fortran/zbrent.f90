  ! ----------------------------------------------------------------------
  ! ------------solves non-linear equation (numerical recipes)------------
  ! ----------------------------------------------------------------------
  
  function zbrent(func,x1,x2,tol,nb,nss,is,h,grid,cph,planner,shadow)
  use nrtype

  implicit none
  real(sp), intent(in) :: x1,x2,tol
  integer ,intent(in) :: is   ,h ,nb                  ,planner    ,nss
  real(dp), intent(in) :: grid(nb),cph(nb,nss) ,shadow(nb,nss)

  real(sp) :: zbrent
  integer(i4b), parameter :: itmax=200
  real(dp), parameter :: eps=epsilon(x1)
  integer(i4b) :: iter
  real(dp) :: a,b,c,dd,e,fa,fb,fc,p,q,r,s,tol1,xm

  interface
  function func(nb,nss,is,h,b,cph,cguess,planner,shadow)
  use common_par
  use nrtype
  implicit none
  integer ,intent(in) :: is ,h ,nb ,planner ,nss
  real(dp), intent(in) :: b(nb),cph(nb,nss),cguess,shadow(nb,nss)

  real(dp) :: func
  end function func
  end interface

  external func

  a=x1
  b=x2

  fa=func(nb,nss,is,h,grid,cph,a,planner,shadow )

  fb=func(nb,nss,is,h,grid,cph,b,planner,shadow )

  if ((fa > 0.0 .and. fb > 0.0) .or. (fa < 0.0 .and. fb < 0.0)) then
  print*, 'root must be bracketed for zbrent'
  stop
  return
  endif

  c=b
  fc=fb
  do iter=1,itmax
  if ((fb > 0.0 .and. fc > 0.0) .or. (fb < 0.0 .and. fc < 0.0)) then
  c=a
  fc=fa
  dd=b-a
  e=dd
  end if
  if (abs(fc) < abs(fb)) then
  a=b
  b=c
  c=a
  fa=fb
  fb=fc
  fc=fa
  end if
  tol1=2.0_dp*eps*abs(b)+0.5_dp*tol
  xm=0.5_dp*(c-b)
  if (abs(xm) <= tol1 .or. fb == 0.0) then
  zbrent=b

  return
  end if
  if (abs(e) >= tol1 .and. abs(fa) > abs(fb)) then
  s=fb/fa
  if (a == c) then
  p=2.0_dp*xm*s
  q=1.0_dp-s
  else
  q=fa/fc
  r=fb/fc
  p=s*(2.0_dp*xm*q*(q-r)-(b-a)*(r-1.0_dp))
  q=(q-1.0_dp)*(r-1.0_dp)*(s-1.0_dp)
  end if
  if (p > 0.0) q=-q
  p=abs(p)
  if (2.0_dp*p  <  min(3.0_dp*xm*q-abs(tol1*q),abs(e*q))) then
  e=dd
  dd=p/q
  else
  dd=xm
  e=dd
  end if
  else
  dd=xm
  e=dd
  end if
  a=b
  fa=fb
  b=b+merge(dd,sign(tol1,xm), abs(dd) > tol1 )
  fb=func(nb,nss,is,h,grid,cph,b,planner,shadow )
  end do
  print*,    'zbrent: exceeded maximum iterations'

  zbrent=b

  end function zbrent
