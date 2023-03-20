  ! Javier Bianchi 2010, "Overborrowing and Systemic Externalities in the Business Cycle"
  program main
  use nrtype;  use common_par
  implicit none
  real(dp), dimension(nb1,nss1)  :: bph ,plast,cph ,m,shadoww ,v,vsp   , taxx
  real(dp), dimension(nb1)  :: b
  integer ::    h,is      ,planner ,binding
  real (4) :: elapt, ta(2)
  real(dp) :: interp_value

  open(unit=111,file='output\parameters.txt')     ! stores parameters2
  write(111,*) sig,betay,mu,omega,rbar,kappat,kappan
  close(unit=111)

  call process(nss1)      ! stochastic process

  binding=1
  planner=0
  call grid_generate(nb1,b)

  call   bindingc(planner,nb1,nss1,b,v)              ! solve for DE
  call  plannervf(nb1,nss1,b,v)                     ! solve planner

  PRINT*,'--------------------------------------------------'
  	elapt = etime(ta)
  	write (*,*), ' '
  	write (*,'(a17, x, f8.2, x, a20)'), ' Program has used', elapt/60, 'minutes of CPU time,'
  	write (*,'(a5, x, f8.2, x, a24, x, f6.3, x, a23)'), ' with', ta(1)/60, 'minutes  of user time and', ta(2),&
    &  'seconds of system time.'
  PRINT*,'--------------------------------------------------'

  end program
