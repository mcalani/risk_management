  subroutine process(ns)        !  specification of markov process of exogenous shocks

  use common_par
  use global
  implicit none

  integer,intent(in) :: ns

  real(dp) :: ymat2(ns,nvar)
  integer is,iss

 call tauchenhussey(p,ymat2)                   ! calls tauchen hussey routine

 p=transpose(p)

 se=1+ymat2(:,1)
 senn=A*(1+ymat2(:,nvar))

! prints shock process

  open(unit=513,file='output\equil\shockprocess.txt')

  do is=1,ns
  do iss=1,ns
  write(513,*) p(iss,is)
  enddo
  enddo
  do is=1,ns
  write(513,*) se(is)
  enddo
  do is=1,ns
  write(513,*) senn(is)
  enddo
  close(unit=513)

  end subroutine

