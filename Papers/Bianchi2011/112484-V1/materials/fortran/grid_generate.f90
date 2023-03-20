  subroutine grid_generate(nb,b)                ! generates a grid of size nb

  use common_par
  use nrtype
  implicit none

  integer, intent(in) :: nb
  real(dp), intent(out) :: b(nb)
  integer  h

  b(1)=minb;
  do h=2,nb;
  b(h)=b(h-1)+(largestb-b(h-1))/(nb-h+1)**gsp
  enddo


  open(unit=51,file='output\gridd.txt')
  do h=1,nb
  write(51,*) b(h)
  enddo
  close(unit=51)

  end subroutine

