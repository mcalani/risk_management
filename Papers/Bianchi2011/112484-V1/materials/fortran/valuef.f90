
  subroutine valuef(bph,b,v,nb,nss)
  use common_par

  use nrtype
  implicit none
  integer, intent(in):: nb    ,nss
  real(dp),intent(in) :: bph(nb,nss),b(nb)
  real(dp), intent(out) :: v(nb,nss)

  real(dp) :: vnew(nb,nss), y0(nss),diff_v ,tempb,util,expvv ,tempc,tempy,tempbf ,tempct ,tempcn
  integer :: iter,h,is,iss

  v=0.
  diff_v=100

  do iter=1,3000

  if (diff_v .lt. 1e-13) exit

  do is=1,nss
  do h=1,nb

  tempy=se(is)
  tempb=b(h)
  tempbf=bph(h,is)
  tempct= tempb*tempr - tempbf + tempy
  tempcn=A*senn(is)

  tempc=(omega*(tempct**(-mu))+(1-omega)*(tempcn**(-mu)))**(-1/mu)

  util=1./(1.-sig)*tempc**(1.-sig)

  do iss=1,nss
  call linint(b,v(:,iss),nb,tempbf,y0(iss),1)

  enddo

  expvv=0
  do  iss=1,nss
  expvv=p(iss,is)*y0(iss)+expvv
  enddo
  vnew(h,is)=util+betay*expvv
  

  enddo
  enddo

  diff_v= maxval(abs(v-vnew)/(1+abs(v)))
  v=vnew

  enddo

  if (diff_v .gt. 1.e-7) print*,'vf didnt converge'

  end subroutine

