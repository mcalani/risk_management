subroutine    tax(h,nb,nss,is,B,tempbf,shadow,cph,taxtotal)
use common_par

use nrtype
implicit none
integer ,intent(in) :: h,is    ,nb    ,nss
real(dp), intent(in) :: B(NB),tempbf ,shadow(nb,nss),cph(nb,nss)
real(dp),intent(out) :: taxtotal

real(dp) ::  tmutil2,tempcn2,tempc2 ,tempct ,tempcc2    ,Etmutil2    ,extern   , tempct2
integer :: iss


   Etmutil2=0.

  tempct=  B(h)*tempr    + se(is) -tempbf

  do iss=1,nss

  call linint(B,cph(:,iss),nb,tempbf,tempct2,1)
  
  tempcn2=A*senn(iss)
  tempc2=omega*tempct2**(-mu)+(1-omega)*tempcn2**(-mu)
  tempcc2=tempc2**(-1/mu)
  tmutil2 = tempcc2**(-sig)* tempc2**(-1/mu-1)  *omega*tempct2**(-mu-1)
  Etmutil2= P(iss,is)*tmutil2 + Etmutil2

   enddo

  call  externalit(h,nb,nss,is,B,tempbf,shadow,cph,extern)

  taxtotal=  extern/(betay*Etmutil2*tempr)


   end subroutine


