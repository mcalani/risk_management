subroutine    externalit(h,nb,nss,is,B,tempbf,shadow,cph,extern)
use common_par

use nrtype
implicit none
integer ,intent(in) :: is    ,nb       ,h          ,nss
real(dp), intent(in) :: B(NB),tempbf ,shadow(nb,nss),cph(nb,nss)
real(dp),intent(out) :: extern
real(dp) ::  mu1(nss)  ,tempct2  , extern0  ,tempct   ,psi0,psi1(nss)    ,extern1  ,  externality1(nss)
integer :: iss


  tempct=  B(h)*tempr  + se(is)  -tempbf
  psi0= (1-omega)/omega*(mu+1)  *(tempct/(A*senn(is)))**mu  *kappan
  extern0=shadow(h,is)*psi0

  extern1=0.

  do iss=1,nss

  call linint(B,cph(:,iss),nb,tempbf,tempct2,1)
  psi1(iss)= (1-omega)/omega*(mu+1)*(tempct2/(A*senn(iss)))**mu *kappan

  call linint(B,shadow(:,iss),nb,tempbf,mu1(iss),1)

  externality1(iss)= betay*tempr* mu1(iss)* psi1(iss)
  extern1=  P(iss,is) *  externality1(iss)+extern1

  enddo

  extern=        extern1       - extern0

   end subroutine

   

