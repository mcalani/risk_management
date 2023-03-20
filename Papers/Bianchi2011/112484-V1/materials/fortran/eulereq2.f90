 function eulereq2(nb,nss,is,h,b,cph,cguess,planner,shadow)
use common_par; use nrtype

implicit none
integer ,intent(in) :: is   ,h ,nb ,planner    ,nss
real(dp), intent(in) :: b(nb),cph(nb,nss),cguess ,shadow(nb,nss)
real(dp) ::  tempbf,lhs,rhs ,tempyt  ,mgutilf,tmutil,tempct  ,tempb,tempct2 ,tmutil2,extern
real(dp) :: eulereq2   ,tempcn2 ,tempcn,tempc2   ,tempc
integer :: iss

  rhs=0.

  tempyt=se(is)
  tempb=b(h)
  tempct = cguess
  tempcn=A*senn(is)
  tempc=(omega*(tempct**(-mu))+(1-omega)*(tempcn**(-mu)))**(-1/mu)
  tmutil = tempc**(1-sig)*omega*tempct**(-mu-1)/ (omega*tempct**(-mu)+(1-omega)*tempcn**(-mu))
  lhs =  tmutil


  if (tempct<0) lhs = infinit

  tempbf= tempb*tempr    + tempyt-tempct
  mgutilf=0.
  do iss=1,nss
  call linint(b,cph(:,iss),nb,tempbf,tempct2,1)
  tempcn2=A*senn(iss)
  
  tempc2=(omega*(tempct2**(-mu))+(1-omega)*(tempcn2**(-mu)))**(-1/mu)
  
  tmutil2=tempc2**(1-sig)*omega*tempct2**(-mu-1)/ (omega*tempct2**(-mu)+(1-omega)*tempcn2**(-mu))
  
  if (tempct2<0 .or. tempcn2 <0) tmutil2=infinit
	mgutilf=p(iss,is)*tmutil2 + mgutilf
  enddo


  rhs=betay*tempr*mgutilf


  if (planner .eq. 1) then
  call externalit(h,nb,nss,is,B,tempbf,shadow,cph,extern)
  rhs = rhs + extern
  endif

  eulereq2=lhs-rhs


   end function
