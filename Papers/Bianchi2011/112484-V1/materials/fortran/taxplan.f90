subroutine taxplan(planner,nb,nss,b,bph,cph,plast,shadow,v)

  use common_par;  use nrtype    ; use storedata

  implicit none

  integer, intent(in) :: nb  ,planner ,nss
  real (dp), intent(in) :: b(nb)
  real (dp), intent(inout) :: bph(nb,nss),plast(nb,nss),cph(nb,nss),shadow(nb,nss)
    real (dp), intent(out) ::  v(nb,nss)

  integer is ,iter ,h, bbr(nb,nss)   , bbrnew
  real(dp) bpnext(nb,nss), borrlimit(nb,nss)  ,cpnext(nb,nss)  ,pnext(nb,nss)  ,taxx(nb,nss)
  real(dp) tempdif  ,bnew ,cnew   ,tempdifp  , casdp,bm
  real(dp)   bmax ,pnew,munew,munext(nb,nss)     ,tempd(4) ,temp  ,taxtotal
  real(dp) ::tempdifmu   ,m(nb,nss)  ,xmid ,tempdifc



  bpnext=bph
  cpnext=cph
  munext=shadow
  pnext =plast


  bbr=0
  xmid=0.1
  temp=10

  do iter=1,maxitpolicy

  if (temp .lt. tolpol )  exit

  bph = lambda*bpnext +(1-lambda)*bph
  cph = lambda*cpnext +(1-lambda)*cph
  plast = lambda*pnext + (1-lambda)*plast
  shadow= lambda*munext + (1-lambda)*shadow

  do  is=1,nss
  do  h=1,nb

  call sub_optim(planner,nb,nss,is,h,b,cph,bph,shadow,plast,cnew,pnew,bnew,bbrnew,bmax,&
  &munew,xmid)

  borrlimit(h,is) = bmax
  bpnext(h,is) = bnew
  cpnext(h,is) = cnew
  bbr(h,is) = bbrnew
  munext(h,is) = munew
  pnext(h,is)=pnew

  enddo
  enddo

  tempdifmu = maxval(abs(shadow-munext)/(1.+abs(shadow)))
  tempdifc = maxval(abs(cph-cpnext)/(1.+abs(cph)))
  tempdifp = maxval(abs(plast-pnext)/(1.+abs(plast)))
  tempdif = maxval(abs(bph-bpnext)/(1.+abs(bph)))

  tempd(1)=tempdifc
  tempd(2)=tempdifp
  tempd(3)=tempdifmu
  tempd(4)=tempdif

  temp=maxval(tempd)

  enddo

  if(temp .gt. tolpol ) then
  write (*,'(a35, x,i3, x, a24, x,e10.3, x, a23)'),'no convergence', iter-1, 'iterations', temp, 'norm'
  elseif (temp .lt. tolpol ) then
  write (*,'(a23, x, i3, a24)'), ' algorithm converged'
  endif

  do is=1,nss
  do h=1,nb
  call tax(h,nb,nss,is,B,bph(h,is),shadow,cph,taxtotal)
  if (bbr(h,is) .eq. 0) then
  taxx(h,is)=taxtotal
  else
  taxx(h,is) =0.
  endif

  write(41,*) taxx(h,is)

  enddo
  enddo

  


  end subroutine
