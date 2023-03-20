subroutine bindingc(planner,nb,nss,b,v)

  use common_par;  use nrtype  ; use storedata

  implicit none

  integer, intent(in) :: nb  ,planner ,nss
  real (dp), intent(in) :: b(nb)
  real (dp), intent(out) ::  v(nb,nss)
  
  real (dp) :: bph(nb,nss),plast(nb,nss),cph(nb,nss),shadow(nb,nss)
  integer is ,iter ,h, bbr(nb,nss)   , bbrnew ,hnfs
  real(dp) bpnext(nb,nss), borrlimit(nb,nss)  ,cpnext(nb,nss)   ,pnext(nb,nss)
  real(dp) tempdif  ,bnew ,cnew   ,tempdifp  , casdp,bm
  real(dp)   bmax ,pnew,munew,munext(nb,nss)     ,tempd(4) ,temp
  real(dp) ::tempdifmu   ,m(nb,nss)  ,xmid ,tempdifc ,tempcn


   call   openfile(planner)
   

     do  is=1,nss
     do  h=1,nb
     cph(h,is)=se(is)
     tempcn= A*senn(is)
     bph(h,is)=  cph(h,is)  + b(h)*tempr  - cph(h,is)
     plast(h,is)=(1-omega)/omega * (cph(h,is)/tempcn)**(mu+1)
     enddo
     enddo

  bpnext=bph
  cpnext=cph
  munext=shadow
  pnext =plast
  hnfs=0

  tempdif=10
  bbr=0
  tempdifp=10
  xmid=0.1
  temp=10


 munext=shadow

  do iter=1,maxitpolicy

  if (temp .lt. tolpol )  exit

  if (mod(iter,50) ==0 .or. iter .eq. 2 ) print*,'metric',temp

  bph = lambda*bpnext +(1-lambda)*bph
  cph = lambda*cpnext +(1-lambda)*cph
  plast = lambda*pnext + (1-lambda)*plast
  shadow= lambda*munext + (1-lambda)*shadow


  do  is=1,nss
  do  h=1,nb

  call sub_optim(planner,nb,nss,is,h,b,cph,bph,shadow,plast,cnew,pnew,bnew,bbrnew,bmax,munew,xmid)

  borrlimit(h,is) = bmax
  bpnext(h,is) = bnew
  cpnext(h,is) = cnew
  bbr(h,is) = bbrnew
  pnext(h,is)=pnew
  munext(h,is)=munew

  enddo
  enddo

  tempdif = maxval(abs(bph-bpnext)/(1.+abs(bph)))
  tempdifmu = maxval(abs(shadow-munext)/(1.+abs(shadow)))
  tempdifc = maxval(abs(cph-cpnext)/(1.+abs(cph)))
  tempdifp = maxval(abs(plast-pnext)/(1.+abs(plast)))

  tempd(1)=tempdif
  tempd(2)=tempdifmu
  tempd(3)=tempdifc
  tempd(4)=tempdifp

  temp=maxval(tempd)


  enddo


  if(temp .gt. tolpol ) then
  write (*,'(a35, x,i3, x, a24, x,e10.3, x, a23)'),'WARNING: dec rules dont converge', iter-1, 'iterations', temp, 'norm'
  print*, "location",maxloc(abs(bph-bpnext)/(1+abs(bph)))
  stop
  elseif (tempdif .lt. tolpol ) then
  write (*,'(a23, x, i3, a24)'), 'dec rules converge in ',iter-1,'iterations'
  endif

  call  valuef(bpnext,b,v,nb,nss)
  call writefile(pnext,cpnext,bpnext,nb,nss,B,v,bbr,borrlimit)
  call simulation(v,b,nb,nss,bph,bm,casdp)
  call sub_timeseries(nb,nss,b,bph,borrlimit,planner,bm ,casdp)
  call closefile



  end subroutine bindingc



