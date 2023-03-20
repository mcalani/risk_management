  subroutine plannervf(nb,nss,b,v)
  use nrtype;  use common_par       ; use storedata
  implicit none

  integer, intent(in) :: nb,nss
  real(dp), intent(in) :: b(nb)
  real(dp), intent(inout) ::  v(nb,nss)

  real(dp) :: plast(nb,nss) ,taxx(nb,nss)  ,borrlimit(nb,nss) ,metric,&
  & bpnext(nb,nss) ,maxutil,brent ,shadow(nb,nss),  bph(nb,nss), cph(nb,nss) ,shadoww(nb,nss)

  integer h,is ,iter   ,bbr(nb,nss),planner      ,binding

  real(dp) ::  va(nb,nss)       ,tol
  real(dp) ::vnew,xmin   ,tempcn

  real(dp)  :: ax,bx,cx,fa,fb,fc ,temppn,blimit    ,tempb,bnew ,cnew ,bm ,casdp ,tempc
  real(dp)  difbind,bnew0 ,tempdifmu,tempdifc,tempdif ,tempdifp   ,diff (nb,nss) ,v0(nb,nss)
  integer j

  planner=1
  call openfile(planner)

  open(unit=41,file='output\planner\tax.txt')
  open(unit=42,file='output\planner\ktight.txt')



  metric=10.


            
  do  iter = 1,maxit

  if (mod(iter,100) ==0 .or. iter .eq. 2 ) print*,'metric',metric

  if (  metric < tolvfi ) exit

  bbr=0

  do is=1,nss
  do  h=1,nb

  ax=0.4
  bx=2.

  call mnbrak(ax,bx,cx,fa,fb,fc,maxutil,h,is,nb,nss,V,B)

  vnew=brent(ax,bx,cx,maxutil,tolbrent,xmin,h,is,nb,nss,V,B)
                          !
  cnew=xmin
  va(h,is)= -vnew

  tempb=b(h)
  tempcn=senn(is)*A
  temppn=(1-omega)/omega*(cnew/tempcn)**(mu+1)
  bnew=   se(is)-cnew+tempb*tempr
  blimit=-(kappat*se(is)+ kappan*temppn*tempcn)
  difbind=10.

  if (bnew .lt. blimit) then
  do j=1,1000
  if (difbind .le. 1.e-10 ) exit
  difbind=ABS(bnew0-bnew)
  bnew0=bnew
  bnew=   se(is)-cnew+tempb*tempr
  blimit=-(kappat*se(is)+ kappan* tempcn*temppn)
  bnew=max(bnew,blimit)
  cnew=  se(is)-bnew+tempb*tempr
  tempcn=A*senn(is)
  temppn=(1-omega)/omega*(cnew/tempcn)**(mu+1)
  bnew=   se(is)-cnew+tempb*tempr
  blimit=-(kappat*se(is)+ kappan*temppn*tempcn)
  enddo
  vnew=maxutil(cnew,h,is,nb,nss,V,B)
  bbr(h,is)=1
  va(h,is)= -vnew
  planner=1

  else
  bbr(h,is)=0
  endif

  if (bnew .lt. B(1) .and. bnew .gt. blimit ) then
  bnew=  B(1)
  cnew=  se(is)-bnew+tempb*tempr
  tempcn=A*senn(is)
  temppn=(1-omega)/omega*(cnew/tempcn)**(mu+1)
  vnew=maxutil(cnew,h,is,nb,nss,V,B)
  bbr(h,is)=1
  va(h,is)= -vnew
  endif

  plast(h,is)=temppn
  cph(h,is)=cnew
  bph(h,is)=bnew

  enddo
  enddo

  diff=abs(va-v)
  metric=maxval(diff)
  v=va

  enddo   ! this finishes VFI

  if ( metric .gt. tolvfi) then
  write (*,'(a35, x,i3, x, a24, x,e10.3, x, a23)'),'WARNING: VF dont converge', iter-1, 'iterations',metric, 'norm'
  else
  write (*,'(a23, x, i3, a24,e10.3)'), 'VF converge in ',iter-1,'iterations'
  endif

  do is=1,nss
  do h=1,nb
  borrlimit(h,is)=-( kappat*se(is) + kappan*A*senn(is)*plast(h,is))
  write(42,*) bph(h,is)/borrlimit(h,is)
  enddo
  enddo

  shadoww=0.
  binding=1

  call writefile(plast,cph,bph,nb,nss,B,v,bbr,borrlimit)
  call simulation(v,b,nb,nss,bph,bm,casdp)
  call sub_timeseries(nb,nss,b,bph,borrlimit,planner,bm, casdp)
  call  closefile
  call taxplan(planner,nb,nss,b,bph,cph,plast,shadoww,v)         ! PLANNER EULER EQUATION TAX

    close(unit=41)
  close(unit=42)


  end subroutine

  function maxutil(tempct,h,is,nb,nss,V,B)
  use nrtype
  use common_par
  implicit none

  integer, intent(in) :: h,is,nb,nss
  real(dp), intent(in) :: tempct,v(nb,nss),b(nb)
  integer :: iss
  real(dp) :: tempcn,tempc,tempbf,vv ,expvv0   ,maxutil  ,util,tempb

  tempb=b(h)
  tempcn=A*senn(is)
  tempc=omega*tempct**(-mu)+(1-omega)*tempcn**(-mu)
  tempc=tempc**(-1/mu)
  util=(1./(1.-sig))*tempc**(1.-sig)

  if (tempc .le. 0)  util=-10000000.

  tempbf=se(is)-tempct+tempb*tempr
  expvv0=0.
  do iss=1,nss
  call linint(b,v(:,iss),nb,tempbf,vv,1)
  expvv0=p(iss,is)*vv+expvv0
  enddo

  maxutil=   util+betay*expvv0
  maxutil=-maxutil

  end function
