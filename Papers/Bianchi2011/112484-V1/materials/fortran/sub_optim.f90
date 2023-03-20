 !--------------------------------------------------------
 ! solves optimization problem at each point in the grid
 !--------------------------------------------------------

  subroutine sub_optim(planner,nb,nss,is,h,b,cph,bph,shadow,plast,cnew,pnew,bnew,bbrnew,bmax,&
  & munew,xmid)
  use nrtype
  use common_par

  implicit none
  integer ,intent(in) :: is ,h ,nb ,planner ,nss
  real(dp), intent(in) :: b(nb),cph(nb,nss) ,plast(nb,nss) ,bph(nb,nss),shadow(nb,nss)
  real(dp ), intent(out) ::  cnew,pnew,bnew ,bmax,munew
  real(dp), intent(inout) :: xmid
  integer, intent(out) :: bbrnew
  real(dp) ::  tempbf,tempyt   ,tempb  ,eulereq2 ,zbrent   ,tempm
  real(dp) ::  x1 ,x2,cmax,tempyn
  logical(lgt) succes

  succes = .true.
  munew=0.
  bbrnew =0

  tempb=b(h)
  tempbf=bph(h,is)
  tempm=0.

  tempyt=  se(is)

  tempyn=a*senn(is)
  bmax=-(kappat*tempyt+ kappan*tempyn*plast(h,is))

  if (bmax .ge. b(1)) then

  cmax=tempyt  + tempb*tempr  - bmax        ! calculate consumption assuming constraint binds
  munew=  eulereq2(nb,nss,is,h,b,cph,cmax,planner,shadow)


  if (munew .gt. 0) then      ! binding constraint
  bnew=bmax
  bbrnew=1
  cnew=tempyt  + tempb*tempr  - bnew
  pnew=(1-omega)/omega*((cnew/tempyn)**(mu+1))
  xmid=cnew
  return
  endif

  elseif    (bmax .le. b(1)) then

  cmax=tempyt  + tempb*tempr  - b(1)
  munew=  eulereq2(nb,nss,is,h,b,cph,cmax,planner,shadow)


  if (munew .gt. 0.) then      !  constraint   binds
  bnew=b(1)
  bbrnew=1
  cnew=tempyt  + tempb*tempr  - bnew
  pnew=(1-omega)/omega*(cnew/tempyn)**(mu+1)
  return

  endif

  endif

  munew=0.

  if (h .eq. 1  ) then
  x1=cph(h,is)
  else
  x1=xmid
  endif

  x2 =x1 + 0.05

  call zbrac(eulereq2,x1,x2,succes,nb,nss,is,h,b,cph,planner,shadow)
  xmid=  zbrent(eulereq2,x1,x2,tolzbrent,nb,nss,is,h,b,cph,planner,shadow)

  cnew=xmid
  bnew=tempb*tempr +tempyt-cnew  

  if (bnew .gt. b(nb))  then
  bnew= b(nb )
  cnew=tempyt  + tempb*tempr  - bnew
  munew=  eulereq2(nb,nss,is,h,b,cph,cnew,planner,shadow)
  endif


  pnew=(1-omega)/omega*(cnew/tempyn)**(mu+1)

  return

  end subroutine
