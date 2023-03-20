  subroutine simulation(vv,b,nb,nss,bph,bm,casdp)
  use nrtype
  use common_par
  implicit none

  integer, intent(in) :: nb    ,nss
  real (dp), intent(in) :: vv(nb,nss),b(nb)  ,bph(nb,nss)
  real (dp), intent(out) ::   bm  ,casdp
  integer(dp) hu,hl  ,huu(nb,nss),hll(nb,nss) , cz,mm ,h ,is,iss ,jjj
  real(dp) ::  pka1(nb,nss),pka2(nb,nss),pka(nb,nss) ,pb(nb),casd
  real(dp) ::  sss,cam  ,ca(nb,nss) ,bbp(nb) ,tempbf     ,tempb

  jjj=(nss*nb)-count(vv.le.-999.0e+09)

  do is=1,nss
  do  h=1,nb
  pka2(h,is)=0.0
  if(vv(h,is).gt.-999.0e+09) then
  pka2(h,is)=1./jjj                        !initial guess for ergodic distribution
  endif
  enddo
  enddo

  cz=0
  pka1=0.
  sss=10.0

  do  is=1,nss
  do h=1,nb
  bbp=(b-bph(h,is))
  hu=minloc(bbp,dim=1,mask=bbp.ge.0.0)
  hl=maxloc(bbp,dim=1,mask=bbp.le.0.0)
  if(bbp(nb).lt.0) then
  hu=nb
  endif
  if(bbp(1).ge.0) then
  hl=1
  endif
  huu(h,is)=hu
  hll(h,is)=hl
  enddo
  enddo

  do 290 mm=1,5001
  sss=0.0
  cz=cz+1
  do  is=1,nss
  do  h=1,nb
  hu=huu(h,is)
  hl=hll(h,is)
  do 295 iss=1,nss
  if(hu.eq.hl)        pka1(hl,iss)=pka1(hl,iss)+p(iss,is)*pka2(h,is)
  if(hu.ne.hl) then
  pka1(hl,iss)=pka1(hl,iss)+p(iss,is)*pka2(h,is)*(b(hu)-bph(h,is))/(b(hu)-b(hl))
  pka1(hu,iss)=pka1(hu,iss)+p(iss,is)*pka2(h,is)*(bph(h,is)-b(hl))/(b(hu)-b(hl))
  endif
  295 continue
  enddo
  enddo
  sss=sum(abs(pka2-pka1))
  if(sss.le.0.00000005) go to 340
  pka2=pka1
  pka1=0.0
  if(mm.le.5000) go to 290
  print *,'limiting distribution does not converge'
  290 continue
  340 pka=pka1
  if(mm.gt.2000) pka=pka2

  bm=0.0

  do h=1,nb
  pb(h)=sum(pka(h,1:nss))
  enddo

  do  is=1,nss
  do  h=1,nb
  tempb=b(h)
  tempbf=bph(h,is)
  ca(h,is)=tempbf-b(h)
  enddo
  enddo

  cam=0.

  do is=1,nss
  do h=1,nb
  cam=pka(h,is)*ca(h,is) + cam
  enddo
  enddo

  casd=0.
  do h=1,nb
  do is=1,nss

  casd= pka(h,is)*(ca(h,is)-cam)**2. + casd
  enddo
  enddo

  casdp= (casd**0.5)

  do  is=1,nss
  do  h=1,nb
  write(23,*)  pka(h,is)

  bm=pka(h,is)*b(h) + bm
  enddo
  enddo

  close(unit=23)
  end subroutine simulation



