  subroutine sub_timeseries(nb,nss,b,bph,borrlimit,planner,bm,threshold1)

  use common_par
  implicit none

  integer, intent(in) :: nb ,planner     ,nss
  real (dp), intent(in) :: b(nb), bph(nb,nss), bm,borrlimit(nb,nss)    ,threshold1

  integer :: t  , state(TT)     , productivity  ,is  ,simConstraint(TT)
  real(dp) ::       transitioncdf(nss,nss)  ,pp(nss,nss)
  real(dp) :: simB(TT+1),simC(TT),simY(TT+1),simPrice(TT),simLimit(TT)   ,simyt(TT),simRER(TT) ,simCN(TT),simYN(TT),&
  &  simBu(TT+1) ,simCA(TT),simCAY(TT),shockt(TT) ,shocknt(TT)  ,simCT(TT)
  real(dp) ::     randomnumber   ,threshold

  integer, dimension(1) :: j
  integer ::      simSS(TT)

   open(unit=24,file='output\thres.txt')


  if (planner .eq. 0) then
  write(24,*) threshold1
  threshold=threshold1
  else
  read(24,*) threshold
  endif


  open(unit=23,file='output\state_simTT.txt')
  open(unit=25,file='output\state_sim.txt')

  ! repeat shock

  pp=transpose(p)
  transitioncdf(:,1)=pp(:,1)
  do is= 2,nss
  transitioncdf(:,is) = transitioncdf(:,is-1)+pp(:,is)
  enddo


  productivity=2       !initial state
  simBu(1)=bm     ! initial value of bonds
  simSS=0
  simConstraint=0

   
  do t=1,TT

  state (t)= productivity
  shockt(t)=se(state(t))
  shocknt(t)=senn(state(t))

  call linint(B,bph(:,state(t)),nb,simB(t),simBu(t+1),1)

  simYT(t) = shockt(t)
  simCN(t)=A*shocknt(t)
  simCT(t)= simB(t)*(1+rbar)-simBu(t+1) +simYT(t)
  simPrice(t)=(1-omega)/omega*(simCT(t)/simCN(t))**(1+mu)

  simLimit(t)=-(kappan*simPrice(t)*simCN(t)+kappat*simYT(t))


  if (simBu(t+1) .le. simLimit(t)+tollimit) then    !
  simConstraint(t) = 1
  simB(t+1) = simLimit(t)
  else
  simB(t+1)=simBu(t+1)
  endif


  if (simB(t+1) .le.  B(1)+0.00000001 ) then
  simConstraint(t)=1
  simB(t+1)=B(1)
  endif


  simCT(t)= simB(t)*(1+rbar)-simB(t+1) +simYT(t)
  simPrice(t)=(1-omega)/omega*(simCT(t)/simCN(t))**(1+mu)
  simRER(t)=   (  omega**(1/(1+mu))  +  (1-omega)**(1/(1+mu))*simPrice(t)**(mu/(1+mu)))  **((1+mu)/mu);
  simYN(t)= A*shocknt(t)*simPrice(t)

  simY(t)=simYT(t)+ simYN(t)
  simCA(t)=simB(t+1)-simB(t)
  simCAY(t)=simCA(t)/simY(t)
  simC(t) =(omega*(simCT(t)**(-mu))+(1-omega)*(simCN(t)**(-mu)))**(-1/mu)

  if (simConstraint(t) .eq. 1 .and. simCA(t) .gt. threshold) then
  simSS(t)=1
  endif

  
  
  if (t>tcut)  then
   write(21,*)  simB(t),simY(t),simCA(t), simCAY(t),simC(t),simRER(t) ! ,
  write(22,*) simConstraint(t),simSS(t)
  write(23,*) productivity
   endif


  if (planner .eq. 0 ) then
  call random_number(randomnumber)
  j = minloc(transitioncdf(productivity,:), mask=transitioncdf(productivity,:)>randomnumber)         ! simulating markov chain
  productivity=j(1)
  write(25,*) productivity
  else
  read(25,*) productivity
  endif

  enddo


  close(unit=21)
  close(unit=22)
  close(unit=23)
  close(unit=24)
  close(unit=25)

  end subroutine



