module common_par
use nrtype
implicit none

!'-------------------------------------------------------------'
! -------------------model parameters --------------------------
!'-------------------------------------------------------------'

real(dp),parameter  ::  mu=.205,  rbar= 0.04 ,tempr=1+rbar , betay= 0.906  ,sig=2   ,A=1.
real(dp),parameter  ::   kappan=0.3235 ,  kappat=kappan  , omega= 0.307000693252802

! --------------grid  ----------------------------------------
!'-------------------------------------------------------------'

real(dp) ,parameter  :: minb= -1.02
real(dp) ,parameter  :: largestb=-0.4
integer ,	parameter  :: nb1=800
real(dp) ,parameter  :: gsp=1.1

!'-------------------------------------------------------------'
! --------------numerical  ------------------------------
!'-------------------------------------------------------------'

integer ,parameter  ::  maxitpolicy=1000
real(dp) ,parameter  :: tolpol=1.0e-7
real(dp),parameter  ::  tolvfi=1.000e-7
integer ,	parameter ::  ntry   =200
integer ,	parameter ::  maxitbis=2000
real(dp), parameter ::  xacc=1.e-10
real(dp),parameter ::   infinit= HUGE(0.)
integer, parameter  ::  maxit=1000
real(dp),parameter  ::  lambda=1.
real(dp),parameter  ::  tolzbrent=1.0e-8
 real(dp),parameter ::  tolbrent=1e-8
!'-------------------------------------------------------------'
! --------------stochastic process ------------------------------
!'-------------------------------------------------------------'

integer, parameter  ::  nip=4 ! no. of T shocks  and nt shocks
integer, parameter  ::  nipp=nip ! no. of Nt shocks       (must equal nvaln defined below)
integer , parameter  :: nss1=nip*nipp

!'-------------------------------------------------------------'
! --------------others ------------------------------
!'-------------------------------------------------------------'

integer, parameter  ::  tcut=1000           ! observation deleted
integer, parameter  ::  tt=80000      ! number of simulations
real(dp),parameter  :: tollimit =   0.00001
! global
real(dp) :: senn(nss1), se(nss1),p(nss1,nss1)
real(dp) :: ktight(nb1,nss1)


end module common_par


module global

 ! supporting module for tauchen hussey

  implicit none

      integer zznvar,zznlag,zzns,zzthet,zzpmat,zziw,zznst,zziw2, &
     &   zziw3,zzw  ,nvaln ,nlag,nrquad ,nvar
      parameter(nvaln=4)         ! NUMBER OF DISCRETE POINT FOR EACH VARIABLE
      parameter(nvar=2)
      parameter(nlag=1)
      parameter(nrquad=210)
      parameter(zznvar=2)
      parameter(zznlag=1)
      parameter(zzns=nvaln*nvaln)
      parameter(zzthet=zznvar+(zznvar**2)*(zznlag+2)+1)
      parameter(zzpmat=(zzns**zznlag)*zzns)
      parameter(zziw=2*zzns+2*zznlag)
      parameter(zznst=zzns**zznlag)
      parameter(zzw=4*zzns+2*zznvar+4*zzns*zznvar+2*zznvar*zznlag)
      parameter(zziw2=zznst*zznlag)
      parameter(zziw3=zznst*zzns)

end module   global


