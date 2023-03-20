module storedata

contains


subroutine writefile(pnext,cpnext,bpnext,nb,nss,B,v,bbr,borrlimit)
  use nrtype
  use common_par
  implicit none
  integer, intent(in) :: nb,nss
  real(dp), dimension(nb,nss), intent(in)::  pnext ,cpnext,bpnext,v ,borrlimit
  integer, intent(in) :: bbr(nb,nss)
  real(dp),intent(in) :: B(nb)

  integer :: h,is
  real(dp), dimension(nb,nss) :: gdpnt  ,ca,y  ,cay
  real(dp) :: tempcn,tempc


  do is=1,nss
  do h=1,nb


  tempcn=A*senn(is)
  gdpnt(h,is)=tempcn*pnext(h,is)
  y(h,is)= se(is) + gdpnt(h,is)
  ca(h,is)=bpnext(h,is)-b(h)
  cay(h,is)=ca(h,is)/y(h,is)

  tempc=(omega*(cpnext(h,is)**(-mu))+(1-omega)*(tempcn**(-mu)))**(-1/mu)

  write(1,*)  pnext(h,is)
  write(3,*)  bpnext(h,is)
  write(40,*) v(h,is)
  write(9,*) bbr(h,is)
  write(8,*)  cpnext(h,is)
  write(5,*)  borrlimit(h,is)


  enddo
  enddo


   end subroutine  writefile



  subroutine openfile(planner)
     implicit none
  integer, intent(in) :: planner



  if (planner .eq. 0) then

  print*,'-------------------------------------------------------------'
  print*,'----------- compute DE -----------'
  print*,'-------------------------------------------------------------'


  open(unit=1,file='output\equil\guesspn.txt')
  open(unit=2,file='output\equil\guessm.txt')
  open(unit=3,file='output\equil\guessp.txt')
  open(unit=4,file='output\equil\guessmu.txt')
  open(unit=5,file='output\equil\blimit.txt')
  open(unit=8,file='output\equil\guessc.txt')
  open(unit=9,file='output\equil\constraint.txt')
   open(unit=21,file='output\equil\simul.txt')
  open(unit=23,file='output\equil\probsb.txt')
  open(unit=22,file='output\equil\suddens.txt')
  open(unit=40,file='output\equil\valuef.txt')
  
  elseif (planner .eq. 1) then


  print*,'-------------------------------------------------------------'
  print*,'-----------compute planner VFI ------------'
  print*,'-------------------------------------------------------------'

  open(unit=1,file='output\planner\guesspn.txt')
  open(unit=2,file='output\planner\guessm.txt')
  open(unit=3,file='output\planner\guessp.txt')
  open(unit=4,file='output\planner\guessmu.txt')
  open(unit=5,file='output\planner\blimit.txt')
  open(unit=8,file='output\planner\guessc.txt')
  open(unit=9,file='output\planner\constraint.txt')
   open(unit=21,file='output\planner\simul.txt')
  open(unit=23,file='output\planner\probsb.txt')
  open(unit=22,file='output\planner\suddens.txt')
  open(unit=40,file='output\planner\valuef.txt')

   endif

  end subroutine

     subroutine closefile

  close(unit=1)
  close(unit=2)
  close(unit=4)
  close(unit=8)
  close(unit=5)
  close(unit=40)
  close(unit=9)
  close(unit=23)
  close(unit=9)

  end subroutine closefile

    end module storedata

