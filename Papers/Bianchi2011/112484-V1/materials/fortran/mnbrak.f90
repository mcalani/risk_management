	SUBROUTINE mnbrak(ax,bx,cx,fa,fb,fc,func,h,is,nb,nss,V,BB)
	USE nrtype; USE nrutil, ONLY : swap
  	use common_par
	IMPLICIT NONE
			integer, intent(in) :: h,is , nb,nss
	real(dp), intent(in) :: v(nb,nss),Bb(nb)
	REAL(SP), INTENT(INOUT) :: ax,bx
	REAL(SP), INTENT(OUT) :: cx,fa,fb,fc
	INTERFACE
		FUNCTION func(x,h,is,nb,nss,V,BB)
		USE nrtype
		use common_par
		IMPLICIT NONE
    		integer, intent(in) :: h,is , nb,nss
	real(dp), intent(in) :: v(nb,nss),Bb(nb)
    REAL(SP), INTENT(IN) :: x
		REAL(SP) :: func
		END FUNCTION func
	END INTERFACE
	REAL(SP), PARAMETER :: GOLD=0.0618034_sp,GLIMIT=100.0_sp,TINY=1.0e-20_sp
	REAL(SP) :: fu,q,r,u,ulim
	fa=func(ax,h,is,nb,nss,V,BB)
	fb=func(bx,h,is,nb,nss,V,BB)

	if (fb > fa) then
		call swap(ax,bx)
		call swap(fa,fb)
	end if
	if (bx .gt. ax) then
	cx=bx+GOLD*(bx-ax)
	else
	cx=ax+GOLD*(bx-ax)
	endif
	fc=func(cx,h,is,nb,nss,V,BB)

	do

    if (fc .ne. fc)  fc=10000

    if (fb < fc) RETURN
		r=(bx-ax)*(fb-fc)
		q=(bx-cx)*(fb-fa)
		u=bx-((bx-cx)*q-(bx-ax)*r)/(2.0_sp*sign(max(abs(q-r),TINY),q-r))
		ulim=bx+GLIMIT*(cx-bx)
		if ((bx-u)*(u-cx) > 0.0) then
			fu=func(u,h,is,nb,nss,V,BB)
			if (fu < fc) then
				ax=bx
				fa=fb
				bx=u
				fb=fu
				RETURN
			else if (fu > fb) then
				cx=u
				fc=fu
				RETURN
			end if
			u=cx+GOLD*(cx-bx)
			fu=func(u,h,is,nb,nss,V,BB)
		else if ((cx-u)*(u-ulim) > 0.0) then
			fu=func(u,h,is,nb,nss,V,BB)
			if (fu < fc) then
				bx=cx
				cx=u
				u=cx+GOLD*(cx-bx)
				call shft(fb,fc,fu,func(u,h,is,nb,nss,V,BB))
			end if
		else if ((u-ulim)*(ulim-cx) >= 0.0) then
			u=ulim
			fu=func(u,h,is,nb,nss,V,BB)
		else
			u=cx+GOLD*(cx-bx)
			fu=func(u,h,is,nb,nss,V,BB)
		end if
		call shft(ax,bx,cx,u)
		call shft(fa,fb,fc,fu)
	end do
	CONTAINS

	SUBROUTINE shft(a,b,c,d)
	REAL(SP), INTENT(OUT) :: a
	REAL(SP), INTENT(INOUT) :: b,c
	REAL(SP), INTENT(IN) :: d
	a=b
	b=c
	c=d
	END SUBROUTINE shft
	END SUBROUTINE mnbrak
