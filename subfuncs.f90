module subfuncs
	use const
! 	use globals
	implicit none

contains








! ===== function : convert coeff number=============================================

subroutine set(r, h, x1, x2, info)
	real(8), intent(in) :: r, h
	real(8), intent(out) :: x1, x2
	integer, intent(out) :: info
	real(8) tmp 

	if(abs(h) >= r) then 
		x1 = 0.d0
		x2 = 0.d0 
		info = 0

	else 
		tmp = (r**2.d0 -h**2.d0)**0.5d0
		x1 = -tmp 
		x2 =  tmp 
		info = 1

	end if
end subroutine set





! ===== function : basis =============================================

complex(8) function basis_1(f, dir, k, x)

	! f type 0: f, 1: dx f
	! direc type 1, 2

	integer, intent(in) :: f, dir
	complex(8), intent(in) :: k
	real(8), intent(in) :: x
	real(8) :: tmp, sg, kk

	if(f /= 0 .and. f /= 1) stop "ERROR function basis_1: check f value is 0 or 1"

	sg = 3.d0 -2.d0*dble(dir)
	if(abs(sg) /= 1.d0) stop "ERROR function basis_1: check dir value is 1 or 2"

	tmp = real(k**2.d0)
	kk = aimag(k**2.d0)
! 	write(*, *) kk
! 	if(kk /= 0.d0) stop "ERROR function basis_1: check k value"

	kk = abs(k)

	if(f == 0) then 
		if(tmp < 0.d0) then 
			basis_1 = exp(sg*mt_i *kk*x)

		else if(tmp == 0.d0) then 
			if(dir == 1) then 
				basis_1 = x
			else if(tmp == 2) then 
				basis_1 = 1.d0
			end if

		else if(tmp > 0.d0) then 
			basis_1 = exp(-sg     *kk*x)

		end if

	else if(f == 1) then 
		if(tmp < 0.d0) then 
			basis_1 = exp(sg*mt_i *kk*x) &
						 *sg*mt_i *kk

		else if(tmp == 0.d0) then 
			if(dir == 1) then 
				basis_1 = 1.d0
			else if(tmp == 2) then 
				basis_1 = 0.d0
			end if

		else if(tmp > 0.d0) then 
			basis_1 = exp(-sg     *kk*x) &
						*(-sg)    *kk
		end if
	end if
end function basis_1








! ===== function : equation d/dx y(x) = f(x,y)=============================================

real(8) function poten(c, r, h, x)
	real(8), intent(in) :: c, r, h, x

	poten = c/(x**2.d0 +h**2.d0)**0.5d0
! 	poten = 0.d0 
end function poten




end module subfuncs









