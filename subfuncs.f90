module subfuncs
	use const
! 	use globals
	implicit none

contains










! ! ===== function : xy_to_polar =============================================

! real(8) function to_theta(input_x, input_y)
! 	real(8), intent(in) :: input_x, input_y
! 	to_theta = atan2(input_y, input_x)
! 	if(to_theta < 0.d0) to_theta = to_theta +2.d0 * mt_pi
! end function to_theta



! real(8) function to_r(input_x, input_y)
! 	real(8), intent(in) :: input_x, input_y
! 	to_r = abs(input_x +mt_i * input_y)
! end function to_r



! real(8) function to_x(input_r, input_theta)
! 	real(8), intent(in) :: input_r, input_theta
! 	to_x = input_r * cos(input_theta)
! end function to_x



! real(8) function to_y(input_r, input_theta)
! 	real(8), intent(in) :: input_r, input_theta
! 	to_y = input_r * sin(input_theta)
! end function to_y










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



! complex(8) function basis_xy(type, input_k1, input_k2, input_x1, input_x2)

! 		! type 0: f, 1:dx f, 2:dy

! 	integer, intent(in) :: type
! 	complex(8), intent(in) :: input_k1, input_k2
! 	real(8), intent(in) :: input_x1, input_x2

! 	basis_xy = exp(input_k1*input_x1)*exp(input_k2*input_x2)

! 	if(type == 0) then 
! 	else if(type == 1) then 
! 		basis_xy = input_k1*basis_xy
! 	else if(type == 2) then 
! 		basis_xy = input_k2*basis_xy
! 	else 
! 		stop "ERROR function basis_xy: check type value"
! 	end if
! end function basis_xy



! complex(8) function basis_polar(type, basis, input_l, input_k, input_r, input_theta)

! 		! type 0: f, 1:dr f, 2:dtheta f
! 		! real 1: I, 2: K
! 		! img  1: J, 2: Y

! 	integer, intent(in) :: type, basis, input_l
! 	complex(8), intent(in) :: input_k
! 	real(8), intent(in) :: input_r, input_theta

! ! 	external bessel_Jn, bessel_Yn, bessel_In, bessel_Kn
! 	real(8) f_r
! 	complex(8) f_theta

! 	if(type == 0 .or. type == 2) then 
! 		if(real(input_k**2.d0) < 0.d0) then 
! 			if(basis == 1) then 
! 				f_r = bessel_Jn(input_l, abs(input_k) * input_r)
! 			else if(basis == 2) then 
! 				f_r = bessel_Yn(input_l, abs(input_k) * input_r)
! 			else 
! 				stop "ERROR function basis_polar: check basis value"
! 			end if 

! 		else if(real(input_k**2.d0) == 0.d0) then 
! 			stop "ERROR function basis_polar: zero"

! 		else if(real(input_k**2.d0) > 0.d0) then 
! 			stop "ERROR function basis_polar: zero"
! ! 			if(basis == 1) then 
! ! 				f_r = bessel_In(input_l, abs(input_k) * input_r)
! ! 			else if(basis == 2) then 
! ! 				f_r = bessel_Kn(input_l, abs(input_k) * input_r)
! ! 			else 
! ! 				stop "ERROR function basis_polar: check basis value"
! ! 			end if 

! 		end if 
! 	else if(type == 1) then 
! 		if(real(input_k**2.d0) < 0.d0) then 
! 			if(basis == 1) then 
! 				f_r = abs(input_k) * &
! 						0.5d0 * (bessel_Jn(input_l -1, abs(input_k) * input_r) &
! 									-bessel_Jn(input_l +1, abs(input_k) * input_r))
! 			else if(basis == 2) then 
! 				f_r = abs(input_k) * &
! 						0.5d0 * (bessel_Yn(input_l -1, abs(input_k) * input_r) &
! 									-bessel_Yn(input_l +1, abs(input_k) * input_r))
! 			else 
! 				stop "ERROR function basis_polar: check basis value"
! 			end if 

! 		else if(real(input_k**2.d0) == 0.d0) then 
! 			stop "ERROR function basis_polar: zero"

! 		else if(real(input_k**2.d0) > 0.d0) then 
! 			stop "ERROR function basis_polar: zero"
! ! 			if(basis == 1) then 
! ! 				f_r = 0.5d0 * (bessel_In(input_l -1, abs(input_k) * input_r) &
! ! 									-bessel_In(input_l +1, abs(input_k) * input_r))
! ! 			else if(basis == 2) then 
! ! 				f_r = 0.5d0 * (bessel_Kn(input_l -1, abs(input_k) * input_r) &
! ! 									-bessel_Kn(input_l +1, abs(input_k) * input_r))
! ! 			else
! ! 				stop "ERROR function basis_polar: check basis value"
! ! 			end if 

! 		end if 
! 	else 
! 		stop "ERROR function basis_polar: check type value"
! 	end if

! 	f_theta = exp(mt_i * dble(input_l) * input_theta)
! 	if(type == 0) then 
! 	else if(type == 1) then 
! 	else if(type == 2) then 
! 		f_theta = mt_i * dble(input_l) * f_theta 
! 	else 
! 		stop "ERROR function basis_polar: check type value"
! 	end if

! 	basis_polar = f_r * f_theta

! end function basis_polar








! ! ===== function : convert coeff number=============================================

subroutine set(r, y, x1, x2, info)
	real(8), intent(in) :: r, y
	real(8), intent(out) :: x1, x2
	integer, intent(out) :: info
	real(8) tmp 

	if(abs(y) >= r) then 
		x1 = 0.d0
		x2 = 0.d0 
		info = 0

	else 
		tmp = (r**2.d0 -y**2.d0)**0.5d0
		x1 = -tmp 
		x2 =  tmp 
		info = 1

	end if
end subroutine set



! integer function x_num(n_l, type, l_num, basis)

! 		! type 0: inner, 1: ref, 2: trans
! 		! l_num, basis 

! 	integer, intent(in) :: n_l, type, l_num, basis
! 	integer n_total

! 	if(l_num < -n_l .or. l_num > n_l) stop "ERROR function x_num: check n_l, l_num"
! 	if(type < 0 .or. type > 2) stop "ERROR function x_num: check type"
! 	if(basis < 1 .or. basis > 2) stop "ERROR function x_num: check basis"

! ! 	n_total = 2 * n_l +1
! 	n_total = 2 * n_l
! 	x_num = l_num 

! 	if(l_num > 0) then 
! 		x_num = l_num 
! 	else if(l_num == 0) then 
! 		stop "ERROR function x_num: check l_num"
! 	else if(l_num < 0) then 
! 		x_num = n_l -l_num
! 	end if

! 	x_num = 2 * n_total * type +n_total * (basis -1) +x_num 

! end function x_num



! integer function b_num(n_r, n_theta, type, region, num)

! 		! type 0: f, 1: dx f
! 		! region 1: circle1, 2: circle2, 3: outer1, 4: outer2

! 	integer, intent(in) :: n_r, n_theta, type, region, num
! 	integer n_total

! 	b_num = num 

! 	if(region == 1 .or. region == 2) then 
! 		if(num < 1 .or. num > n_theta) stop "ERROR function b_num: check num"
! 		b_num = b_num +n_theta * (region -1)

! 	else if(region == 3 .or. region == 4) then 
! 		if(num < 1 .or. num > n_r) stop "ERROR function b_num: check num"
! 		b_num = b_num +n_theta * 2 +n_r * (region -3)

! 	else 
! 		stop "ERROR function b_num: check region"
! 	end if

! 	if(type == 0) then 
! 	else if(type == 1) then 
! 		b_num = b_num +n_r * 2 +n_theta * 2
! 	else 
! 		stop "ERROR function b_num: check type"
! 	end if

! end function b_num











! ! ===== subroutine : solve eq =============================================

! subroutine gauss_jordan(input_a, output_x, input_b, n)
! 	integer, intent(in) :: n
! 	complex(8), intent(in) :: input_a(1:n_max, 1:n_max), input_b(1:n_max)
! ! 	complex(16), intent(in) :: input_a(1:n_max, 1:n_max), input_b(1:n_max)
! 	complex(8), intent(out) :: output_x(1:n_max)
! ! 	complex(16), intent(out) :: output_x(1:n_max)

! 	integer i, k, m
! 	real(8) am
! 	complex(8) ar, t, a(1:n_max, 1:n_max), w(n_max)
! ! 	complex(16) ar, t, a(1:n_max, 1:n_max), w(n_max)



! 	a(:, :) = input_a(:, :)
! 	output_x(:) = input_b(:)

! 	do k = 1, n

! 		m = k 
! 		am = abs(a(k, k))
! 		do i = k +1, n
! 			if(abs(a(i, k)) > am) then
! 				am = abs(a(i, k))
! 				m = i
! 			end if
! 		end do

! 		if(am == 0.d0) stop "A is singular"
! 		if(k /= m) then 
! 			w(k:n) = a(k, k:n)
! 			a(k, k:n) = a(m, k:n)
! 			a(m, k:n) = w(k:n)
! 			t = output_x(k)
! 			output_x(k) = output_x(m)
! 			output_x(m) = t
! 		end if

! 		ar = 1.d0/a(k, k)
! 		a(k, k) = 1.d0
! 		a(k, k +1:n) = ar*a(k, k +1:n)
! 		output_x(k) = ar*output_x(k)

! 		do i = 1, n
! 			if(i /= k) then 
! 				a(i, k +1:n) = a(i, k +1:n) -a(i, k)*a(k, k +1:n)
! 				output_x(i) = output_x(i) -a(i, k)*output_x(k)
! 				a(i, k) = 0.d0
! 			end if
! 		end do

! 	end do

! end subroutine gauss_jordan



subroutine solve(n, A, x, b)
	integer, intent(in) :: n
	complex(8), intent(in) :: A(1:, 1:), b(1:)
	complex(8), intent(out) :: x(1:)

	complex(8), allocatable :: tmp(:, :)
	integer, allocatable :: ipiv(:)
	integer info
	external zgesv

	complex(8) check



	allocate(tmp(1:n, 1:n), ipiv(1:n))
	x(1:) = b(1:)
	tmp(1:, 1:) = A(1:n, 1:n)



	check = det(n, A)
	write(*, *) "RUN subroutine solve: mat A determinant is ", abs(check)
! 	if(check == 0.d0) stop "ERROR external subroutine: not able to solve eq"



	call zgesv(n, 1, tmp(1:, 1:), n, ipiv(1:), x(1:), size(x(1:)), info)
	if(info == 0) then 
		write(*, *) "RUN subroutine solve: solved"
	else 
		stop "ERROR subroutine solve: failure"
	endif
end subroutine solve



complex(8) function det(n, mat) 
	integer, intent(in) :: n
	complex(8), intent(in) :: mat(1:, 1:)

	integer i, info, count
	integer, allocatable :: ipiv(:)
	complex(8), allocatable :: tmp(:, :)
	real(8) sgn

	external zgetrf 



	allocate(ipiv(1:n), tmp(1:n, 1:n))
	ipiv = 0.d0
	tmp = mat

	call zgetrf(n, n, tmp, n, ipiv, info)

	det = 1.d0
	count = 0
	do i = 1, n
		if(abs(tmp(i, i)) == 0.d0) then 
			count = count +1
! 			write(*, *) i, "ERROR external subroutine: not able to solve eq"
! 			stop "ERROR external subroutine: not able to solve eq"
		end if
		det = det*tmp(i, i)
	end do
	write(*, *) count, n

	sgn = 1.d0
	do i = 1, n
		if(ipiv(i) /= i) sgn = -sgn 
	end do 

	det = sgn*det 
	if(abs(det) == 0.d0) write(*, *) "ERROR function det: det is zero"

end function det 










! ! ===== function : phase =============================================

! complex(8) function phase(n_l, n_eq, eq_x, r0, k_input, k_inner, k_outer, input_r, input_theta) 
! 	integer, intent(in) :: n_l, n_eq
! 	complex(8), intent(in) :: eq_x(1:n_eq), k_input(1:2), k_inner, k_outer
! 	real(8), intent(in) :: r0, input_r, input_theta

! 	integer type, n_total, j, bs, sg, jj

! 	n_total = 4 * n_l

! 	if(input_r <= r0) then 
! 		type = 0
! 	else if(input_theta < 0.5d0 * mt_pi .and. input_theta >= 1.5d0 * mt_pi) then 
! 		type = 1
! 	else if(input_theta <= 0.d0 .and. input_theta > 2.d0 * mt_pi) then 
! 		stop "ERROR function phase: check input_theta"
! ! 	else if(input_r > r1) then 
! ! 		stop "ERROR function phase: check input_r"
! 	else 
! 		type = 2
! 	end if

! 	if(type == 0) then 
! 		phase = 0.d0
! 		do bs = 1, 2
! 		do sg = 1, -1, -2
! 		do jj = 1, n_l

! 			j = sg * jj

! 			phase = phase &
! 				+eq_x(x_num(n_l, 0, j, bs)) &
! 					*basis_polar(0, bs, j, k_inner, input_r, input_theta)

! 		end do 
! 		end do 
! 		end do 

! 	else if(type == 1) then 
! 		phase = basis_xy(0, k_input(1), k_input(2), to_x(input_r, input_theta), to_y(input_r, input_theta))
! 		phase = 0.d0
! 		do bs = 1, 2
! 		do sg = 1, -1, -2
! 		do jj = 1, n_l

! 			j = sg * jj

! 			phase = phase &
! 				+eq_x(x_num(n_l, 1, j, bs)) &
! 					*basis_polar(0, bs, j, k_inner, input_r, input_theta)

! 		end do 
! 		end do 
! 		end do 

! 	else if(type == 2) then 
! 		phase = 0.d0
! 		phase = 0.d0
! 		do bs = 1, 2
! 		do sg = 1, -1, -2
! 		do jj = 1, n_l

! 			j = sg * jj

! 			phase = phase &
! 				+eq_x(x_num(n_l, 2, j, bs)) &
! 					*basis_polar(0, bs, j, k_inner, input_r, input_theta)

! 		end do 
! 		end do 
! 		end do 
! 	end if
! end function phase 

end module subfuncs









