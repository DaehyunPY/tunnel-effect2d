module subfuncs
	use const
! 	use globals
	implicit none

contains










! ===== function : xy_to_polar =============================================

real(8) function to_theta(input_x, input_y)
	real(8), intent(in) :: input_x, input_y
	to_theta = atan2(input_y, input_x)
	if(to_theta < 0.d0) to_theta = to_theta +2.d0 * mt_pi
end function to_theta



real(8) function to_r(input_x, input_y)
	real(8), intent(in) :: input_x, input_y
	to_r = abs(input_x +mt_i * input_y)
end function to_r



real(8) function to_x(input_r, input_theta)
	real(8), intent(in) :: input_r, input_theta
	to_x = input_r * cos(input_theta)
end function to_x



real(8) function to_y(input_r, input_theta)
	real(8), intent(in) :: input_r, input_theta
	to_y = input_r * sin(input_theta)
end function to_y










! ===== function : basis_xy =============================================

complex(8) function basis_xy(type, input_k1, input_k2, input_x1, input_x2)

		! type 0: f, 1:dx f, 2:dy

	integer, intent(in) :: type
	complex(8), intent(in) :: input_k1, input_k2
	real(8), intent(in) :: input_x1, input_x2

	basis_xy = exp(input_k1*input_x1)*exp(input_k2*input_x2)

	if(type == 0) then 
	else if(type == 1) then 
		basis_xy = input_k1*basis_xy
	else if(type == 2) then 
		basis_xy = input_k2*basis_xy
	else 
		stop "ERROR function basis_xy: check type value"
	end if
end function basis_xy



complex(8) function basis_polar(type, basis, input_l, input_k, input_r, input_theta)

		! type 0: f, 1:dr f, 2:dtheta f
		! real 1: +I, 2: -I, 3: +K, 4: -K
		! img  1: +J, 2: -J, 3: +Y, 4: -Y

	integer, intent(in) :: type, basis, input_l
	complex(8), intent(in) :: input_k
	real(8), intent(in) :: input_r, input_theta

! 	external bessel_Jn, bessel_Yn, bessel_In, bessel_Kn
	real(8) f_r
	complex(8) f_theta

	if(type == 0 .or. type == 2) then 
		if(real(input_k**2.d0) < 0.d0) then 
			if(basis == 1 .or. basis == 2) then 
				f_r = bessel_Jn(input_l, abs(input_k) * input_r)
			else if(basis == 3 .or. basis == 4) then 
				f_r = bessel_Yn(input_l, abs(input_k) * input_r)
			else 
				stop "ERROR function basis_polar: check basis value"
			end if 

		else if(real(input_k**2.d0) == 0.d0) then 
			stop "ERROR function basis_polar: zero"

		else if(real(input_k**2.d0) > 0.d0) then 
			stop "ERROR function basis_polar: zero"
! 			if(basis == 1) then 
! 				f_r = bessel_In(input_l, abs(input_k) * input_r)
! 			else if(basis == 2) then 
! 				f_r = bessel_Kn(input_l, abs(input_k) * input_r)
! 			else 
! 				stop "ERROR function basis_polar: check basis value"
! 			end if 

		end if 
	else if(type == 1) then 
		if(real(input_k**2.d0) < 0.d0) then 
			if(basis == 1 .or. basis == 2) then 
				f_r = abs(input_k) * &
						0.5d0 * (bessel_Jn(input_l -1, abs(input_k) * input_r) &
									-bessel_Jn(input_l +1, abs(input_k) * input_r))
			else if(basis == 3 .or. basis == 4) then 
				f_r = abs(input_k) * &
						0.5d0 * (bessel_Yn(input_l -1, abs(input_k) * input_r) &
									-bessel_Yn(input_l +1, abs(input_k) * input_r))
			else 
				stop "ERROR function basis_polar: check basis value"
			end if 

		else if(real(input_k**2.d0) == 0.d0) then 
			stop "ERROR function basis_polar: zero"

		else if(real(input_k**2.d0) > 0.d0) then 
			stop "ERROR function basis_polar: zero"
! 			if(basis == 1) then 
! 				f_r = 0.5d0 * (bessel_In(input_l -1, abs(input_k) * input_r) &
! 									-bessel_In(input_l +1, abs(input_k) * input_r))
! 			else if(basis == 2) then 
! 				f_r = 0.5d0 * (bessel_Kn(input_l -1, abs(input_k) * input_r) &
! 									-bessel_Kn(input_l +1, abs(input_k) * input_r))
! 			else
! 				stop "ERROR function basis_polar: check basis value"
! 			end if 

		end if 
	else 
		stop "ERROR function basis_polar: check type value"
	end if

	if(type == 0 .or. type == 1) then 
		if(basis == 1 .or. basis == 3) then 
			f_theta = exp( mt_i * dble(input_l) * input_theta)
		else if(basis == 2 .or. basis == 4) then 
			f_theta = exp(-mt_i * dble(input_l) * input_theta)
		else 
			stop "ERROR function basis_polar: check	basis value"
		end if 

	else if(type == 2) then 
		if(basis == 1 .or. basis == 3) then 
			f_theta =  mt_i * dble(input_l) * exp( mt_i * dble(input_l) * input_theta)
		else if(basis == 2 .or. basis == 4) then 
			f_theta = -mt_i * dble(input_l) * exp(-mt_i * dble(input_l) * input_theta)
		else 
			stop "ERROR function basis_polar: check	basis value"
		end if 
		
	else 
		stop "ERROR function basis_polar: check type value"
	end if

	basis_polar = f_r * f_theta

end function basis_polar








! ! ===== function : convert coeff number=============================================

! integer function k_num(n_inner, n_outer, k, region)

! 		! region 1: inner, 2: ref, 3: trans

! 	integer, intent(in) :: n_inner, n_outer, k, region



! 	if(region == 1) then 
! 		k_num = k
! 	else if(region == 2) then 
! 		k_num = n_inner +k
! 	else if(region == 3) then 
! 		k_num = n_inner +n_outer +k
! 	else
! 		stop "ERROR function k_num: check region value"
! 	end if

! end function k_num









! ===== function : convert coeff number=============================================

integer function x_num(n_l, type, l_num, basis)

		! type 0: inner, 1: ref, 2: trans
		! l_num, basis 

	integer, intent(in) :: n_l, type, l_num, basis
	integer n_total

	if(l_num < -n_l .or. l_num > n_l) stop "ERROR function x_num: check n_l, l_num"
	if(type < 0 .or. type > 2) stop "ERROR function x_num: check type"
	if(basis < 1 .or. basis > 4) stop "ERROR function x_num: check basis"

! 	n_total = 2 * n_l +1
	n_total = 2 * n_l
	x_num = l_num 

	if(l_num > 0) then 
		x_num = l_num 
	else if(l_num == 0) then 
		stop "ERROR function x_num: check l_num"
	else if(l_num < 0) then 
		x_num = n_l -l_num
	end if

	x_num = 4 * n_total * type +n_total * (basis -1) +x_num 

end function x_num



integer function b_num(n_r, n_theta, type, region, num)

		! type 0: f, 1: dx f
		! region 1: circle1, 2: circle2, 3: outer1, 4: outer2

	integer, intent(in) :: n_r, n_theta, type, region, num
	integer n_total

	b_num = num 

	if(region == 1 .or. region == 2) then 
		if(num < 1 .or. num > n_theta) stop "ERROR function b_num: check num"
		b_num = b_num +n_theta * (region -1)

	else if(region == 3 .or. region == 4) then 
		if(num < 1 .or. num > n_r) stop "ERROR function b_num: check num"
		b_num = b_num +n_theta * 2 +n_r * (region -3)

	else 
		stop "ERROR function b_num: check region"
	end if

	if(type == 0) then 
	else if(type == 1) then 
		b_num = b_num +n_r * 2 +n_theta * 2
	else 
		stop "ERROR function b_num: check type"
	end if

end function b_num










! ! ===== function : convert coeff number=============================================

! integer function x_num(type, n_x, n_y, x, region)

! 		! type 1: f, 2: dx f
! 		! region 1: x1, 2: x2, 3: x3, 4:x4, 5:y1, 6:y2

! 	integer, intent(in) :: type, n_x, n_y, x, region 

! 	if(region > 0 .and. region < 5) then 
! 		x_num = n_x*dble(region -1) +x
! 	else if(region > 4 .and. region < 7) then 
! 		x_num = 4*n_x +n_y*dble(region -5) +x
! 	else 
! 		stop "ERROR function x_num: check region value"
! 	end if

! 	if(type == 1) then 
! 	else if(type == 2) then 
! 		x_num = 4*n_x +2*n_y +x_num
! 	else if(type == 3) then 
! 		x_num = 2*(4*n_x +2*n_y) +x_num
! 	else
! 		stop "ERROR function x_num: check type value"
! 	end if

! end function x_num



! integer function f_num(type, n_x, x)

! 		! type 1: x1, 2: x2

! 	integer, intent(in) :: type, n_x, x
! 	integer tmp

! 	tmp = (x -1)/n_x +1
! 	if(tmp == 1) then 
! 		f_num = 1
! 	else if(tmp == 2) then 
! 		f_num = 2
! 	else if(tmp == 3) then 
! 		f_num = 1
! 	else if(tmp == 4) then 
! 		f_num = 2
! 	else 
! 		stop "ERROR function f_num: check n_x, x value"
! 	end if

! 	if(type == 1) then 
! 	else if(type == 2) then 
! 		f_num = 3 -f_num
! 	else
! 		stop "ERROR function f_num: check type value"
! 	end if

! end function f_num








! ! ===== subroutine : solve eq by gauss jordan method =============================================

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










! ===== subroutine : solve eq by lapack =============================================

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










! ===== function : calculate determinant =============================================

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









! ===== function : phase =============================================

complex(8) function phase(n_l, n_eq, eq_x, r0, k_input, k_inner, k_outer, input_r, input_theta) 
	integer, intent(in) :: n_l, n_eq
	complex(8), intent(in) :: eq_x(1:n_eq), k_input(1:2), k_inner, k_outer
	real(8), intent(in) :: r0, input_r, input_theta

	integer type, n_total, j, bs, sg, jj

	n_total = 4 * (2 * n_l)

	if(input_r <= r0) then 
		type = 0
	else if(input_theta < 0.5d0 * mt_pi .and. input_theta >= 1.5d0 * mt_pi) then 
		type = 1
	else if(input_theta <= 0.d0 .and. input_theta > 2.d0 * mt_pi) then 
		stop "ERROR function phase: check input_theta"
! 	else if(input_r > r1) then 
! 		stop "ERROR function phase: check input_r"
	else 
		type = 2
	end if

	if(type == 0) then 
		phase = 0.d0
		do bs = 1, 2
		do sg = 1, -1, -2
		do jj = 1, n_l

			j = sg * jj

			phase = phase &
				+eq_x(x_num(n_l, 0, j, bs)) &
					*basis_polar(0, bs, j, k_inner, input_r, input_theta)

		end do 
		end do 
		end do 

	else if(type == 1) then 
		phase = basis_xy(0, k_input(1), k_input(2), to_x(input_r, input_theta), to_y(input_r, input_theta))
		phase = 0.d0
		do bs = 1, 2
		do sg = 1, -1, -2
		do jj = 1, n_l

			j = sg * jj

			phase = phase &
				+eq_x(x_num(n_l, 1, j, bs)) &
					*basis_polar(0, bs, j, k_inner, input_r, input_theta)

		end do 
		end do 
		end do 

	else if(type == 2) then 
		phase = 0.d0
		phase = 0.d0
		do bs = 1, 2
		do sg = 1, -1, -2
		do jj = 1, n_l

			j = sg * jj

			phase = phase &
				+eq_x(x_num(n_l, 2, j, bs)) &
					*basis_polar(0, bs, j, k_inner, input_r, input_theta)

		end do 
		end do 
		end do 
	end if
end function phase 

end module subfuncs









