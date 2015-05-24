module subfuncs
	use const
! 	use globals
	implicit none

contains










! ! ===== function : set =============================================

! complex(8) function set(input_k, input_x)

! 		! type 1: f, 2:dx f

! 	complex(8), intent(in) :: input_k
! 	real(8), intent(in) :: input_x



! 	set = exp(input_k*input_x)

! end function set



! ===== function : phase =============================================

complex(8) function phase(type, input_k1, input_k2, input_x1, input_x2, region)

		! type 1: f, 2:dx f
		! region 1: x1, 2: x2, 3: x3, 4:x4, 5:y1, 6:y2

	integer, intent(in) :: type, region
	complex(8), intent(in) :: input_k1, input_k2
	real(8), intent(in) :: input_x1, input_x2

	integer tmp 



	tmp = mod(region, 2)
! 	write(*, *) tmp
! 	phase = set(input_k1,input_x1)*set(input_k2,input_x2)
	phase = exp(input_k1*input_x1)*exp(input_k2*input_x2)

	if(type == 1) then 
	else if(type == 2) then 

		if(tmp == 0) then 
! 			phase = input_k1*phase
			phase = input_k2*phase
		else if(tmp == 1) then 
! 			phase = input_k2*phase 
			phase = input_k1*phase
		else 
			stop "ERROR function phase: check region value"
		end if

	else 
		stop "ERROR function phase: check type value"
	end if

end function phase










! ===== function : convert coeff number=============================================

integer function k_num(n_inner, n_outer, k, region)

		! region 1: inner, 2: ref, 3: trans

	integer, intent(in) :: n_inner, n_outer, k, region



	if(region == 1) then 
		k_num = k
	else if(region == 2) then 
		k_num = n_inner +k
	else if(region == 3) then 
		k_num = n_inner +n_outer +k
	else
		stop "ERROR function k_num: check region value"
	end if

end function k_num









! ===== function : convert coeff number=============================================

integer function x_num(type, n_x, n_y, x, region)

		! type 1: f, 2: dx f
		! region 1: x1, 2: x2, 3: x3, 4:x4, 5:y1, 6:y2

	integer, intent(in) :: type, n_x, n_y, x, region 



	if(region > 0 .and. region < 5) then 
		x_num = n_x*dble(region -1) +x
	else if(region > 4 .and. region < 7) then 
		x_num = 4*n_x +n_y*dble(region -5) +x
	else 
		stop "ERROR function x_num: check region value"
	end if

	if(type == 1) then 
	else if(type == 2) then 
		x_num = 4*n_x +2*n_y +x_num
	else
		stop "ERROR function x_num: check type value"
	end if

end function x_num









! ===== function : convert coeff number=============================================

integer function f_num(type, n_x, x)

		! type 1: x1, 2: x2

	integer, intent(in) :: type, n_x, x
	integer tmp



	tmp = (x -1)/n_x +1
	if(tmp == 1) then 
		f_num = 1
	else if(tmp == 2) then 
		f_num = 2
	else if(tmp == 3) then 
		f_num = 1
	else if(tmp == 4) then 
		f_num = 2
	else 
		stop "ERROR function f_num: check n_x, x value"
	end if



	if(type == 1) then 
	else if(type == 2) then 
		f_num = 3 -f_num
	else
		stop "ERROR function f_num: check type value"
	end if

end function f_num








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

end module subfuncs









