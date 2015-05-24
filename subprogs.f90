module subprogs
	use globals
	use subfuncs
	use sde
	implicit none

contains










! ===== process 1: make globals =============================================

subroutine proc1
	write(*, *) "PROCESS 1"
	call input
end subroutine proc1









! ===== process 2: setting =============================================

! subroutine proc2
! 	integer i, j, f



! ! 		! make matrix A
! ! 	write(*, *) "PROCESS 2"
! ! 	mat_xk(:, :) = 0.d0 
! ! 	do f = 0, 1
! ! 		mat_xk(2*f +1, 1) = -basis_1(f, 2, k_outer, x1)
! ! ! 		mat_xk(2*f +1, 2) =  0.d0
! ! 		mat_xk(2*f +1, 3) =  basis_1(f, 1, k_inner, x1)
! ! 		mat_xk(2*f +1, 4) =  basis_1(f, 2, k_inner, x1)

! ! ! 		mat_xk(2*f +2, 1) =  0.d0
! ! 		mat_xk(2*f +2, 2) = -basis_1(f, 1, k_outer, x2)
! ! 		mat_xk(2*f +2, 3) =  basis_1(f, 1, k_inner, x2)
! ! 		mat_xk(2*f +2, 4) =  basis_1(f, 2, k_inner, x2)
! ! 	end do 

! ! 		! make vector x
! ! 	write(*, *) "PROCESS 2: setting vector x"
! ! 	vec_coeff(:) = 0.d0



! ! 		! make vector b
! ! 	write(*, *) "PROCESS 2: setting vector b"
! ! 	vec_input(:) = 0.d0
! ! 	vec_input(1) = basis_1(0, 1, k_outer, x1)
! ! 	vec_input(3) = basis_1(1, 1, k_outer, x1)
	
! end subroutine proc2










! ===== process 3: solve eq =============================================

function f(input_x, input_y) result(result) ! dx y = f(x, y) y
	real(8), intent(in) :: input_x
	complex(8), intent(in) :: input_y(1:2)
	complex(8) result(1:2)

	result(1) = k_outer*input_y(2)
	result(2) = k_outer*(1.d0 - poten(b_charge, r, y, input_x)/p_energy)*input_y(1)
end function f



subroutine proc3
	real(8) width, dx
	complex(8) phase, tem(1:2)
	integer i



! 	write(*, *) "PROCESS 3: slove Ax = b"
	width = x2 -x1 	
	dx = width/dble(n_cal)
	do i = 0, n_cal 
		x(i) = x1 +dx*dble(i)
	enddo

	inner(:, :) = 0.d0
	inner(:, n_cal) = 1.d0
	refle = 0.d0 
	trans = 1.d0 

	do i = n_cal -1, 0, -1
		tem(:) = inner(:, i +1)
! 		call rk8_vec(f, dx, x(i +1), tem(:))
		call rk8_vec(f, dx, x(i +1), tem(:))
		inner(:, i) = tem(:)
	enddo

	phase = (inner(1, 0) +inner(2, 0))/2.d0
	inner(:, :) = inner(:, :)/phase
	phase = basis_1(0, 1, k_outer, x1)	

! 	input = (inner(1, 0) +inner(2, 0))/2.d0 	*exp(psi_k_outter*r_zero) 
! 	refle = (inner(1, 0) -inner(2, 0))/2.d0 	*exp(psi_k_outter*r_zero) 
! 	trans = inner(1, n_cal) 					*exp(psi_k_outter*r_zero)
! 	inner(:, :) = inner(:, :) 					*exp(psi_k_outter*r_zero)

! 	input = (inner(1, 0) +inner(2, 0))/2.d0 	*phase/phase
	refle = (inner(1, 0) -inner(2, 0))/2.d0 	*phase*phase
	trans = inner(1, n_cal) 					*phase*phase
	inner(:, :) = inner(:, :) 					*phase 

	write(*, *) abs((inner(1, 0) -inner(2, 0))/2.d0)**2.d0 +abs(inner(1, n_cal))**2.d0
! 	write(*, *) abs((inner(1, 0) -inner(2, 0))/2.d0)**2.d0, abs(inner(1, n_cal))**2.d0

end subroutine proc3









! ! ===== process 4: plot =============================================

subroutine proc4
	integer :: file_output = 200, file_poten = 300

	complex(16) tmp 
	integer i, j, count



	write(*, *) "PROCESS 4"

	do i = 0, n_grid

		if(info == 0 .or. x1 >= x2) then 
			tmp = basis_1(0, 1, k_outer, grid(1, i))
			write(file_output, '(99(E19.12, 1X))') grid(1, i)*au_aa, y*au_aa, abs(tmp)**2.d0 
! 			write(file_output, '(99(E19.12, 1X))') grid(1, i)*au_aa, y*au_aa, real(tmp)**2.d0 

		else if(grid(1, i) < x1) then 
			tmp = basis_1(0, 1, k_outer, grid(1, i))
			tmp = tmp +refle*basis_1(0, 2, k_outer, grid(1, i))
			write(file_output, '(99(E19.12, 1X))') grid(1, i)*au_aa, y*au_aa, abs(tmp)**2.d0 
! 			write(file_output, '(99(E19.12, 1X))') grid(1, i)*au_aa, y*au_aa, real(tmp)**2.d0 

		else if(grid(1, i) >= x1 .and. grid(1, i) <= x2) then 
		count = 0 
		do j = 0, n_cal 
			if(i /= n_grid) then 
				if(grid(1, i) <= x(j) .and. grid(1, i +1) >= x(j) .and. count == 0) then 
					count = count +1
					tmp = inner(1, j)
					write(file_output, '(99(E19.12, 1X))') x(j)*au_aa, y*au_aa, abs(tmp)**2.d0 
! 					write(file_output, '(99(E19.12, 1X))') x(j)*au_aa, y*au_aa, real(tmp)**2.d0 
				end if
			end if
		end do 

		else if(grid(1, i) > x2) then 
			tmp = 0.d0 
			tmp = tmp +trans*basis_1(0, 1, k_outer, grid(1, i))
			write(file_output, '(99(E19.12, 1X))') grid(1, i)*au_aa, y*au_aa, abs(tmp)**2.d0 
! 			write(file_output, '(99(E19.12, 1X))') grid(1, i)*au_aa, y*au_aa, real(tmp)**2.d0 

! 		else 
! 			stop "something problem"
		end if

		write(file_poten, *) grid(1, i)*au_aa, y*au_aa, poten(b_charge, r, y, grid(1, i))

	end do 
	write(file_output, '(99(E19.12, 1X))')


end subroutine proc4

end module subprogs




