module subprogs
	use globals
	use subfuncs
	implicit none

contains










! ===== process 1: make globals =============================================

subroutine proc1
	write(*, *) "PROCESS 1"
	call input
end subroutine proc1









! ===== process 2: setting =============================================

subroutine proc2
	integer i, j, f



		! make matrix A
	write(*, *) "PROCESS 2"
	mat_xk(:, :) = 0.d0 
	do f = 0, 1
		mat_xk(2*f +1, 1) = -basis_1(f, 2, k_outer, x1)
! 		mat_xk(2*f +1, 2) =  0.d0
		mat_xk(2*f +1, 3) =  basis_1(f, 1, k_inner, x1)
		mat_xk(2*f +1, 4) =  basis_1(f, 2, k_inner, x1)

! 		mat_xk(2*f +2, 1) =  0.d0
		mat_xk(2*f +2, 2) = -basis_1(f, 1, k_outer, x2)
		mat_xk(2*f +2, 3) =  basis_1(f, 1, k_inner, x2)
		mat_xk(2*f +2, 4) =  basis_1(f, 2, k_inner, x2)
	end do 

		! make vector x
	write(*, *) "PROCESS 2: setting vector x"
	vec_coeff(:) = 0.d0



		! make vector b
	write(*, *) "PROCESS 2: setting vector b"
	vec_input(:) = 0.d0
	vec_input(1) = basis_1(0, 1, k_outer, x1)
	vec_input(3) = basis_1(1, 1, k_outer, x1)
	
end subroutine proc2










! ===== process 3: solve eq =============================================

subroutine proc3
	write(*, *) "PROCESS 3: slove Ax = b"
	call solve(n_total, mat_xk(1:, 1:), vec_coeff(1:), vec_input(1:))
end subroutine proc3









! ===== process 4: plot =============================================

subroutine proc4
	integer :: file_output = 200

	complex(16) tmp 
	integer i, j



	write(*, *) "PROCESS 4"
! 	open(file_output, file = "output.d")

	do i = 0, n_grid

		if(info == 0 .or. x1 >= x2) then 
			tmp = basis_1(0, 1, k_outer, grid(1, i))
			write(file_output, '(99(E19.12, 1X))') grid(1, i)*au_aa, y*au_aa, abs(tmp)**2.d0 

		else if(grid(1, i) < x1) then 
			tmp = basis_1(0, 1, k_outer, grid(1, i))
			tmp = tmp +vec_coeff(1)*basis_1(0, 2, k_outer, grid(1, i))
			write(file_output, '(99(E19.12, 1X))') grid(1, i)*au_aa, y*au_aa, abs(tmp)**2.d0 

		else if(grid(1, i) >= x1 .and. grid(1, i) <= x2) then 
			tmp = 0.d0 
			tmp = tmp +vec_coeff(3)*basis_1(0, 1, k_inner, grid(1, i))
			tmp = tmp +vec_coeff(4)*basis_1(0, 2, k_inner, grid(1, i))
			write(file_output, '(99(E19.12, 1X))') grid(1, i)*au_aa, y*au_aa, abs(tmp)**2.d0 

		else if(grid(1, i) > x2) then 
			tmp = 0.d0 
			tmp = tmp +vec_coeff(2)*basis_1(0, 1, k_outer, grid(1, i))
			write(file_output, '(99(E19.12, 1X))') grid(1, i)*au_aa, y*au_aa, abs(tmp)**2.d0 

! 		else 
! 			stop "something problem"
		end if
	end do 
	write(file_output, '(99(E19.12, 1X))')

end subroutine proc4
end module subprogs




