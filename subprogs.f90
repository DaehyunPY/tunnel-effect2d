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

	do i = 1, n_x
		do f = 1, 2
		do j = 1, n_inner

				!        f,         x type,                        k type (inner)
			mat_xk(x_num(f, n_x, n_y, i, 1), k_num(n_inner, n_outer, j, 1)) = &
				+phase(f, k_inner(1, j), k_inner(2, j), x1(1, i), x1(2, i), 1)

			mat_xk(x_num(f, n_x, n_y, i, 2), k_num(n_inner, n_outer, j, 1)) = &
				+phase(f, k_inner(1, j), k_inner(2, j), x2(1, i), x2(2, i), 2)

			mat_xk(x_num(f, n_x, n_y, i, 3), k_num(n_inner, n_outer, j, 1)) = &
				+phase(f, k_inner(1, j), k_inner(2, j), x3(1, i), x3(2, i), 3)

			mat_xk(x_num(f, n_x, n_y, i, 4), k_num(n_inner, n_outer, j, 1)) = &
				+phase(f, k_inner(1, j), k_inner(2, j), x4(1, i), x4(2, i), 4)

! 			mat_xk(x_num(f, n_x, n_y, i, 5), k_num(n_inner, n_outer, j, 1)) = &
! 				+phase(f, k_inner(1, j), k_inner(2, j), y1(1, i), y1(2, i), 5)

! 			mat_xk(x_num(f, n_x, n_y, i, 6), k_num(n_inner, n_outer, j, 1)) = &
! 				+phase(f, k_inner(1, j), k_inner(2, j), y2(1, i), y2(2, i), 6)

		end do
		end do 
	end do 
	do i = 1, n_x
		do f = 1, 2
		do j = 1, n_outer 

				!        f,         x type,                        k type (ref)
			mat_xk(x_num(f, n_x, n_y, i, 1), k_num(n_inner, n_outer, j, 2)) = &
				-phase(f, k_ref(1, j), k_ref(2, j), x1(1, i), x1(2, i), 1)

! 			mat_xk(x_num(f, n_x, n_y, i, 2), k_num(n_inner, n_outer, j, 2)) = &
! 				-phase(f, k_ref(1, j), k_ref(2, j), x2(1, i), x2(2, i), 2)

! 			mat_xk(x_num(f, n_x, n_y, i, 3), k_num(n_inner, n_outer, j, 2)) = &
! 				-phase(f, k_ref(1, j), k_ref(2, j), x3(1, i), x3(2, i), 3)

			mat_xk(x_num(f, n_x, n_y, i, 4), k_num(n_inner, n_outer, j, 2)) = &
				-phase(f, k_ref(1, j), k_ref(2, j), x4(1, i), x4(2, i), 4)

		end do 
		end do 
	end do 
	do i = 1, n_y
		do f = 1, 2
		do j = 1, n_outer

			mat_xk(x_num(f, n_x, n_y, i, 5), k_num(n_inner, n_outer, j, 2)) = &
				-phase(f, k_ref(1, j), k_ref(2, j), y1(1, i), y1(2, i), 5)

			mat_xk(x_num(f, n_x, n_y, i, 6), k_num(n_inner, n_outer, j, 2)) = &
				-phase(f, k_ref(1, j), k_ref(2, j), y2(1, i), y2(2, i), 6)

		end do 
		end do 
	end do 
	do i = 1, n_x
		do f = 1, 2
		do j = 1, n_outer

				!        f,         x type,                        k type (trans)
! 			mat_xk(x_num(f, n_x, n_y, i, 1), k_num(n_inner, n_outer, j, 3)) = &
! 				+phase(f, k_trans(1, j), k_trans(2, j), x1(1, i), x1(2, i), 1)

			mat_xk(x_num(f, n_x, n_y, i, 2), k_num(n_inner, n_outer, j, 3)) = &
				-phase(f, k_trans(1, j), k_trans(2, j), x2(1, i), x2(2, i), 2)

			mat_xk(x_num(f, n_x, n_y, i, 3), k_num(n_inner, n_outer, j, 3)) = &
				-phase(f, k_trans(1, j), k_trans(2, j), x3(1, i), x3(2, i), 3)

! 			mat_xk(x_num(f, n_x, n_y, i, 4), k_num(n_inner, n_outer, j, 3)) = &
! 				+phase(f, k_trans(1, j), k_trans(2, j), x4(1, i), x4(2, i), 4)

		end do 
		end do 
	end do 
	do i = 1, n_y
		do f = 1, 2
		do j = 1, n_outer

			mat_xk(x_num(f, n_x, n_y, i, 5), k_num(n_inner, n_outer, j, 3)) = &
				+phase(f, k_trans(1, j), k_trans(2, j), y1(1, i), y1(2, i), 5)

			mat_xk(x_num(f, n_x, n_y, i, 6), k_num(n_inner, n_outer, j, 3)) = &
				+phase(f, k_trans(1, j), k_trans(2, j), y2(1, i), y2(2, i), 6)

		end do
		end do 
	end do



		! make vector x
	write(*, *) "PROCESS 2: setting vector x"
	vec_coeff(:) = 0.d0



		! make vector b
	write(*, *) "PROCESS 2: setting vector b"
	vec_input(:) = 0.d0

	do i = 1, n_x
		do f = 1, 2 

				!           f,           x type, k type (inner)
			vec_input(x_num(f, n_x, n_y, i, 1)) = &
				+phase(f, k_input(1), k_input(2), x1(1, i), x1(2, i), 1)

! 			vec_input(x_num(f, n_x, n_y, i, 2)) = &
! 				+phase(f, k_input(1), k_input(2), x2(1, i), x2(2, i), 2)

! 			vec_input(x_num(f, n_x, n_y, i, 3)) = &
! 				+phase(f, k_input(1), k_input(2), x3(1, i), x3(2, i), 3)

			vec_input(x_num(f, n_x, n_y, i, 4)) = &
				+phase(f, k_input(1), k_input(2), x4(1, i), x4(2, i), 4)

		end do 
	end do 
	do i =1, n_y
		do f = 1, 2

			vec_input(x_num(f, n_x, n_y, i, 5)) = &
				+phase(f, k_input(1), k_input(2), y1(1, i), y1(2, i), 5)

			vec_input(x_num(f, n_x, n_y, i, 6)) = &
				+phase(f, k_input(1), k_input(2), y2(1, i), y2(2, i), 6)

		end do 
	end do 

end subroutine proc2










! ===== process 3: solve eq =============================================

subroutine proc3
	write(*, *) "PROCESS 3: slove Ax = b"
	call solve(n_total, mat_xk(1:, 1:), vec_coeff(1:), vec_input(1:))
end subroutine proc3









! ===== process 4: standardization =============================================

subroutine proc4
	real(8) tmp, tmp1, tmp2
	integer i

		! inner
	tmp = 0.d0
	do i = 1, n_inner
		tmp = tmp +abs(vec_coeff(k_num(n_inner, n_outer, i, 1)))**2.d0
	end do 
	write(*, *) "inner", tmp**0.5d0

	do i = 1, n_inner
		vec_coeff(k_num(n_inner, n_outer, i, 1)) = vec_coeff(k_num(n_inner, n_outer, i, 1))/tmp**0.5d0
	end do 

		! ref & trans
	tmp1 = 0.d0
	tmp2 = 0.d0
	tmp = 0.d0
	do i = 1, n_outer
		tmp1 = tmp1 +abs(vec_coeff(k_num(n_inner, n_outer, i, 2)))**2.d0
		tmp2 = tmp2 +abs(vec_coeff(k_num(n_inner, n_outer, i, 3)))**2.d0
	end do 
	tmp = tmp1 +tmp2

	write(*, *) "ref", tmp1**0.5d0
	write(*, *) "trans", tmp2**0.5d0
	write(*, *) "norm", tmp

! 	do i = 1, n_outer
! 		vec_coeff(k_num(n_inner, n_outer, i, 2)) = vec_coeff(k_num(n_inner, n_outer, i, 2))/tmp**0.5d0
! 		vec_coeff(k_num(n_inner, n_outer, i, 3)) = vec_coeff(k_num(n_inner, n_outer, i, 3))/tmp**0.5d0
! 	end do 
end subroutine proc4









! ===== process 5: plot =============================================

subroutine proc5
	integer :: file_output = 100

	complex(8) tmp 
	integer i, j, n, m 



	write(*, *) "PROCESS 4"
	open(file_output, file = "output.d")



		! region inner
	do n = 1, n_grid
	do m = 1, n_grid

		tmp = 0.d0
		do j = 1, n_inner
			tmp = tmp +vec_coeff(k_num(n_inner, n_outer, j, 1)) &
				*phase(1, k_inner(1, j), k_inner(2, j), grid(1, n), grid(2, m), 1)
		end do 
		write(file_output, *) grid(1, n)*au_aa, grid(2, m)*au_aa, abs(tmp)**2.d0

	end do 
	write(file_output, *)	
	end do 
	write(file_output, *)	



		! region ref
	do n = -n_grid, 0
	do m = -n_grid, 2*n_grid

		tmp = phase(1, k_input(1), k_input(2), grid(1, n), grid(2, m), 1)
! 		tmp = 0.d0
		do j = 1, n_outer
			tmp = tmp +vec_coeff(k_num(n_inner, n_outer, j, 2)) &
				+phase(1, k_ref(1, j), k_ref(2, j), grid(1, n), grid(2, m), 1)
		end do 
		write(file_output, *) grid(1, n)*au_aa, grid(2, m)*au_aa, abs(tmp)**2.d0

	end do 
	write(file_output, *)	
	end do 
	write(file_output, *)	
	do n = 1, 2*n_grid
	do m = -n_grid, 0

		tmp = phase(1, k_input(1), k_input(2), grid(1, n), grid(2, m), 1)
! 		tmp = 0.d0
		do j = 1, n_outer
			tmp = tmp +vec_coeff(k_num(n_inner, n_outer, j, 2)) &
				+phase(1, k_ref(1, j), k_ref(2, j), grid(1, n), grid(2, m), 1)
		end do 
		write(file_output, *) grid(1, n)*au_aa, grid(2, m)*au_aa, abs(tmp)**2.d0

	end do 
	write(file_output, *)	
	end do 
	write(file_output, *)	



		! region trans
	do n = 1, 2*n_grid
	do m = n_grid +1, 2*n_grid

		tmp = 0.d0
		do j = 1, n_outer
			tmp = tmp +vec_coeff(k_num(n_inner, n_outer, j, 3)) &
				*phase(1, k_trans(1, j), k_trans(2, j), grid(1, n), grid(2, m), 1)
		end do 
		write(file_output, *) grid(1, n)*au_aa, grid(2, m)*au_aa, abs(tmp)**2.d0

	end do 
	write(file_output, *)	
	end do 
	write(file_output, *)	
	do n = n_grid +1, 2*n_grid
	do m = 1, n_grid

		tmp = 0.d0
		do j = 1, n_outer
			tmp = tmp +vec_coeff(k_num(n_inner, n_outer, j, 3)) &
				*phase(1, k_trans(1, j), k_trans(2, j), grid(1, n), grid(2, m), 1)
		end do 
		write(file_output, *) grid(1, n)*au_aa, grid(2, m)*au_aa, abs(tmp)**2.d0

	end do 
	write(file_output, *)
	end do 
	write(file_output, *)	

	close(file_output)

end subroutine proc5
end module subprogs




