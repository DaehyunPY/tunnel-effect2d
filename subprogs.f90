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
	integer i, j, f, rg, bs, sg, jj 



		! matrix A
	eq_A = 0.d0

	do f = 0, 1
	do i = 1, n_theta
		do bs = 1, 4
		do sg = 1, -1, -2
		do jj = 1, n_l 

			j = sg * jj

				! circle1(1) & inner(0): +
			eq_A(b_num(n_r, n_theta, f, 1, i), x_num(n_l, 0, j, bs)) = &
				 basis_polar(f, bs, j, k_inner, r(0), 0.5d0 * mt_pi +theta(i))

				! circle1(1) & ref(1): -
			eq_A(b_num(n_r, n_theta, f, 1, i), x_num(n_l, 1, j, bs)) = &
				-basis_polar(f, bs, j, k_outer, r(0), 0.5d0 * mt_pi +theta(i))

				! circle1(2) & inner(0): +
			eq_A(b_num(n_r, n_theta, f, 2, i), x_num(n_l, 0, j, bs)) = &
				 basis_polar(f, bs, j, k_inner, r(0), 1.5d0 * mt_pi +theta(i))

				! circle1(2) & trans(2): -
			eq_A(b_num(n_r, n_theta, f, 2, i), x_num(n_l, 2, j, bs)) = &
				-basis_polar(f, bs, j, k_outer, r(0), 1.5d0 * mt_pi +theta(i))

		end do 
		end do 
		end do 
	end do 
	end do 

	do i = 1, n_r 
		do bs = 1, 4
		do sg = 1, -1, -2
		do jj = 1, n_l 

			j = sg * jj

				! f(0) & outer1(3) & trans(2): +
			eq_A(b_num(n_r, n_theta, 0, 3, i), x_num(n_l, 2, j, bs)) = &
				 basis_polar(0, bs, j, k_outer, r(i), 0.5d0 * mt_pi)

				! f(0) & outer1(3) & ref(1): -
			eq_A(b_num(n_r, n_theta, 0, 3, i), x_num(n_l, 1, j, bs)) = &
				-basis_polar(0, bs, j, k_outer, r(i), 0.5d0 * mt_pi)

				! f(0) & outer2(4) & trans(2): +
			eq_A(b_num(n_r, n_theta, 0, 4, i), x_num(n_l, 2, j, bs)) = &
				 basis_polar(0, bs, j, k_outer, r(i), 1.5d0 * mt_pi)

				! f(0) & outer2(4) & ref(1): -
			eq_A(b_num(n_r, n_theta, 0, 4, i), x_num(n_l, 1, j, bs)) = &
				-basis_polar(0, bs, j, k_outer, r(i), 1.5d0 * mt_pi)


				! dx f(1) & outer1(3) & trans(2): +
			eq_A(b_num(n_r, n_theta, 1, 3, i), x_num(n_l, 2, j, bs)) = &
				 cos(0.5d0 * mt_pi) &
				 	* basis_polar(1, bs, j, k_outer, r(i), 0.5d0 * mt_pi) &
				+cos(0.5d0 * mt_pi) / r(i) &
					* basis_polar(2, bs, j, k_outer, r(i), 0.5d0 * mt_pi)

				! dx f(1) & outer1(3) & ref(1): -
			eq_A(b_num(n_r, n_theta, 1, 3, i), x_num(n_l, 1, j, bs)) = &
				 cos(0.5d0 * mt_pi) &
				 	* basis_polar(1, bs, j, k_outer, r(i), 0.5d0 * mt_pi) &
				+cos(0.5d0 * mt_pi) / r(i) &
					* basis_polar(2, bs, j, k_outer, r(i), 0.5d0 * mt_pi)

				! dx f(1) & outer2(4) & trans(2): +
			eq_A(b_num(n_r, n_theta, 1, 4, i), x_num(n_l, 2, j, bs)) = &
				 cos(1.5d0 * mt_pi) &
				 	* basis_polar(1, bs, j, k_outer, r(i), 1.5d0 * mt_pi) &
				+cos(1.5d0 * mt_pi) / r(i) &
					* basis_polar(2, bs, j, k_outer, r(i), 1.5d0 * mt_pi)

				! dx f(1) & outer2(4) & ref(1): -
			eq_A(b_num(n_r, n_theta, 1, 4, i), x_num(n_l, 1, j, bs)) = &
				 cos(1.5d0 * mt_pi) &
				 	* basis_polar(1, bs, j, k_outer, r(i), 1.5d0 * mt_pi) &
				+cos(1.5d0 * mt_pi) / r(i) &
					* basis_polar(2, bs, j, k_outer, r(i), 1.5d0 * mt_pi)

		end do 
		end do 
		end do 
	end do 

		! vector x 
	eq_x = 0.d0

		! vector b
	eq_b = 0.d0

	do i = 1, n_theta
		eq_b(b_num(n_r, n_theta, 0, 1, i)) = &
			basis_xy(0, k_input(1), k_input(2), to_x(r(0), 0.5d0 * mt_pi +theta(i)), to_y(r(0), 0.5d0 * mt_pi +theta(i)))

		eq_b(b_num(n_r, n_theta, 1, 1, i)) = &
			 cos(0.5d0 * mt_pi +theta(i)) &
			 	* basis_xy(1, k_input(1), k_input(2), to_x(r(0), 0.5d0 * mt_pi +theta(i)), to_y(r(0), 0.5d0 * mt_pi +theta(i))) &
			+sin(0.5d0 * mt_pi +theta(i)) &
				* basis_xy(2, k_input(1), k_input(2), to_x(r(0), 0.5d0 * mt_pi +theta(i)), to_y(r(0), 0.5d0 * mt_pi +theta(i)))
	end do 

	do f = 0, 1
	do i = 1, n_r
		eq_b(b_num(n_r, n_theta, f, 3, i)) = & 
			basis_xy(f, k_input(1), k_input(2), to_x(r(i), 0.5d0 * mt_pi), to_y(r(i), 0.5d0 * mt_pi))
			
		eq_b(b_num(n_r, n_theta, f, 4, i)) = & 
			basis_xy(f, k_input(1), k_input(2), to_x(r(i), 1.5d0 * mt_pi), to_y(r(i), 1.5d0 * mt_pi))
	end do
	end do

end subroutine proc2










! ===== process 3: solve eq =============================================

subroutine proc3
	write(*, *) "PROCESS 3: slove Ax = b"
	call solve(n_eq, eq_A(1:, 1:), eq_x(1:), eq_b(1:))
end subroutine proc3









! ===== process 4: plot =============================================

subroutine proc4
	integer :: file_output = 100

	complex(8) tmp 
	integer i, j



	write(*, *) "PROCESS 4"
	open(file_output, file = "output.d")



	do i = 0, n_grid
	do j = 0, n_grid

		tmp = phase(n_l, n_eq, eq_x, r(0), k_input, k_inner, k_outer, grid(3, i, j), grid(4, i, j))
! 		write(*, *) "writing output file"
		write(file_output, *) grid(1, i, j)*au_aa, grid(2, i, j)*au_aa, abs(tmp)**2.d0

	end do 
	write(file_output, *)
	end do 
	write(file_output, *)

	close(file_output)

end subroutine proc4

end module subprogs




