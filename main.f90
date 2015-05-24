program main
	use subprogs
	implicit none

	integer, parameter :: file_test1 = 101, file_test2 = 102, file_test3 = 103

	complex(8) tmp
! 	complex(16) tmp
	integer i, j



	call proc1 ! input
	call proc2 ! make

! 		test 1
	open(file_test1, file = "test1.d")
		write(file_test1, *) 
		write(file_test1, *) "eq_A: "
		write(file_test1, *)
		write(file_test1, '(9999i54, "  ")')	(/ (j, j = 1, size(eq_A(1, :))) /)
! 		write(file_test1, '(9999i26, "  ")')	(/ (j, j = 1, size(eq_A(1, :))) /)
! 		write(file_test1, *) ((abs(eq_A(i, j)), j = 1, n_eq), i = 1, n_eq)
		do i = 1, n_eq
			write(file_test1, *) (eq_A(i, j), j = 1, n_eq)
		enddo
		write(file_test1, *)
		write(file_test1, *)
		write(file_test1, *) "eq_x: "
		write(file_test1, *)
! 		write(file_test1, *) (abs(eq_x(j)), j = 1, n_eq)
		write(file_test1, *) (eq_x(j), j = 1, n_eq)
		write(file_test1, *)
		write(file_test1, *)
		write(file_test1, *) "vec_input: "
		write(file_test1, *)
! 		write(file_test1, *) (abs(eq_b(j)), j = 1, n_eq)
		write(file_test1, *) (eq_b(j), j = 1, n_eq)
		write(file_test1, *)
	close(file_test1)

	call proc3 ! solve

! 		! test 2
! 	open(file_test2, file = "test2.d")
! 		write(file_test2, '(9999i26, "  ")')	(/ (j, j = 1, size(mat_xk(1, :))) /)
! 		write(file_test2, *)
! 		write(file_test2, *) "vec_coeff: "
! 		write(file_test2, *)
! 		write(file_test2, *) (abs(vec_coeff(j)), j = 1, size(vec_coeff(:)))
! 		write(file_test2, *)
! 		write(file_test2, *)
! 		write(file_test2, *) "vec_input: "
! 		write(file_test2, *)
! 		write(file_test2, *) (abs(vec_input(j)), j = 1, size(vec_input(:)))
! 	close(file_test2)		

		! test 3 
	do i = 1, n_eq 

		tmp = 0.d0
		do j = 1, n_eq
			tmp = tmp +eq_A(i, j)*eq_x(j)
		enddo 
		tmp = tmp -eq_b(i)

		write(*, *) "TEST 3: ", i, abs(tmp), abs(eq_b(i))
	enddo

	call proc4 ! plot

end program main


