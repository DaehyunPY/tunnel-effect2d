program main
	use subprogs
	implicit none

	integer, parameter :: file_test1 = 101, file_test2 = 102, file_test3 = 103

	complex(8) tmp
! 	complex(16) tmp
	integer i, j



	call proc1 ! input
	call proc2 ! make

! ! 		test 1
! 	open(file_test1, file = "test1.d")
! 		write(file_test1, *) 
! 		write(file_test1, *) "mat_xk: "
! 		write(file_test1, *)
! ! 		write(file_test1, '(9999i54, "  ")')	(/ (j, j = 1, size(mat_xk(1, :))) /)
! ! 		write(file_test1, '(9999i26, "  ")')	(/ (j, j = 1, size(mat_xk(1, :))) /)
! ! 		write(file_test1, *) ((abs(mat_xk(i, j)), j = 1, n_total), i = 1, n_total)
! 		do i = 1, n_total
! 			write(file_test1, *) (mat_xk(i, j), j = 1, n_total)
! 		enddo
! 		write(file_test1, *)
! 		write(file_test1, *)
! 		write(file_test1, *) "vec_coeff: "
! 		write(file_test1, *)
! ! 		write(file_test1, *) (abs(vec_coeff(j)), j = 1, n_total)
! 		write(file_test1, *) (vec_coeff(j), j = 1, n_total)
! 		write(file_test1, *)
! 		write(file_test1, *)
! 		write(file_test1, *) "vec_input: "
! 		write(file_test1, *)
! ! 		write(file_test1, *) (abs(vec_input(j)), j = 1, n_total)
! 		write(file_test1, *) (vec_input(j), j = 1, n_total)
! 		write(file_test1, *)
! 	close(file_test1)

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
	do i = 1, n_total 

		tmp = 0.d0
		do j = 1, n_total
			tmp = tmp +mat_xk(i, j)*vec_coeff(j)
		enddo 
		tmp = tmp -vec_input(i)

		write(*, *) "TEST 3: ", i, abs(tmp), abs(vec_input(i))
	enddo

	call proc4 ! standardization


! 		! test 3 
! 	do i = 1, n_total 

! 		tmp = 0.d0
! 		do j = 1, n_total
! 			tmp = tmp +mat_xk(i, j)*vec_coeff(j)
! 		enddo 
! 		tmp = tmp -vec_input(i)

! 		write(*, *) "TEST 3: ", i, abs(tmp), abs(vec_input(i))
! 	enddo

	call proc5 ! plot

end program main


