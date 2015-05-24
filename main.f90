program main
	use subprogs
	implicit none

	integer, parameter :: file_output = 200, file_poten = 300
	integer i, j

	open(file_output, file = "output.d")
	open(file_poten, file = "poten.d")
	call proc1 ! input

	do j = 0, n_grid
		y = grid(2, j)
		call set(r, y, x1, x2, info)

		if(info /= 0 .or. x1 >= x2) then 
! 			call proc2 ! make
			call proc3 ! solve
		end if 

		call proc4 ! plot
	end do 
	close(file_output)
	close(file_poten)
end program main


