module globals
	use const
	implicit none

	integer, parameter :: n_total = 2*12*2**8, n_grid = 40, times_y = 1
! 	integer, parameter :: n_total = 12*8, n_grid = 40
	integer, parameter :: n_inner = n_total/2, n_outer = n_total/4
	integer, parameter :: n_x = n_total/(2*(4 +2*times_y)), n_y = n_total*times_y/(2*(4 +2*times_y))
! 	real(8), parameter :: c_outer = 1.0d0
! 	real(8), parameter :: c_outer = 5.0d0

	complex(8) k_inner(1:2, 1:n_inner), k_input(1:2), k_ref(1:2, 1:n_outer), k_trans(1:2, 1:n_outer)
	real(8)    x1(1:2, 1:n_x), x2(1:2, 1:n_x), x3(1:2, 1:n_x), x4(1:2, 1:n_x), y1(1:2, 1:n_y), y2(1:2, 1:n_y)
	complex(8) mat_xk(1:n_total, 1:n_total), vec_coeff(1:n_total), vec_input(1:n_total) ! Ax = b
	real(8)    grid(-n_grid:2*n_grid, -n_grid:2*n_grid)

	real(8), parameter :: au_aa = 1.d10*au_bh
	
contains





subroutine input 
	real(8) p_mass,  p_energy, p_angle
	real(8) b_poten, b_width1, b_width2

	integer, parameter :: file_input = 100

	character dum_ch
	complex(8) inner, outer
	real(8) tmp, d_theta, d_x, d_grid1, d_grid2
	integer i



	write(*, *) "RUN subroutine input: reading"
	open(file_input, file = 'input.d')

	read(file_input, *) dum_ch				! '==================================================' 
	read(file_input, *) dum_ch				! 'PARTICLE'
	read(file_input, *) dum_ch				! '=================================================='
	read(file_input, *) dum_ch, p_mass		! 'MASS 				[au]'					1.d0
	read(file_input, *) dum_ch, p_energy	! 'ENERGY 			[eV]'					20.d0
	read(file_input, *) dum_ch, p_angle		! 'INPUT ANGLE		[o]'					60.d0
	read(file_input, *)
	read(file_input, *)
	read(file_input, *) dum_ch				! '=================================================='
	read(file_input, *) dum_ch				! 'POTENTIAL TYPE I: SQUARE'
	read(file_input, *) dum_ch				! '=================================================='
	read(file_input, *) dum_ch, b_poten		! 'POTENTIAL ENERGY 	[eV]'					30.d0
	read(file_input, *) dum_ch, b_width1	! 'WIDTH[x] 			[aa]'					10.d0
	read(file_input, *) dum_ch, b_width2	! 'WIDTH[y]			[aa]'					10.d0
	read(file_input, *)
! 	read(file_input, *)
! 	read(file_input, *) dum_ch				! '=================================================='
! 	read(file_input, *) dum_ch				! 'CALCULATE'
! 	read(file_input, *) dum_ch				! '=================================================='
! 	read(file_input, *) dum_ch, n_total		! 'CALCULATE NUMBER (POWER OF TWO)'			2**8

	close(file_input)	



! 	if(n_total > n_max) stop "ERROR subroutine input: n_total is too big"

! 	p_mass   = p_mass   /si_e_evt/sc_c**2.d0/au_m 	! eV to au
	p_energy = p_energy /si_e_evt/au_hr 			! eV to au
	p_angle  = p_angle  /360.d0*(2.d0*mt_pi) 		! do to radian
	b_poten  = b_poten  /si_e_evt/au_hr 			! eV to au	
	b_width1 = b_width1 /1.d10/au_bh 				! aa to au
	b_width2 = b_width2 /1.d10/au_bh 				! aa to au

! 	n_k = n_total/2
! 	n_x = n_total/8
	write(*, *) n_total



	write(*, *) "RUN subroutine input: setting k"

		! outer k
! 	outer = -2.d0*p_mass*p_energy
	tmp = -2.d0*p_mass*p_energy
	outer = mt_i*(-tmp)**0.5

		! outer input k
	k_input(1)  = outer*cos(p_angle)
	k_input(2)  = outer*sin(p_angle)

		! outer ref k
! 	d_theta = 1.5d0*mt_pi/dble(n_outer -1)
! 	do i = 1, n_outer
! 		k_ref(1, i) = -outer*cos(-0.5d0*mt_pi +d_theta*dble(i -1))
! 		k_ref(2, i) = -outer*sin(-0.5d0*mt_pi +d_theta*dble(i -1))
! 	end do
	d_theta = 1.5d0*mt_pi/dble(n_outer +1)
	do i = 1, n_outer
		k_ref(1, i) = -outer*cos(-0.5d0*mt_pi +d_theta*dble(i))
		k_ref(2, i) = -outer*sin(-0.5d0*mt_pi +d_theta*dble(i))
	end do

		! outer trans k
	d_theta = 0.5d0*mt_pi/dble(n_outer +1)
	do i = 1, n_outer
		k_trans(1, i) = outer*cos(d_theta*dble(i))
		k_trans(2, i) = outer*sin(d_theta*dble(i))
	end do

		! inner k 
	tmp = -2.d0*p_mass*(p_energy -b_poten)
	if(tmp < 0.d0) then 
		inner = mt_i*(-tmp)**0.5
	else if(tmp > 0.d0) then 
		inner = -tmp**0.5
	else
		stop "p_energy and b_poten is same, it is not supported"
	end if

	d_theta = 2.d0*mt_pi/dble(n_inner)
	do i = 1, n_inner
		k_inner(1, i) = inner*cos(d_theta*dble(i))
		k_inner(2, i) = inner*sin(d_theta*dble(i))
	end do
! 	write(*, *) k_inner




	write(*, *) "RUN subroutine input: setting x"

		! down side (count right)
	d_x = b_width1/dble(n_x)
	do i = 1, n_x
		x1(1, i) = -0.5d0*b_width1 +d_x*dble(i)
		x1(2, i) = -0.5d0*b_width2
	end do

		! right side (count up)
	d_x = b_width2/dble(n_x)
	do i = 1, n_x
		x2(1, i) =  0.5d0*b_width1
		x2(2, i) = -0.5d0*b_width2 +d_x*dble(i)
	end do 

		! up side (count left)
	d_x = b_width1/dble(n_x)
	do i = 1, n_x
		x3(1, i) =  0.5d0*b_width1 -d_x*dble(i)
		x3(2, i) =  0.5d0*b_width2
	end do 

		! left side (count down)
	d_x = b_width2/dble(n_x)
	do i = 1, n_x
		x4(1, i) = -0.5d0*b_width1 
		x4(2, i) =  0.5d0*b_width2 -d_x*dble(i)
	end do 

		! non box, down side (count right)
	d_x = b_width1/dble(n_x)
	do i = 1, n_y
		y1(1, i) =  0.5d0*b_width1 +d_x*dble(i)
		y1(2, i) = -0.5d0*b_width2
	end do 
! 	y1(:, :) = c_outer*y1(:, :)

		! non box, left side (count up)
	d_x = b_width2/dble(n_x)
	do i = 1, n_y
		y2(1, i) = -0.5d0*b_width1 
		y2(2, i) =  0.5d0*b_width2 +d_x*dble(i)
	end do 
! 	y2(:, :) = c_outer*y2(:, :)




	write(*, *) "RUN subroutine input: setting grid"
	d_grid1 = b_width1/dble(n_grid)
	d_grid2 = b_width2/dble(n_grid)

		! grid
	do i = -n_grid, 2*n_grid
		grid(1, i) = -0.5d0*b_width1 +d_grid1*dble(i)
		grid(2, i) = -0.5d0*b_width2 +d_grid2*dble(i)
	end do

end subroutine input
end module globals

