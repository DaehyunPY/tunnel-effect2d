module globals
	use const
	use subfuncs
	implicit none

	integer, parameter :: n_l = 2**8
! 	integer, parameter :: n_eq = 3 * 2 * (2 * n_l +1)
	integer, parameter :: n_eq = 3 * 4 * (2 * n_l)
	integer, parameter :: n_r = n_eq / 8, n_theta = n_eq / 8
! 	integer, parameter :: n_inner = 2**8
	integer, parameter :: n_inner = n_l / 2**5
! 	integer, parameter :: n_inner = n_l
	real(8), parameter :: rg_outer = dble(n_r) / dble(n_inner)

! 	integer, parameter :: n_l = 2**9
! 	integer, parameter :: n_eq = 2 * 2 * (2 * n_l +1)
! 	integer, parameter :: n_r = n_eq / 3, n_theta = n_eq / 3
! 	integer, parameter :: n_inner = 2**8
! 	real(8), parameter :: rg_outer = dble(n_r) / dble(n_inner)

! 	integer, parameter :: n_r = 2**10, n_theta = 2**10
! 	real(8), parameter :: rg_outer = 2.d0
! 	integer, parameter :: n_inner = n_r, n_outer = rg_outer * n_r
! 	integer, parameter :: n_eq = 2 * (n_theta +2 * n_outer)
	
	complex(8) k_inner, k_outer, k_input(1 : 2)
	real(8)    r(0 : n_r), theta(0 : n_theta)
! 	real(8)    r(0 : n_inner +n_outer), theta(0 : n_theta)

	complex(8) eq_A(1 : n_eq, 1 : n_eq), eq_x(1 : n_eq), eq_b(1 : n_eq) ! Ax = b
! 	real(8)    eq_rg1_r, eq_rg1_theta(1 : n_theta)
! 	real(8)    eq_rg2_r(1 : n_r), eq_rg2_theta
! 	real(8)    eq_rg3_r(1 : n_r), eq_rg3_theta

! 	complex(8) c_inner(-n_l : n_l, 1 : 2), c_input, c_ref(-n_l : n_l, 1 : 2), c_trans(-n_l : n_l, 1 : 2)
! 	complex(8) c_inner(-n_l : n_l, 1 : 2), c_input, c_ref(-n_l : n_l), c_trans(-n_l : n_l)

	integer, parameter :: n_grid = 100
	real(8) grid(1 : 4, 0 : n_grid, 0 : n_grid)

	real(8), parameter :: au_aa = 1.d10*au_bh
	
contains





subroutine input 
	real(8) p_mass,  p_energy, p_angle
	real(8) b_poten, b_width1, b_width2, b_radius

	integer, parameter :: file_input = 100

	character dum_ch
	complex(8) inner, outer
	real(8) tmp, d_r, d_theta, d_grid
	integer i, j





! ===== reading =============================================

	write(*, *) "RUN subroutine input: reading"
	write(*, *) n_eq, 2 * 2 * (n_r +n_theta)

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
	read(file_input, *) dum_ch, b_radius	! 'WIDTH[x] 			[aa]'					10.d0
	read(file_input, *)
! 	read(file_input, *)
! 	read(file_input, *) dum_ch				! '=================================================='
! 	read(file_input, *) dum_ch				! 'CALCULATE'
! 	read(file_input, *) dum_ch				! '=================================================='
! 	read(file_input, *) dum_ch, n_total		! 'CALCULATE NUMBER (POWER OF TWO)'			2**8

	close(file_input)	
	p_angle = 0.d0 ! POTENTIAL TYPE II set p_angle zero



! 	if(n_total > n_max) stop "ERROR subroutine input: n_total is too big"

! 	p_mass   = p_mass   /si_e_evt/sc_c**2.d0/au_m 	! eV to au
	p_energy = p_energy /si_e_evt/au_hr 			! eV to au
	p_angle  = p_angle  /360.d0*(2.d0*mt_pi) 		! do to radian
	b_poten  = b_poten  /si_e_evt/au_hr 			! eV to au	
	b_width1 = b_width1 /1.d10/au_bh 				! aa to au
	b_width2 = b_width2 /1.d10/au_bh 				! aa to au
	b_radius = b_radius /1.d10/au_bh 				! aa to au

! 	n_k = n_total/2
! 	n_x = n_total/8





! ===== k, r, theta =============================================

	write(*, *) "RUN subroutine input: setting k, r, theta"

		! k
	tmp = -2.d0*p_mass*p_energy
	k_outer = mt_i*(-tmp)**0.5d0

	tmp = -2.d0*p_mass*(p_energy -b_poten)
	if(tmp < 0.d0) then 
		k_inner = mt_i*(-tmp)**0.5d0
	else if(tmp > 0.d0) then 
		k_inner = -tmp**0.5d0
	else
		stop "p_energy and b_poten is same, it is not supported"
	end if

	k_input(1)  = k_outer*cos(p_angle)
	k_input(2)  = k_outer*sin(p_angle)



		! r, theta
	r = 0.d0
	d_r = b_radius / dble(n_inner)
	do i = 0, n_r
		r(i) = b_radius +d_r * dble(i)
	end do 
	write(*, *) "outer region: ", rg_outer
! 	d_r = b_radius / dble(n_inner)
! 	do i = 0, n_inner +n_outer
! 		r(i) = d_r * dble(i)
! 	end do 

	theta = 0.d0
	d_theta = mt_pi / dble(n_theta)
	do i = 0, n_theta
		theta(i) = d_theta * dble(i)
	end do 





! ===== eq, grid =============================================

	write(*, *) "RUN subroutine input: setting eq, grid"

		! eq 
	eq_A = 0.d0
	eq_x = 0.d0
	eq_b = 0.d0

! 	eq_rg1_r = b_radius
! 	eq_rg1_theta = 0.d0
! 	eq_rg1_theta(1:) = theta(1:)

! 	eq_rg2_r = 0.d0
! 	eq_rg2_r(1:) = r(n_inner +1:)
! 	eq_rg2_theta = 0.5d0 * mt_pi

! 	eq_rg3_r = 0.d0
! 	eq_rg3_r(1:) = r(n_inner +1:)
! 	eq_rg3_theta = 1.5d0 * mt_pi



		! grid
	grid = 0.d0
	d_grid = 6.d0 * b_radius / dble(n_grid)

	do i = 0, n_grid
		do j = 0, n_grid
			grid(1, i, j) = -3.d0 * b_radius +d_grid * dble(i)
			grid(2, i, j) = -3.d0 * b_radius +d_grid * dble(j)
			grid(3, i, j) = to_r    (grid(1, i, j), grid(2, i, j))
			grid(4, i, j) = to_theta(grid(1, i, j), grid(2, i, j))
		end do
	end do

end subroutine input
end module globals




