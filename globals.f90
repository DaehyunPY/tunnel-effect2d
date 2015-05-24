module globals
	use const
	implicit none

	integer, parameter :: n_total = 4
	complex(8) k_inner, k_outer
	complex(8) mat_xk(1:n_total, 1:n_total), vec_coeff(1:n_total), vec_input(1:n_total) ! Ax = b
	real(8)    r, y, x1, x2
	integer    info

	integer, parameter :: n_grid = 100
	real(8)    grid(1:2, 0:n_grid)

	real(8), parameter :: au_aa = 1.d10*au_bh
	
contains





subroutine input 
	real(8) p_mass,  p_energy, p_angle
	real(8) b_poten, b_radius

	integer, parameter :: file_input = 100

	character dum_ch
	real(8) tmp, d_theta, d_x, d_grid
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
	read(file_input, *) dum_ch, b_radius	! 'WIDTH[x] 			[aa]'					10.d0
	read(file_input, *)

	close(file_input)	



! 	p_mass   = p_mass   /si_e_evt/sc_c**2.d0/au_m 	! eV to au
	p_energy = p_energy /si_e_evt/au_hr 			! eV to au
	p_angle  = p_angle  /360.d0*(2.d0*mt_pi) 		! do to radian
	b_poten  = b_poten  /si_e_evt/au_hr 			! eV to au	
	b_radius = b_radius /1.d10/au_bh 				! aa to au

! 	write(*, *) n_total



	write(*, *) "RUN subroutine input: setting k"

		! r 
	r = b_radius 

		! k_outer k
! 	k_outer = -2.d0*p_mass*p_energy
	tmp = -2.d0*p_mass*p_energy
	k_outer = mt_i*(-tmp)**0.5

		! k_inner k 
	tmp = -2.d0*p_mass*(p_energy -b_poten)
	if(tmp < 0.d0) then 
		k_inner = mt_i*(-tmp)**0.5
	else if(tmp > 0.d0) then 
		k_inner = -tmp**0.5
	else
! 		stop "p_energy and b_poten is same, it is not supported"
		k_inner = 0.d0
	end if



	write(*, *) "RUN subroutine input: setting x"



	write(*, *) "RUN subroutine input: setting grid"
	d_grid = 3.d0*b_radius/dble(n_grid)

		! grid
	do i = 0, n_grid
		grid(1, i) = -1.5d0*b_radius +d_grid*dble(i)
		grid(2, i) = -1.5d0*b_radius +d_grid*dble(i)
	end do

end subroutine input
end module globals

