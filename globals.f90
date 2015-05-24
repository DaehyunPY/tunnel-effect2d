module globals
	use const
	implicit none

	integer, parameter :: n_cal = 10**5
	complex(8) refle, trans
	complex(8) k_outer, inner(1:2, 0:n_cal)
	   real(8) p_energy, b_charge 
	   real(8) r, y, x1, x2, x(0:n_cal)
	   integer info

	integer, parameter :: n_grid = 101
	real(8) grid(1:2, 0:n_grid)

	real(8), parameter :: au_aa = 1.d10*au_bh
	
contains





subroutine input 
	real(8) p_mass, p_angle, b_poten, b_radius
	integer, parameter :: file_input = 100

	character dum_ch
	real(8) tmp, d_grid
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
	read(file_input, *) dum_ch, b_charge	! 'CHARGE			 	[e]'					1.d0
	read(file_input, *) dum_ch, b_radius	! 'WIDTH[x] 			[aa]'					10.d0
	read(file_input, *)

	close(file_input)	



! 	p_mass   = p_mass   /si_e_evt/sc_c**2.d0/au_m 	! eV to au
	p_energy = p_energy /si_e_evt/au_hr 			! eV to au
	p_angle  = p_angle  /360.d0*(2.d0*mt_pi) 		! do to radian
	b_poten  = b_poten  /si_e_evt/au_hr 			! eV to au
	b_charge = b_charge	
	b_radius = b_radius /1.d10/au_bh 				! aa to au

! 	write(*, *) n_total



	write(*, *) "RUN subroutine input: setting k"
	r = b_radius 
	tmp = -2.d0*p_mass*p_energy
	k_outer = mt_i*(-tmp)**0.5

	write(*, *) "RUN subroutine input: setting grid"
	d_grid = 3.d0*b_radius/dble(n_grid)
	do i = 0, n_grid
		grid(1, i) = -1.5d0*b_radius +d_grid*dble(i)
		grid(2, i) = -1.5d0*b_radius +d_grid*dble(i)
	end do

end subroutine input
end module globals

