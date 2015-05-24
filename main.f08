module input
	use const
	use sde
	implicit none
	
	! input c
	!! real 8 -> 16

	real(8), parameter :: p_mass = 1.d0 ! au
! 	real(8), parameter :: p_mass = 0.51099906d6/si_e_evt/sc_c**2.d0/au_m ! eV to au
	real(8), parameter :: p_energy = 20.d0/si_e_evt/au_hr ! eV to au
	real(8), parameter :: p_angle = 60.d0/360.d0*(2.d0*mt_pi) ! do to radian (0~90, not 0)

! 	real(8), parameter :: b_poten = 0.d0/si_e_evt/au_hr ! eV to au	
	real(8), parameter :: b_poten = 30.d0/si_e_evt/au_hr ! eV to au	
	real(8), parameter :: b_width_x = 10.d0/1.d10/au_bh ! aa to au
	real(8), parameter :: b_width_y = 10.d0/1.d10/au_bh ! aa to au
	real(8), parameter :: b_width = (b_width_x**2.d0 +b_width_y**2.d0)**0.5d0 

	! calculate c

	integer, parameter :: cal_n = 2**8
! 	integer, parameter :: cal_e = 2**10
	complex(8), parameter :: psi_k_outter = mt_i*(2.d0*p_mass*p_energy)**0.5d0
	integer, parameter :: psi_n = 10**2
	integer, parameter :: psi_m = 10**1

! 	real(8), save :: f_energy
! 	real(8), save :: psi_e(0:cal_e)
! 	real(8), parameter :: psi_de = 3.d0*b_poten/cal_e

contains


	! equation d/dx y(x) = f(x, y)

function f(x, y) result(result) 
	real(8), intent(in) :: x
	complex(8), intent(in) :: y(1:2)   
	complex(8) :: result(1:2)

	result(1) = psi_k_outter*1.d0*y(2)
	result(2) = psi_k_outter*(1.d0 - b_poten/p_energy)*y(1)
! 	result(2) = psi_k_outter*(1.d0 - b_poten/f_energy)*y(1)

end function f


	! cal 1 ! input-ref-tran 

subroutine cal(width, r_zero, input, ref, tran, fun)
	real(8), intent(in) :: width, r_zero
	complex(8), intent(out) :: input, ref, tran
	complex(8), intent(out) :: fun(1:2, 0:cal_n)
	
	real(8), save :: r(0:cal_n), dr
	complex(8), save :: tem(1:2)
	integer :: i 

	dr = width/dble(cal_n)
	do i = 0, cal_n
		r(i) = 0.d0 +dr*dble(i)
	enddo

	fun(:, :) = 0.d0
	fun(:, cal_n) = 1.d0

	do i = cal_n -1, 0, -1
		tem(:) = fun(:, i +1)
		call rk8_vec(f, dr, r(i +1), tem(:))
		fun(:, i) = tem(:)
	enddo

	input = (fun(1, 0) +fun(2, 0))/2.d0 ! @
	fun(:, :) = fun(:, :)/input ! @

	input = (fun(1, 0) +fun(2, 0))/2.d0 	*exp(psi_k_outter*r_zero) 
	  ref = (fun(1, 0) -fun(2, 0))/2.d0 	*exp(psi_k_outter*r_zero) 
	 tran = fun(1, cal_n) 					*exp(psi_k_outter*r_zero)
	fun(:, :) = fun(:, :) 					*exp(psi_k_outter*r_zero)

! 	write(*, *) '       input : ', abs((fun(1, 0) +fun(2, 0))/2.d0)**2.d0
! 	write(*, *) '  reflection : ', abs((fun(1, 0) -fun(2, 0))/2.d0)**2.d0
! 	write(*, *) 'transmission : ', abs(fun(1, cal_n))**2.d0
! 	write(*, *) '       total : ', abs((fun(1, 0) -fun(2, 0))/2.d0)**2.d0 +abs(fun(1, cal_n))**2.d0
	write(*, *) abs((fun(1, 0) -fun(2, 0))/2.d0)**2.d0 +abs(fun(1, cal_n))**2.d0
	write(*, *) abs((fun(1, 0) -fun(2, 0))/2.d0)**2.d0, abs(fun(1, cal_n))**2.d0

! 	if(width == 0.d0 .or. r < 0.d0) then
! 		result = c_input*exp(psi_k_outter*psi_x(i)) +c_refle*exp(-psi_k_outter*psi_x(i))
! 	else if(r > 0.d0 .and. r < width) then
! 		result = fun(1, i)
! 		result = c_trans*exp(psi_k_outter*(psi_x(i)-width))

end subroutine cal
end module input










program main
	use const
	use input
	implicit none

	real(8), save :: x(-psi_n:psi_n, -psi_n:psi_n), dx, tmp_x
	real(8), save :: y(-psi_n:psi_n, -psi_n:psi_n), dy, tmp_y
	real(8), save :: l(-psi_n:psi_n), dl
	real(8), save :: psi_width, psi_zero 
	complex(8), save :: c_input, c_refle, c_trans
	complex(8), save :: psi_fun(1:2, 0:cal_n)
	integer :: i, j, u 




	! zahyo ! 

	dl = b_width/dble(psi_n)
	dx = dl
	dy = dl
	do j = -psi_n, psi_n
		do i = -psi_n, psi_n
			x(i, j) = dx*dble(i)*cos(p_angle) -dy*dble(j)*sin(p_angle)
			y(i, j) = dx*dble(i)*sin(p_angle) +dy*dble(j)*cos(p_angle)
		enddo
		l(j) = dl*dble(j)
	enddo
	dx = dl*cos(p_angle)
	dy = dl*sin(p_angle)

	write(*, *) 'zahyo'





	! calculate ! 
	open(03, file = 'psi_square.d')
	do j = -psi_n, psi_n

		tmp_y = tan(p_angle)*(0.d0 -x(0 ,j)) +y(0, j)
		tmp_x = 1.d0/tan(p_angle)*(0.d0 -y(0, j)) +x(0, j)
		if(tmp_y > 0.d0 .and. tmp_y < b_width_y) then 
		write(*, *) 'calculate', 'case 1'
! 		write(*, *) tmp_x, tmp_y, b_width_x, b_width_y
			psi_zero = ((x(0, j) -0.d0)**2.d0 +(y(0, j) -tmp_y)**2.d0)**0.5d0

			tmp_x = 1.d0/tan(p_angle)*(b_width_y -y(0, j)) +x(0, j)
			if(tmp_x > 0.d0 .and. tmp_x < b_width_x) then 

				psi_width = ((tmp_x)**2.d0 +(b_width_y -tmp_y)**2.d0)**0.5d0

			else 

				tmp_x = tan(p_angle)*(b_width_x -x(0 ,j)) +y(0, j)
				psi_width = ((b_width_x)**2.d0 +(tmp_x -tmp_y)**2.d0)**0.5d0

			endif
		else if(tmp_x > 0.d0 .and. tmp_x < b_width_x) then
		write(*, *) 'calculate', 'case 2'

			psi_zero = ((x(0, j) -tmp_x)**2.d0 +(y(0, j) -0.d0)**2.d0)**0.5d0

			tmp_y = tan(p_angle)*(b_width_y -x(0 ,j)) +y(0, j)
			if(tmp_y > 0.d0 .and. tmp_y < b_width_y) then 

				psi_width = ((b_width_x -tmp_x)**2.d0 +(tmp_y)**2.d0)**0.5d0

			else 

				tmp_y = 1.d0/tan(p_angle)*(b_width_y -y(0, j)) +x(0, j)
				psi_width = ((tmp_x - tmp_y)**2.d0 +(b_width_y)**2.d0)**0.5d0

			endif 
		else 
		write(*, *) 'calculate', 'case 3'
! 		write(*, *) x(0, j), y(0, j)
! 		write(*, *) tmp_x, tmp_y, b_width_x, b_width_y

			psi_width = 0.d0
			psi_zero = 0.d0

		endif





! 	! cal 1 ! input-ref-trans

! 	do i = -cal_n, 2*cal_n
! 		psi_x(i) = 0.d0 +psi_dx*dble(i)
! 	enddo

! 	f_energy = p_energy
! 	psi_fun(:, :) = 0.d0
! 	psi_fun(:, cal_n) = 1.d0

! 	do i = cal_n -1, 0, -1
! 		psi_tem(:) = psi_fun(:, i +1)
! 		call rk8_vec(f, -psi_dx, psi_x(i +1), psi_tem(:))
! 		psi_fun(:, i) = psi_tem(:)
! 	enddo

! 	c_input = (psi_fun(1, 0) +psi_fun(2, 0))/2.d0
! 	psi_fun(:, :) = psi_fun(:, :)/c_input

! 	c_input = (psi_fun(1, 0) +psi_fun(2, 0))/2.d0
! 	c_refle = (psi_fun(1, 0) -psi_fun(2, 0))/2.d0
! 	c_trans = psi_fun(1, cal_n)

! 	write(*, *) '       input : ', abs((psi_fun(1, 0) +psi_fun(2, 0))/2.d0)**2.d0
! 	write(*, *) '  reflection : ', abs((psi_fun(1, 0) -psi_fun(2, 0))/2.d0)**2.d0
! 	write(*, *) 'transmission : ', abs(psi_fun(1, cal_n))**2.d0
! 	write(*, *) '       total : ', abs((psi_fun(1, 0) -psi_fun(2, 0))/2.d0)**2.d0 +abs(psi_fun(1, cal_n))**2.d0





	! plot square poten ! 

	if(psi_width == 0.d0) then 
		write(*, *) 'plot', 'case 1'
		c_input = 1.d0 
		c_refle = 0.d0
		c_trans = 1.d0

		do i = -psi_n, psi_n
			write(03, *) x(i, j)*1.d10*au_bh, y(i, j)*1.d10*au_bh, &
							abs(c_input*exp(psi_k_outter*l(i)))**2.d0
		enddo
		write(03, *)

	else 
		write(*, *) 'plot', 'case 2'
		call cal(psi_width, psi_zero, c_input, c_refle, c_trans, psi_fun)

		do i = -psi_n, psi_n
			if(l(i) < psi_zero) then 

				write(03, *) x(i, j)*1.d10*au_bh, y(i, j)*1.d10*au_bh, & 
! 								abs(c_input*exp(psi_k_outter*l(i))) ! @
								abs(c_input*exp(psi_k_outter*l(i)) +c_refle*exp(-psi_k_outter*l(i)))**2.d0

			else if(l(i) >= psi_zero .and. l(i) < psi_zero +psi_width) then

				do u = 0, cal_n, cal_n/psi_m
					tmp_x = (psi_zero +psi_width/dble(cal_n)*dble(u))*cos(p_angle) -dl*dble(j)*sin(p_angle)
					tmp_y = (psi_zero +psi_width/dble(cal_n)*dble(u))*sin(p_angle) +dl*dble(j)*cos(p_angle)
! 					tmp_x = (psi_width/dble(psi_n)*dble(u))*cos(p_angle) -dl*dble(j)*sin(p_angle)
! 					tmp_y = (psi_width/dble(psi_n)*dble(u))*sin(p_angle) +dl*dble(j)*cos(p_angle)

					write(03, *) tmp_x*1.d10*au_bh, &
								 tmp_y*1.d10*au_bh, &
									abs(psi_fun(1, u))**2.d0
! 					write(03, *)
				enddo

			else 

				write(03, *) x(i, j)*1.d10*au_bh, y(i, j)*1.d10*au_bh, &
! 								abs(0.d0*exp(psi_k_outter*(l(i)-b_width)))**2.d0 ! @
								abs(c_trans*exp(psi_k_outter*(l(i)-b_width)))**2.d0

			endif
		enddo
		write(03, *)

	endif 
	write(*, *) 'hito roof owari', j
	enddo

	close(03)
end program main
