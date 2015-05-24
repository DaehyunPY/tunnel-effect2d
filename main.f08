module input
	use const
	implicit none
	
	! input c
	real(8), parameter :: p_mass = 1.d0 ! au
! 	real(8), parameter :: p_mass = 0.51099906d6/si_e_evt/sc_c**2.d0/au_m ! eV to au
	real(8), parameter :: p_energy = 40.d0/si_e_evt/au_hr ! eV to au
	real(8), parameter :: b_poten = 30.d0/si_e_evt/au_hr ! eV to au
	real(8), parameter :: b_poten_y = 30.d0/si_e_evt/au_hr ! eV to au
	real(8), parameter :: b_width = 10.d0/1.d10/au_bh ! aa to au

	! calculate c
	integer, parameter :: cal_n = 2**8
	integer, parameter :: cal_2d = 2**3 ! must be smaller than n 
	integer, parameter :: cal_e = 2**10

	complex(8), parameter :: psi_k_outter = mt_i*(2.d0*p_mass*p_energy)**0.5d0
! 	complex(8), parameter :: psi_k_inner  =      (2.d0*p_mass*b_poten )**0.5d0
	complex(8), save :: psi_fun(1:2, 0:cal_n)
	real(8), save :: psi_x(-cal_n:2*cal_n)
	real(8), parameter :: psi_dx = b_width/dble(cal_n)

	real(8), save :: f_energy
	real(8), save :: psi_e(0:cal_e)
	real(8), parameter :: psi_de = 3.d0*b_poten/cal_e
contains

function f(x, y) result(result) ! equation d/dx y(x) = f(x, y)
	real(8), intent(in) :: x
	complex(8), intent(in) :: y(1:2)   
	complex(8) :: result(1:2)

	result(1) = psi_k_outter*1.d0*y(2)
	result(2) = psi_k_outter*(1.d0 - b_poten/f_energy)*y(1)
end function f
end module input

program main
	use const
	use input
	use sde
	implicit none

	complex(8) :: psi_tem(1:2), c_input, c_refle, c_trans
	integer :: i, j










	! cal 1 ! input-ref-trans

	do i = -cal_n, 2*cal_n
		psi_x(i) = 0.d0 +psi_dx*dble(i)
	enddo

	f_energy = p_energy
	psi_fun(:, :) = 0.d0
	psi_fun(:, cal_n) = 1.d0

	do i = cal_n -1, 0, -1
		psi_tem(:) = psi_fun(:, i +1)
		call rk8_vec(f, -psi_dx, psi_x(i +1), psi_tem(:))
		psi_fun(:, i) = psi_tem(:)
	enddo

	c_input = (psi_fun(1, 0) +psi_fun(2, 0))/2.d0
	psi_fun(:, :) = psi_fun(:, :)/c_input

	c_input = (psi_fun(1, 0) +psi_fun(2, 0))/2.d0
	c_refle = (psi_fun(1, 0) -psi_fun(2, 0))/2.d0
	c_trans = psi_fun(1, cal_n)

	write(*, *) '       input : ', abs((psi_fun(1, 0) +psi_fun(2, 0))/2.d0)**2.d0
	write(*, *) '  reflection : ', abs((psi_fun(1, 0) -psi_fun(2, 0))/2.d0)**2.d0
	write(*, *) 'transmission : ', abs(psi_fun(1, cal_n))**2.d0
	write(*, *) '       total : ', abs((psi_fun(1, 0) -psi_fun(2, 0))/2.d0)**2.d0 +abs(psi_fun(1, cal_n))**2.d0










	! plot square poten ! 

! 	open(01, file = 'psi_simple.d')
! 	open(02, file = 'psi_euler.d')
	open(03, file = 'psi_square.d')


	do j = -cal_2d, 0
		do i = -cal_n, 2*cal_n
			write(03, *) psi_x(i)*1.d10*au_bh, psi_x(j)*1.d10*au_bh, & 
							abs(c_input*exp(psi_k_outter*psi_x(i)))**2.d0
		enddo
		write(03, *)
	enddo
	write(03, *)

	do j = 0, cal_2d
		do i = -cal_n, -1
	! 			write(01, *) psi_x(i)*1.d10*au_bh, & 
	! 							abs(c_input*exp(psi_k_outter*psi_x(i)) +c_refle*exp(-psi_k_outter*psi_x(i)))**2.d0
	! 			write(02, *) psi_x(i)*1.d10*au_bh, & 
	! 							real(c_input*exp(psi_k_outter*psi_x(i)) +c_refle*exp(-psi_k_outter*psi_x(i))), &
	! 							imag(c_input*exp(psi_k_outter*psi_x(i)) +c_refle*exp(-psi_k_outter*psi_x(i)))
			write(03, *) psi_x(i)*1.d10*au_bh, psi_x(j)*1.d10*au_bh, & 
							abs(c_input*exp(psi_k_outter*psi_x(i)) +c_refle*exp(-psi_k_outter*psi_x(i)))**2.d0
		enddo
		do i = 0, cal_n
	! 			write(01, *) psi_x(i)*1.d10*au_bh, abs(psi_fun(1, i))**2.d0
	! 			write(02, *) psi_x(i)*1.d10*au_bh, real(psi_fun(1, i)), imag(psi_fun(1, i))
			write(03, *) psi_x(i)*1.d10*au_bh, psi_x(j)*1.d10*au_bh, abs(psi_fun(1, i))**2.d0
		enddo
		do i = cal_n +1, 2*cal_n
	! 			write(01, *) psi_x(i)*1.d10*au_bh, &
	! 							abs(c_trans*exp(psi_k_outter*(psi_x(i)-b_width)))**2.d0
	! 			write(02, *) psi_x(i)*1.d10*au_bh, &
	! 							real(c_trans*exp(psi_k_outter*(psi_x(i)-b_width))), &
	! 							imag(c_trans*exp(psi_k_outter*(psi_x(i)-b_width)))
			write(03, *) psi_x(i)*1.d10*au_bh, psi_x(j)*1.d10*au_bh, &
							abs(c_trans*exp(psi_k_outter*(psi_x(i)-b_width)))**2.d0
		enddo
		write(03, *)
	enddo
	write(03, *)

	do j = cal_2d +1, 2*cal_2d
		do i = -cal_n, 2*cal_n
			write(03, *) psi_x(i)*1.d10*au_bh, psi_x(j)*1.d10*au_bh, & 
							abs(c_input*exp(psi_k_outter*psi_x(i)))**2.d0
		enddo
		write(03, *)
	enddo
	write(03, *)

! 	close(01)
! 	close(02)
	close(03)










! 	! cal 2 ! plot energy-ref/trans

! 	do i = 0, cal_e
! 		psi_e(i) = 0.d0 +psi_de*dble(i)
! 	enddo

! 	open(11, file = 'trans.d')
! 	open(12, file = 'refle.d')
! 	do i = 0, cal_e
! 		f_energy = psi_e(i)
! 		psi_fun(:, :) = 0.d0
! 		psi_fun(:, cal_n) = 1.d-10

! 		do j = cal_n -1, 0, -1
! 			psi_tem(:) = psi_fun(:, j +1)
! 			call rk8_vec(f, -psi_dx, psi_x(j +1), psi_tem(:))
! 			psi_fun(:, j) = psi_tem(:)
! 		enddo

! 		c_input = (psi_fun(1, 0) +psi_fun(2, 0))/2.d0
! 		c_refle = (psi_fun(1, 0) -psi_fun(2, 0))/2.d0
! 		c_trans = psi_fun(1, cal_n)
! ! 		write(11, *) psi_e(i)*si_e_evt*au_hr, abs(c_trans/c_input)**2.d0
! ! 		write(12, *) psi_e(i)*si_e_evt*au_hr, abs(c_refle/c_input)**2.d0
! 		write(11, *) psi_e(i)/b_poten, abs(c_trans/c_input)**2.d0
! 		write(12, *) psi_e(i)/b_poten, abs(c_refle/c_input)**2.d0
! 	enddo
! 	close(11)
! 	close(12)

end program main
