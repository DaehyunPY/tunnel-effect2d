module const
	implicit none 

!basic
	real(8),    parameter :: mt_pi = 2.d0*acos(0.d0)
	real(8),    parameter :: mt_e  = exp(1.d0)
	complex(16), parameter :: mt_i  = cmplx(0.d0, 1.d0, kind(0.d0))

	real(8), parameter :: sc_c  = 2.99792458d8					! (m s-)si
	real(8), parameter :: sc_na = 6.02214129d+23				! (mol-)si
	real(8), parameter :: sc_kb = 1.38065042d-23				! (J K-)si
	real(8), parameter :: sc_ep = 1.d+7/(4*mt_pi*sc_c **2.d0)	! (F m-)si
	real(8), parameter :: sc_mu = 1.d-7*(4*mt_pi)				!(N A-2)si
	real(8), parameter :: sc_m  = 9.10938262d-31				!   (kg)si
	real(8), parameter :: sc_q  = 1.60217653d-19				!    (C)si
	real(8), parameter :: sc_hb = 1.05457173d-34				!  (J s)si

!si
	real(8), parameter :: si_debye = sc_c /1.d-21						!(debye)other si-(C m)
	real(8), parameter :: si_w_pcm = 1.d-2/(2.d0*mt_pi*sc_c )			!  (cm-)other si-(Hz)
	real(8), parameter :: si_e_pcm = 1.d-2/(2.d0*mt_pi*sc_c )/sc_hb	!  (cm-)other si-(J)
	real(8), parameter :: si_e_evt = 1.d0/sc_q							!   (eV)other si-(J)
	
!au
	real(8), parameter :: au_hr = sc_m *sc_q **4.d0/(4.d0*mt_pi*sc_ep*sc_hb)**2.d0	!(J)si au-(hartree)
	real(8), parameter :: au_bh = 4.d0*mt_pi*sc_ep*sc_hb**2.d0/(sc_m *sc_q **2.d0)	!(m)si au-(bohr)
	real(8), parameter :: au_m  = sc_m					!  (kg)si au-
	real(8), parameter :: au_t  = sc_hb/au_hr			!   (s)si au-
	real(8), parameter :: au_tm = au_hr/sc_kb			!   (K)si au-
	real(8), parameter :: au_ec = sc_q					!   (C)si au-
	real(8), parameter :: au_ef = au_hr/au_bh/au_ec	!(V m-)si au-
	real(8), parameter :: au_mq = au_hr/au_ec*au_t		!  (Wb)si au-
	real(8), parameter :: au_mf = au_hr/au_bh/au_mq	!(A m-)si au-

end module const