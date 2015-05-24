module const
	implicit none 

!basic
	double precision,    parameter :: mt_pi = 2.d0*acos(0.d0)
	double precision,    parameter :: mt_e  = exp(1.d0)
	complex(kind(0.d0)), parameter :: mt_i  = cmplx(0.d0, 1.d0, kind(0.d0))

	double precision, parameter :: sc_c  = 2.99792458d8					! (m s-)si
	double precision, parameter :: sc_na = 6.02214129d+23				! (mol-)si
	double precision, parameter :: sc_kb = 1.38065042d-23				! (J K-)si
	double precision, parameter :: sc_ep = 1.d+7/(4*mt_pi*sc_c **2.d0)	! (F m-)si
	double precision, parameter :: sc_mu = 1.d-7*(4*mt_pi)				!(N A-2)si
	double precision, parameter :: sc_m  = 9.10938262d-31				!   (kg)si
	double precision, parameter :: sc_q  = 1.60217653d-19				!    (C)si
	double precision, parameter :: sc_hb = 1.05457173d-34				!  (J s)si

!si
	double precision, parameter :: si_debye = sc_c /1.d-21						!(debye)other si-(C m)
	double precision, parameter :: si_w_pcm = 1.d-2/(2.d0*mt_pi*sc_c )			!  (cm-)other si-(Hz)
	double precision, parameter :: si_e_pcm = 1.d-2/(2.d0*mt_pi*sc_c )/sc_hb	!  (cm-)other si-(J)
	double precision, parameter :: si_e_evt = 1.d0/sc_q							!   (eV)other si-(J)
	
!au
	double precision, parameter :: au_hr = sc_m *sc_q **4.d0/(4.d0*mt_pi*sc_ep*sc_hb)**2.d0	!(J)si au-(hartree)
	double precision, parameter :: au_bh = 4.d0*mt_pi*sc_ep*sc_hb**2.d0/(sc_m *sc_q **2.d0)	!(m)si au-(bohr)
	double precision, parameter :: au_m  = sc_m					!  (kg)si au-
	double precision, parameter :: au_t  = sc_hb/au_hr			!   (s)si au-
	double precision, parameter :: au_tm = au_hr/sc_kb			!   (K)si au-
	double precision, parameter :: au_ec = sc_q					!   (C)si au-
	double precision, parameter :: au_ef = au_hr/au_bh/au_ec	!(V m-)si au-
	double precision, parameter :: au_mq = au_hr/au_ec*au_t		!  (Wb)si au-
	double precision, parameter :: au_mf = au_hr/au_bh/au_mq	!(A m-)si au-

end module const