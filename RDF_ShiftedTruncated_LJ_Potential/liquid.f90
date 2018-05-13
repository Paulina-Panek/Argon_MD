	program liquid
	implicit none
	
! Paulina Panek 03/30/2018
! This program generates data to plot radial distribution function for 1000 Ar liquid atoms
! Only one frame (nstep = 1)
! Bin size 0.05 sigma

! Declaring variables
	real, parameter :: density = 1.395, mw_Ar = 39.948, Na = 6.022e23, sigma = 3.4
	real :: mol_per_cm3, xi, yi, zi, box_dim_cm, dr
	real :: vol_cm3, particle_per_cm3, box_dim_Ang, rij_sq,  r_hi, r_lo, h_id, const, rho
	integer :: i, j, N, k, nk
	real :: arrayi(3, 1000), r(3,1000)
	integer, DIMENSION(:), allocatable :: h
	real, DIMENSION(:), allocatable :: g
	real, DIMENSION(3) :: rij
	real, parameter :: epsilon = 0.997, pi = 4.0*ATAN(1.0)

! File to store coordinates for my plot
	open(2, file = "gr_liq.dat")

! Computing the volume for 1000 particles
	mol_per_cm3 = density/mw_Ar
	particle_per_cm3 = mol_per_cm3*Na
	vol_cm3 = 1000/particle_per_cm3

! Coordinate in any direction (x, y, z) can be between 0 and box_dm in dm
	box_dim_cm = (vol_cm3)**(1.0/3.0)

! Converting dm to A
	box_dim_Ang = box_dim_cm*10**8
	
! Converting to box = 1 units
	dr = 0.05*sigma/box_dim_Ang
	
	nk = FLOOR(0.5/dr)

	N = 1000
	
	allocate (h(nk), g(nk))
	
	h(:) = 0

	do i = 1, N
	
! Generating coordinates of the atom

! Leaving the coordinates as box length = 1 unit for the purpose of this hw
	xi = rand() !*box_dim_Ang
	yi= rand() !*box_dim_Ang
	zi = rand() !*box_dim_Ang

! Saving my coordinates in an array
!	arrayi(1,i) = xi
!	arrayi(2,i) = yi
!	arrayi(3,i) = zi

! Coordinates in converted to box = 1 units
	r(1,i) = xi
	r(2,i) = yi
	r(3,i) = zi
	
	enddo
	 

	do i = 1, N-1
	do j = i + 1, N
	
	rij(:) = r(:,i) - r(:,j)
	rij(:) = rij(:) - ANINT(rij(:))
	rij_sq = SUM(rij**2)
	k = floor(SQRT(rij_sq)/dr) + 1
	
	if (k .le. nk) then
	h(k)  = h(k) + 2

	endif
	
	enddo
	enddo

! Normalization
	rho = REAL(N)
	const = 4.0 * pi * rho / 3.0
	
		do k =1, nk
		g(k) = REAL (h(k))/real(N*1)
		r_lo = REAL (k-1) *dr
		r_hi = r_lo + dr
		h_id = const * (r_hi**3 - r_lo **3)
		g(k) = g(k)/h_id
		enddo
		
		! To get x-axis in terms of sigma
		dr = dr*box_dim_Ang
		
		do k = 1, nk
			write(2,*) (REAL(k)-0.5)*dr, g(k)
		enddo

	print*, "File generated."

	end program liquid