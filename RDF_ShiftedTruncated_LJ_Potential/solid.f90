	program solid
	implicit none
	
! Paulina Panek 

! Declaring variables
	real, parameter :: epsilon_kcalmol = 0.997, sigma = 3.40
	real, parameter :: density = 1.616, mw_Ar = 39.948, Na = 6.022e23
	real :: sigma_sq, sigmar_sq, mol_per_cm3, x, y, z, xx, yy, zz, vcut2TS
	real :: xj, yj, zj, rij_sq, rij, box_dim_cm, v, l, ll, v3_5TS, vcut3_5TS, v2TS
	real :: vol_cm3, particle_per_cm3, box_dim_Ang, vij, penalty3_5TS, penalty2TS
	real :: vLJ, r_sq, vi, X_i, Y_i, Z_i, X_j, Y_j, Z_j, vLJ_rcut, sigma_over_rcut, sigma_over_rcut_sq
	real :: vLJ_rcut2, sigma_over_rcut2, sigma_over_rcut_sq2
	integer :: i, j, N, k, count_x, nk
	real :: arrayi(3, 1000), arrayA(3), r(3, 1000),  rij_array(3)
	real :: rcut2, rcut3_5, vcut2, vcut3_5, penalty2, penalty3_5, v2, v3_5, const
	integer, DIMENSION(:), allocatable :: h
	real, DIMENSION(:), allocatable :: g
	real, parameter :: pi = 4.0*ATAN(1.0)
	real :: dr, h_id, r_hi, r_lo, rho

	rcut2 = 2.0*sigma
	rcut3_5 = 3.5*sigma
	
	N = 1000

! Format for % penalty  & L-J
10	format(A, F5.2, A)
11	format(A, F9.2, A)

! File to store coordinates for my plot
	open(2, file = "gr_sol.dat")

! Computing the volume for 1000 particles

	mol_per_cm3 = density/mw_Ar
	particle_per_cm3 = mol_per_cm3*Na
	vol_cm3 = 1000/particle_per_cm3

! Coordinate in any direction (x, y, z) can be between 0 and box_dm in dm
	box_dim_cm = (vol_cm3)**(1.0/3.0)

! Converting cm to A
	box_dim_Ang = box_dim_cm*10**8

! 1000 atoms in a cube means 10 x 10 x 10 atoms on the "lines" of the cube
! Thus atom spacing in Ang is:
	l = box_dim_Ang/9.0
	
! And atom spacing in box lenght = 1 unit units is:
	ll = l/box_dim_Ang

! Generating coordinates

	count_x = 1

	do i = 1, 10 
	xx = l*(i-1)	

	do j = 1, 10
	yy = l*(j-1)
	 
	do k = 1, 10
	zz = l*(k-1)

	arrayi(1, count_x) = xx
	arrayi(2, count_x) = yy
	arrayi(3, count_x) = zz

	count_x = count_x + 1
	
	enddo
	enddo
	enddo
	
! Changing my generated coordinates for box lenght = 1 unit for g(r) part of the assignment
	r = arrayi/box_dim_Ang
	
	! Converting to box = 1 units
	dr = 0.05*sigma/box_dim_Ang
	
	nk = FLOOR(0.5/dr)

	
	allocate (h(nk), g(nk))
	
	h(:) = 0

! Generating data for g(r) plot
	do i = 1, N-1
	do j = i + 1, N

! note all calculations in box lenght = 1 unit
	rij_array(:) = r(:,i) - r(:,j)
	rij_array(:) = rij_array(:) - (1+ll)*ANINT(rij_array(:))
	rij_sq = SUM(rij_array**2)
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


! Calculating vLJ_rcut at 2 sigma
	sigma_over_rcut2 = 1.0/2.0
	sigma_over_rcut_sq2 = sigma_over_rcut2**2
	vLJ_rcut2 = sigma_over_rcut_sq2**6 - sigma_over_rcut_sq2**3

! Calculating vLJ_rcut at 3.5 sigma
	sigma_over_rcut = 1.0/3.5
	sigma_over_rcut_sq = sigma_over_rcut**2
	vLJ_rcut = sigma_over_rcut_sq**6 - sigma_over_rcut_sq**3


! Setting values to begin my loop
    v = 0
	v2 = 0
	v3_5 = 0
	v3_5TS = 0

	N = 1000

	do i = 1, N-1
! Picking first atom from the array
	X_i = arrayi(1,i)
	Y_i = arrayi(2,i)
	Z_i = arrayi(3,i)

	do j = i + 1, N
! Second atom for the pairwise interaction
	X_j = arrayi(1,j)
	Y_j = arrayi(2,j)
	Z_j = arrayi(3,j) 

! Calculating distances in each direction

	arrayA(1) = X_j - X_i
	arrayA(2) = Y_j - Y_i
	arrayA(3) = Z_j - Z_i

! l is my distance between two atoms
! making sure boundary atoms not on top of each other
	arrayA(:) = arrayA(:) - (box_dim_Ang+l)*ANINT(arrayA(:)/box_dim_Ang)
	

! Calculating distance
	rij_sq = SUM(arrayA(:)**2)
	rij = rij_sq**(1.0/2.0)

	sigma_sq = sigma**2
	r_sq = rij**2
	sigmar_sq = sigma_sq/r_sq
	vij = (sigmar_sq**6)-(sigmar_sq**3)

! potential without truncation
	v = v + vij

! truncation rcut = 2sigma
	if (rij .le. rcut2) then
	v2 = v2 + vij
	v2TS = v2TS + vij - vLJ_rcut2
	endif

! truncation rcut = 3.5sigma
	if (rij .le. rcut3_5) then
	v3_5 = v3_5 + vij
	v3_5TS = v3_5TS + vij - vLJ_rcut
	endif
	

	enddo
	enddo


! Calculating L-J potential without truncation

	vLJ = 4*epsilon_kcalmol*v
	vcut2 = 4*epsilon_kcalmol*v2
	vcut2TS = 4*epsilon_kcalmol*v2TS
	vcut3_5 = 4*epsilon_kcalmol*v3_5
	vcut3_5TS= 4*epsilon_kcalmol*v3_5TS

! Calculating penalties
	penalty2 = (vLJ - vcut2)*100/vLJ
	penalty2TS = (vLJ - vcut2TS)*100/vLJ
	penalty3_5 = (vLJ - vcut3_5)*100/vLJ
	penalty3_5TS = (vLJ - vcut3_5TS)*100/vLJ
	
	

	write(6, 11) "L-J Potential without truncation is: ", vLJ, " kcal per mole."
	print*, ""
	write(6, 11) "If truncation is 2 sigma, L-J potential is: ", vcut2, " kcal per mole"
	write(6, 10) "And the L-J energy penalty is: ", penalty2, "%."
	print*, ""
	write(6, 11) "Trancated (2 sigma) and shifted L-J potential is: ", vcut2TS, " kcal per mole"
	write(6, 10) "And the L-J energy penalty is: ", penalty2TS, "%."
	print*, ""
	write(6, 11) "If truncation is 3.5 sigma, L-J potential is: ", vcut3_5, " kcal per mole"
	write(6, 10) "And the L-J energy penalty is: ", penalty3_5, "%."
	print*, ""
	write(6, 11) "Trancated (3.5 sigma) and shifted L-J potential is: ", vcut3_5TS, " kcal per mole"
	write(6, 10) "And the L-J energy penalty is: ", penalty3_5TS, "%."


	end program solid