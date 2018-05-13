	program solid
	implicit none


! Declaring variables
	real, parameter :: epsilon_kcalmol = 0.997, sigma = 3.40
	real, parameter :: density = 1.616, mw_Ar = 39.948, Na = 6.022e23
	real :: sigma_sq, r, sigmar_sq, mol_per_cm3, x, y, z, xx, yy, zz
	real :: xj, yj, zj, rij_sq, rij, box_dim_cm, v, l
	real :: vol_cm3, particle_per_cm3, box_dim_Ang, vij
	real :: vLJ, r_sq, vi, red_vLJ, X_i, Y_i, Z_i, X_j, Y_j, Z_j 
	integer :: i, j, N, k, count_x
	real :: arrayi(3, 1000)

! File to store coordinates for my plot
	open(2, file = "coordinates_solid.dat")


! Computing the volume for 1000 particles

	mol_per_cm3 = density/mw_Ar
	particle_per_cm3 = mol_per_cm3*Na
	vol_cm3 = 1000/particle_per_cm3

! Coordinate in any direction (x, y, z) can be between 0 and box_dm in dm
	box_dim_cm = (vol_cm3)**(1.0/3.0)

! Converting cm to A
	box_dim_Ang = box_dim_cm*10**8

! 1000 atoms in a cube means 10 x 10 x 10 atoms on the "lines" of the cube
! Thus atom spacing is:
	l = box_dim_Ang/9.0
	
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

! Saving to my file that I will later use for a graph
	write(2, *) xx, yy, zz
	count_x = count_x + 1
	
	enddo
	enddo
	enddo


        v = 0
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


! Calculating distance
	rij_sq = ((X_i-X_j)**2) + ((Y_i-Y_j)**2) + ((Z_i-Z_j)**2)
	rij = rij_sq**(1.0/2.0)

	sigma_sq = sigma**2
	r_sq = rij**2
	sigmar_sq = sigma_sq/r_sq
	vij = (sigmar_sq**6)-(sigmar_sq**3)

	v = v + vij

	enddo

	enddo


! Calculating L-J potential

	vLJ = 4*epsilon_kcalmol*v
	red_vLJ = vLJ/epsilon_kcalmol
	

	print*, "L-J Potential for this system is: ", vLJ, "kcal per mole."
	print*, "In reduced units that is: ", red_vLJ, "epsilons."


	end program solid
