	program liquid
	implicit none

! Declaring variables
	real, parameter :: density = 1.395, mw_Ar = 39.948, Na = 6.022e23
	real :: mol_per_cm3, xi, yi, zi, box_dim_cm
	real :: vol_cm3, particle_per_cm3, box_dim_Ang
	integer :: i, j, N
	real :: arrayi(3, 1000)

! File to store coordinates for my plot
	open(2, file = "coordinates_liquid.dat")


! Computing the volume for 1000 particles

	mol_per_cm3 = density/mw_Ar
	particle_per_cm3 = mol_per_cm3*Na
	vol_cm3 = 1000/particle_per_cm3

! Coordinate in any direction (x, y, z) can be between 0 and box_dm in dm
	box_dim_cm = (vol_cm3)**(1.0/3.0)

! Converting dm to A
	box_dim_Ang = box_dim_cm*10**8

	N = 1000

	do i = 1, N
	
! Generating coordinates of the atom
	xi = rand()*box_dim_Ang
	yi= rand()*box_dim_Ang
	zi = rand()*box_dim_Ang


! Saving my coordinates in an array
	arrayi(1,i) = xi
	arrayi(2,i) = yi
	arrayi(3,i) = zi

! Saving to my file that I will later use for a graph
	write(2, *) xi, yi, zi

	enddo


	print*, "File coordinates_liquid.dat generated."

	end program liquid