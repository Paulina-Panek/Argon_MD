 	program ar_md
	implicit none
	
	! Paulina Panek, April 2018
	! MD simulation of 1000 liquid Ar atoms in 85.0K
	

! Declaring variables
	real, parameter :: epsilon_kcalmol = 0.997,  pi = 4.0*ATAN(1.0), temperature = 85.0, total_t = 50*1E-12
	real, parameter :: density = 1.395, mw_Ar = 39.948, Na = 6.022e23, kb = 1.3807e-23, delta_t = 5*1E-15
	integer, parameter :: N = 1000
	real :: mol_per_cm3, x, y, z, const, sigma_sq, rij_sq, box_dim_m
	real :: box_dim_cm, l, vol_cm3, particle_per_cm3
	real :: sigma, sig_over_r_sq, num_tra, kinetic_energy
	real :: dr, h_id, r_lo, rho, r_sq, epsilon_Joules_molec, force_constant, box
	real :: x1, x2, x3, x4, x5, x6, y1, y2, y3, y4, y5, y6
	real :: vx_cm, vy_cm, vz_cm, kinE, m, v_size, v_sq, v_sq_x, v_sq_y, v_sq_z
	real :: half_step, total_energy, vLJ, J_to_kcalmol, bin_size, r_hi, rdf_range
	real :: r_plot, ro_n, box_gr
	real :: r(3, N),  rij_array(3), vel_before_norm(3, N)
	real :: fij(3), acceleration(3, N), vel_norm(3, N), v_sizes(N)
	integer :: i, j, k, p, nk, counter, counter2, num_bins, num_bins_xyz, b, bins_rdf, rdf_count
	integer :: histo_v(2, 40), histo_vxyz(4, 101)
	real, DIMENSION(:), allocatable :: g, h
	

	
	sigma = 3.40         ! Angstroms
	sigma = sigma*1E-10  ! Meters
	m = mw_Ar/(Na*1000)  ! kg (single atom)

! Computing the volume for 1000 particles
	mol_per_cm3 = density/mw_Ar
	particle_per_cm3 = mol_per_cm3*Na
	vol_cm3 = real(N)/particle_per_cm3

	box_dim_cm = (vol_cm3)**(1.0/3.0)    ! box length (cm)
	box_dim_m = box_dim_cm*0.01          ! box length (m)
	l = box_dim_m/9.0                    ! atom spacing (m)
	

! Generating coordinates of initial lattice
	counter = 1

	do i = 1, 10 
	x = l*(i-1)	
		do j = 1, 10
		y = l*(j-1)
	 		do k = 1, 10
			z = l*(k-1)

			r(1, counter) = x
			r(2, counter) = y
			r(3, counter) = z

			counter = counter + 1
			enddo
		enddo
	enddo
	
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Generating initial velocities!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	! Generating array with normal distribution vectors
	counter2 = 500

	do i =1, 500
	! Generating 1000 vectors, 2 generated in one iteration => i between 1 and 500
		x1 = rand()
		x2 = rand()
		x3 = rand()
		x4 = rand()
		x5 = rand()
		x6 = rand()
		
		y1 = sqrt(-2.0*log(x1))*cos(2.0*pi*x2)
		y3 = sqrt(-2.0*log(x1))*sin(2.0*pi*x2)
		y5 = sqrt(-2.0*log(x3))*cos(2.0*pi*x4)
		y2 = sqrt(-2.0*log(x3))*sin(2.0*pi*x4)
		y4 = sqrt(-2.0*log(x5))*cos(2.0*pi*x6)
		y6 = sqrt(-2.0*log(x5))*sin(2.0*pi*x6)
		
	! First vector
		vel_before_norm(1, i) = y1
		vel_before_norm(2, i) = y2
		vel_before_norm(3, i) = y3
		
		counter2 = counter2 + 1
		
	! Second vector	
		vel_before_norm(1, counter2) = y4
		vel_before_norm(2, counter2) = y5
		vel_before_norm(3, counter2) = y6
	enddo
	
	! Removing center of mass
	vx_cm = sum(vel_before_norm(1,:))/real(N)
	vy_cm = sum(vel_before_norm(2,:))/real(N)
	vz_cm = sum(vel_before_norm(3,:))/real(N)
	
	vel_before_norm(1,:) = vel_before_norm(1,:) - vx_cm
	vel_before_norm(2,:) = vel_before_norm(2,:) - vy_cm
	vel_before_norm(3,:) = vel_before_norm(3,:) - vz_cm
	
	! Scale to avg energy at 85.0K
	kinE = (kb*temperature/m)**(1.0/2.0)
	
! Initial velocities by component
	vel_norm(1,:) = vel_before_norm(1,:)*kinE
	vel_norm(2,:) = vel_before_norm(2,:)*kinE
	vel_norm(3,:) = vel_before_norm(3,:)*kinE
	
! Initial velocities vector size
	do i = 1, 1000
		v_sq_x = vel_norm(1,i) * vel_norm(1,i)
		v_sq_y = vel_norm(2,i) * vel_norm(2,i)
		v_sq_z = vel_norm(3,i) * vel_norm(3,i)
		
		v_sq = v_sq_x + v_sq_y + v_sq_z
		v_size = sqrt(v_sq)
		
		v_sizes(i) = v_size
	
	enddo
	
	
!!!!!!!!!!!!!!!!!!!
! Histograms code!!
!!!!!!!!!!!!!!!!!!!

! Histogram of velocity vectors
	open(2, file = "histo_v.dat")

	bin_size = 40.0
	num_bins = 40	
	
	do k = 1, num_bins
		histo_v(1, k) = int(bin_size*(k-1)) ! labels in first column
		histo_v(2, k) = 0                   ! default values set to 0 in second column
	enddo
	
	do i = 1, 1000
	
		jloop: do j = 1, num_bins
			if (v_sizes(i) .le. real(histo_v(1, j))) then
			histo_v(2, j) = histo_v(2, j) + 1
			exit jloop
			
			endif
			enddo jloop
	enddo
			
		do i = 1, num_bins
			write(2,*) histo_v(1,i), histo_v(2,i)
		enddo

! Histogram of velocity components
	open(22, file = "histo_v_xyz.dat")

	bin_size = 20.0
	num_bins_xyz = 101
	
	do i = 2, 4				! Setting starting values at 0 for the 3 columns with tally
		do k = 1, num_bins_xyz
			histo_vxyz(i, k) = 0.0
		enddo
	enddo

	
	do k = 1, num_bins_xyz
		histo_vxyz(1, k) = int((-1000.0 + bin_size*real(k-1))) ! checking for velocities from -1000 to 1000
	enddo
	

	do k = 1, 3
		do i = 1, 1000
			ploop: do p = 1, num_bins_xyz
				if (vel_norm(k,i) .le. real(histo_vxyz(1, p))) then
					histo_vxyz(k+1, p) = histo_vxyz(k+1, p) + 1
					exit ploop
				endif
				enddo ploop
		enddo
	enddo
			
	do i = 1, num_bins_xyz
		write(22,*) histo_vxyz(1,i), histo_vxyz(2,i), histo_vxyz(3,i), histo_vxyz(4,i)
	enddo
		

!!!!!!!!!!!!!!!!!!!!
! Leap Frog    !!!!!
!!!!!!!!!!!!!!!!!!!!

! Setting up initial variables
	sigma_sq = sigma*sigma	
	epsilon_Joules_molec = (epsilon_kcalmol*4.184*1000)/Na
	force_constant = 24.0*epsilon_Joules_molec
	box = l + box_dim_m ! For PBC (adding increment to avoid boundary atom overlap)
	num_tra = total_t/delta_t  !Number of trajectories
	J_to_kcalmol = 0.000239006*Na  ! adjusting Joules/molecule to kcal/mol
	half_step = delta_t/2.0
	kinetic_energy = 0.0
	
! Initial acceleration calculation before loops t = t
	call force(r, box, sigma_sq, force_constant, m, epsilon_Joules_molec, acceleration, vLJ) !plug in coordinates array, get out acceleration array


! Setup for gr in reduced units
	box_gr = box/sigma		!Reduced units
	dr = 0.05
	ro_n = N/((box_gr)**3.0)           ! Number density
	rdf_range = 5.0
	bins_rdf = int(rdf_range/dr)
	const = 4.0 * pi * ro_n/3.0
	

	allocate(h(bins_rdf), g(bins_rdf))

	do i = 1, bins_rdf
		h(i) = 0.0
		g(i) = 0.0
	enddo


! Now start the 10 000 trajectories
	open(12, file  = "evst.dat")

	do b = 1, int(num_tra)
		print*, b
		
		!  t = t + 1/2 delta t
!		 Calculating velocity from acceleration
		vel_norm(:,:) = vel_norm(:,:) + acceleration(:,:) * half_step
	
		! t = t + delta t
		! Advancing the coordinate according to new velocity

		r(:,:) = r(:,:) + vel_norm(:,:)*delta_t
		
		! gr calculations for last 1000 trajectories
		if (b.ge. 9000) then
			do i = 1, N-1
				do j = i + 1, N
				rij_array(:) = (r(:,i) - r(:,j))
				rij_array(:) = rij_array(:)/sigma	! reduced units		
				rij_array(:) = rij_array(:)-(box_gr)*ANINT(rij_array(:)/box_gr)
				rij_sq = sum(rij_array**2.0)
				rdf_count = floor(sqrt(rij_sq)/dr) + 1
					if (rdf_count .le. bins_rdf) then
					h(rdf_count) = h(rdf_count) + 2
					endif
				enddo
			enddo
		endif
			
		! Next step: t = t + delta t
		! Computing new acceleration according to new coordinate
		
		call force(r, box, sigma_sq, force_constant, m, epsilon_Joules_molec, acceleration, vLJ)
		
		! t = t + (3/2)*delta t
		! Advancing the velocity according to new velocity
		
		vel_norm(:,:) = vel_norm(:,:) + acceleration(:,:) * half_step
		
		! kinetic energy as a function of simulation step
		
		kinetic_energy = sum(0.5*m*(vel_norm(:,:)**2.0))		! Joules/1000molecules
		
		kinetic_energy = kinetic_energy*J_to_kcalmol              ! kcal/mol
		
		total_energy = kinetic_energy + vLJ
		
		write(12, *) b, vLJ, kinetic_energy, total_energy

	enddo

! Data file for plotting g(r)
	open(11, file = "gr_liq.dat")

	do k = 1, bins_rdf
	r_lo = real(k-1) * dr
	r_hi = r_lo + dr
	r_plot = (r_hi + r_lo)/2.0
	h_id = const * (r_hi**3.0 - r_lo**3.0)
	h(k) = h(k)/real(N*1000.0)     ! averaging over 1000 trajectories
	g(k) = real(h(k)/real(h_id))
	write(11,*) r_plot, g(k)
	enddo

	end program ar_md

	
	subroutine force(r, box, sigma_sq, force_constant, m, epsilon_Joules_molec, acceleration, vLJ)
		implicit none
		! calculate acceleration from coordinates of the system, calculate LJ potential
		! 1000 in array size stands for N - change it if you ever model different # of atoms
		
		real, intent(in) :: r(3, 1000)
		real, intent(in) :: box, sigma_sq, m, force_constant, epsilon_Joules_molec
		real, intent(out) :: acceleration(3, 1000), vLJ
		real :: rij_array(3), fij(3), f(3, 1000)
		real :: rij_sq, sig_over_r_sq, vij, v
		real, parameter :: epsilon_kcalmol = 0.997
		integer :: i, j
		
		f(:,:) = 0.0
		v = 0.0
		
		do i = 1, 999
			do j = i + 1, 1000
			rij_array(:) = r(:,i) - r(:,j)		
			rij_array(:) = rij_array(:)-(box)*ANINT(rij_array(:)/box)
			rij_sq = sum(rij_array**2.0)
			sig_over_r_sq = sigma_sq/rij_sq
			vij = (sig_over_r_sq**6.0 - sig_over_r_sq**3.0)	
				

			fij(:) = (2.0*sig_over_r_sq**6.0 - sig_over_r_sq**3.0)*rij_array(:)/rij_sq
	
			f(:,i) = f(:,i) + fij(:)
			f(:,j) = f(:,j) - fij(:)
			
			v = v + vij
			enddo
		enddo
		
		f = force_constant*f
		
		acceleration = f/m
		vLJ = 4.0*epsilon_kcalmol*v

		return
	end subroutine force
	
