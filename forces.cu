#include "globals.h"


__global__ void d_zero_all_ntyp(float*, int, int);
__global__ void d_zero_all_directions_charges(float*, int);
__global__ void d_bonds(int*, int*, int*, float*,
	float*, float*, float*, float*, float*, int, int, int);
__global__ void d_angles(const float*, float*, const float*,
    const float*, const int*, const int*, const int*, const int*,
    const int*, const float*, const float*, const int, const int,
    const int);

__global__ void d_prep_components(float*, float*, float*,
	const int, const int, const int);

__global__ void d_zero_particle_forces(float*, int, int);
__global__ void d_charge_grid_charges(float*, float*, int*, int*, float*, const int*,
	const float*, const float, const int, const int, const int, const int, float*, float*, int);
__global__ void d_charge_grid(float*, float*, int*, int*, float*, const int*,
	const float*, const float, const int, const int, const int, const int);
__global__ void d_real2complex(float*, cufftComplex*, int M);
__global__ void d_add_grid_forces2D(float*, const float*,
	const float*, const float*, const float*,
	const int*, const int*, const float,
	const int, const int, const int, const int);
__global__ void d_add_grid_forces_charges_2D(float*, const float*,
	const float*, const float*, const float*,
	const int*, const float, const float*,
	const int, const int, const int, const int);
__global__ void d_add_grid_forces3D(float*, const float*,
	const float*, const float*, const float*, const float*,
	const int*, const int*, const float,
	const int, const int, const int, const int);
__global__ void d_add_grid_forces_charges_3D(float*, const float*,
	const float*, const float*, const float*, const float*,
	const int*, const float, const float*,
	const int, const int, const int, const int);

void cuda_collect_rho(void);
void write_grid_data(const char*, float*);
void write_lammps_traj(void);
void cuda_collect_x(void);
__global__ void d_make_dens_step(float*, float*, float*, int*, int, int, int);


void forces() {
	
	// Zeros the types*M density field
	d_zero_all_ntyp<<<M_Grid, M_Block>>>(d_all_rho, M, ntypes);

	d_zero_all_ntyp<<<M_Grid, M_Block>>>(d_all_fx, M, ntypes);
	d_zero_all_ntyp<<<M_Grid, M_Block>>>(d_all_fy, M, ntypes);
	if ( Dim == 3 ) d_zero_all_ntyp<<<M_Grid, M_Block>>>(d_all_fz, M, ntypes);

	if (do_charges == 1) {
		d_zero_all_directions_charges<<<M_Grid, M_Block>>>(d_charge_density_field, M);
		d_zero_all_directions_charges<<<M_Grid, M_Block>>>(d_electrostatic_potential, M);

		d_zero_all_directions_charges<<<M_Grid, M_Block>>>(d_all_fx_charges, M);
		d_zero_all_directions_charges<<<M_Grid, M_Block>>>(d_all_fy_charges, M);

		if (Dim == 3) d_zero_all_directions_charges<<<M_Grid, M_Block>>>(d_all_fz_charges, M);
	}

	cudaReturn = cudaGetLastError();
	if (cudaReturn != cudaSuccess) {
		char cherror[90];
		sprintf(cherror, "Cuda failed with error \"%s\" while forces zeroed \n", cudaGetErrorString(cudaReturn));
		die(cherror);
	}

	// Fills the ntypes*M density field
	// and Fills d_charge_density_field if charges is flagged
	if (do_charges == 1) {
		d_charge_grid_charges<<<ns_Grid, ns_Block>>>(d_x, d_grid_W,
			d_grid_inds, d_typ, d_all_rho, d_Nx, d_dx,
			V, ns, pmeorder, M, Dim, d_charge_density_field, d_charges, do_charges);
	}
	else {
		d_charge_grid<<<ns_Grid, ns_Block>>>(d_x, d_grid_W,
			d_grid_inds, d_typ, d_all_rho, d_Nx, d_dx,
			V, ns, pmeorder, M, Dim);
	}


	cudaReturn = cudaGetLastError();
	if (cudaReturn != cudaSuccess) {
		char cherror[90];
		sprintf(cherror, "Cuda failed with error \"%s\" while d_charge_grid ran\n", cudaGetErrorString(cudaReturn));
		die(cherror);
	}



	// Zeros the grid forces and copies the density into
	// its class structure.
	for (int i = 0; i < ntypes; i++) {
		d_prep_components<<<M_Grid, M_Block>>>
			(Components[i].d_force, Components[i].d_rho, d_all_rho, i, M, Dim);
	}

	cudaReturn = cudaGetLastError();
	if (cudaReturn != cudaSuccess) {
		char cherror[90];
		sprintf(cherror, "Cuda failed with error \"%s\" while d_prep_components ran\n", cudaGetErrorString(cudaReturn));
		die(cherror);
	}


  // Gaussian potential forces
	for (Gaussian& Iter : Gausses) {
		Iter.CalcForces();
	}

  // Erf potential forces
	for (Erf& Iter : Erfs) {
		Iter.CalcForces();
	}
	//d_make_dens_step<<<M_Grid, M_Block>>>(d_all_rho, d_L, d_dx, d_Nx, Dim, M, ntypes);
	

  // FieldPhase forces
	for (FieldPhase& Iter : Fields) {
		Iter.CalcForces();
	}

  // Gauss/Erf potential forces
	for (GaussianErf& Iter : GaussErfs) {
		Iter.CalcForces();
	}
	
  // Electrostatic forces
	if (do_charges == 1) (*Charge_Pairstyle).CalcCharges();
	

	// Zeros the forces on the particles
	d_zero_particle_forces<<<ns_Grid, ns_Block>>>(d_f, ns, Dim);

	cudaReturn = cudaGetLastError();
	if (cudaReturn != cudaSuccess) {
		char cherror[90];
		sprintf(cherror, "Cuda failed with error \"%s\" while d_zero_particle_forces ran\n", cudaGetErrorString(cudaReturn));
		die(cherror);
	}


	// Accumulates forces from grid onto particles
	if (Dim == 2) {
		if (do_charges == 1) {
			d_add_grid_forces_charges_2D<<<ns_Grid, ns_Block>>>(d_f, d_all_fx_charges,
				d_all_fy_charges, d_charge_density_field, d_grid_W, d_grid_inds, gvol,
				d_charges, grid_per_partic, ns, M, Dim);
		}
		
		d_add_grid_forces2D<<<ns_Grid, ns_Block>>>(d_f, d_all_fx,
			d_all_fy, d_all_rho, d_grid_W, d_grid_inds, d_typ, gvol,
			grid_per_partic, ns, M, Dim);
		
	}

	else if (Dim == 3) {
		if (do_charges == 1) {
			
			d_add_grid_forces_charges_3D<<<ns_Grid, ns_Block>>>(d_f, d_all_fx_charges,
				d_all_fy_charges, d_all_fz_charges, d_charge_density_field, d_grid_W, d_grid_inds, gvol,
				d_charges, grid_per_partic, ns, M, Dim);
		}
		
		d_add_grid_forces3D<<<ns_Grid, ns_Block>>>(d_f, d_all_fx,
			d_all_fy, d_all_fz, d_all_rho, d_grid_W, d_grid_inds, d_typ, gvol,
			grid_per_partic, ns, M, Dim);
		
	}

	cudaReturn = cudaGetLastError();
	if (cudaReturn != cudaSuccess) {
		char cherror[90];
		sprintf(cherror, "Cuda failed with error \"%s\" while d_add_grid_forces ran\n", cudaGetErrorString(cudaReturn));
		die(cherror);
	}




	// Accumulates bonded forces
	if (n_total_bonds > 0) {
		d_bonds<<<ns_Grid, ns_Block>>>(d_n_bonds, d_bonded_to,
			d_bond_type, d_bond_req, d_bond_k, d_x, d_f,
			d_L, d_Lh, ns, MAX_BONDS, Dim);

    check_cudaError("bond forces");
	}

    if ( n_total_angles > 0 ) {
        d_angles<<<ns_Grid, ns_Block>>>(d_x, d_f, d_angle_k,
            d_angle_theta_eq, d_n_angles, d_angle_type,
            d_angle_first,d_angle_mid, d_angle_end, d_L, d_Lh,
            ns, MAX_ANGLES, Dim);
        
        check_cudaError("angle forces");
    }

    for ( int i=0 ; i<n_extra_forces ; i++ ) 
      ExtraForces[i].AddExtraForce();



}
