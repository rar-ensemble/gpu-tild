#include "globals.h"
#include "timing.h"
#include <string>
#include <sstream>
#include <cmath>

void allocate_device_memory(void);
void allocate_device_potentials(void);
void send_box_params_to_device(void);
void send_3n_to_device(float**, float*);
void allocate_grid_memory(void);
void read_input(void);
void cuda_collect_x(void);
void init_binary_output(void);
__global__ void init_dev_rng(unsigned int, curandState*, int);
__global__ void d_assignValueR(float*, float, int, int);



void initialize() {
	step = 0;
	mem_use = 0;
	
	// Global unit complex value
	I = complex<float>(0.0f, 1.0f);
	
	
	read_input();
	cout << "Input file read!" << endl;


	M = 1;
	grid_per_partic = 1;
	gvol = 1.f;
	for (int j = 0; j < Dim; j++) {
		M *= Nx[j];
		dx[j] = L[j] / float(Nx[j]);
		gvol *= dx[j];
		grid_per_partic *= (pmeorder + 1);
	}
	
	// Define the number of pressure tensor components
	// if Dim == 2: 0=xx, 1=yy, 2=xy
	// If Dim == 3: 0=xx, 1=yy, 2=zz, 3=xy, 4=xz, 5=yz
	if (Dim == 2)
		n_P_comps = 3;
	else if (Dim == 3)
		n_P_comps = 6;

	print_tot_time = 0;
	bond_tot_time = 0;
	compute_tot_time = 0;

	device_mem_use = 0;
	allocate_grid_memory();
	cout << "Grid memory allocated! " << endl;

	allocate_device_memory();



	// Sends box geometry, bond topology, other
	// one-time communication information
	send_box_params_to_device();
	cout << "Box parameters sent to device" << endl;


	allocate_device_potentials();
	cout << "Device potentials initialized" << endl;

	// Send initial box coordinates
	send_3n_to_device(x, d_x);
	cout << "Initial box coordinates sent!" << endl;
	

	cout << "Device memory allocated: " << float(device_mem_use) / powf(10.0f, 6)
		<< " MB\n";

	cuda_collect_x();

	if (bin_freq != 0)
		init_binary_output();

	
}







void allocate_device_memory() {
	device_mem_use = 0;


	ns_Block = threads;
	ns_Grid = (int)ceil((float)(ns) / ns_Block);
	printf("ns: %d, ns_Grid: %d\n", ns, ns_Grid);


	M_Block = threads;
	M_Grid = (int)ceil((float)(M) / M_Block);
	cout << "M: " << M << ", M_Grid: " << M_Grid << endl;


	int size = ns * Dim * sizeof(float);
	cudaMalloc(&d_x, size);
	cudaMalloc(&d_f, size);
	cudaMalloc(&d_xo, size);
	cudaMalloc(&d_v, size);
	cudaMalloc(&d_3n_tmp, size);

	device_mem_use += size * 5;

	if (using_GJF) {
		cudaMalloc(&d_xo, size);
		cudaMalloc(&d_prev_noise, size);

		d_assignValueR << <ns_Grid, ns_Block >> > (d_prev_noise, 0.0f, Dim, ns);

		cudaReturn = cudaGetLastError();
		if (cudaReturn != cudaSuccess) {
			char cherror[90];
			sprintf(cherror, "Cuda failed with error \"%s\" while allocating memory\n", cudaGetErrorString(cudaReturn));
			die(cherror);
		}

		device_mem_use += 2 * size;
	}



	// Initialize random number seeds on each thread
	cudaMalloc(&d_states, ns * Dim * sizeof(curandState));
	init_dev_rng<<<ns_Grid, ns_Block>>>(RAND_SEED, d_states, ns);

	device_mem_use += ns * Dim * sizeof(curandState);

	cudaMalloc(&d_mass, ntypes * sizeof(float));
	cudaMalloc(&d_Diff, ntypes * sizeof(float));

	// Allocate box parameters on device
	cudaMalloc(&d_L, 3 * sizeof(float));
	cudaMalloc(&d_Lh, 3 * sizeof(float));

	device_mem_use += sizeof(float) * (3 + 3);

	cudaMalloc(&d_typ, ns * sizeof(int));
	device_mem_use += sizeof(int) * ns;

	// Bond parameters on device
	cudaMalloc(&d_n_bonds, ns * sizeof(int));
	cudaMalloc(&d_bonded_to, ns * MAX_BONDS * sizeof(int));
	cudaMalloc(&d_bond_type, ns * MAX_BONDS * sizeof(int));
	cudaMalloc(&d_bond_req, nbond_types * sizeof(float));
	cudaMalloc(&d_bond_k, nbond_types * sizeof(float));

	device_mem_use += sizeof(float) * (ns + ns * MAX_BONDS * 2 + 2 * nbond_types);
	
	cudaMalloc(&d_bondE, ns * sizeof(float));
	cudaMalloc(&d_bondVir, ns * n_P_comps * sizeof(float));

	device_mem_use += ns * (n_P_comps + 1) * sizeof(float);

	if (n_total_angles > 0) {
		die("Angles not yet set up on device!");
	}

	cudaMalloc(&d_Nx, 3 * sizeof(int));
	device_mem_use += sizeof(int) * 3;

	cudaMalloc(&d_dx, 3 * sizeof(float));
	cudaMalloc(&d_tmp, M * sizeof(float));
	cudaMalloc(&d_tmp2, M * sizeof(float));
	cudaMalloc(&d_all_rho, ntypes * M * sizeof(float));

	if (do_charges == 1) {
		cudaMalloc(&d_charge_density_field, M * sizeof(float));
		cudaMalloc(&d_electrostatic_potential, M * sizeof(float));
		cudaMalloc(&d_charges, ns * sizeof(float));
		cudaMalloc(&d_electric_field, M * Dim * sizeof(float));
		//cudaMalloc(&d_charges_over_M, ns * sizeof(float));
	}

	cudaMalloc(&d_all_fx, ntypes * M * sizeof(float));
	cudaMalloc(&d_all_fy, ntypes * M * sizeof(float));
	device_mem_use += sizeof(float) * (2 * M + ntypes * M * 3 + 3);

	if (Dim == 3) {
		cudaMalloc(&d_all_fz, ntypes * M * sizeof(float));
		device_mem_use += sizeof(float*) * (ntypes * M);
	}

	if (do_charges == 1) {
		cudaMalloc(&d_all_fx_charges, M * sizeof(float));
		cudaMalloc(&d_all_fy_charges, M * sizeof(float));
		device_mem_use += sizeof(float) * (2 * M);

		if (Dim == 3) {
			cudaMalloc(&d_all_fz_charges, M * sizeof(float));
			device_mem_use += sizeof(float*) * (M);
		}
	}

	cudaMalloc(&d_cpx1, M * sizeof(cufftComplex));
	cudaMalloc(&d_cpx2, M * sizeof(cufftComplex));
	device_mem_use += sizeof(cufftComplex) * (2 * M);

	cudaMalloc(&d_nan, sizeof(bool));
	device_mem_use += sizeof(bool);

	cudaMalloc(&d_grid_W, ns * grid_per_partic * sizeof(float));
	cudaMalloc(&d_grid_inds, ns * grid_per_partic * sizeof(int));
	device_mem_use += ns * grid_per_partic * (sizeof(float) + sizeof(int));

	if (Dim == 2)
		cufftPlan2d(&fftplan, Nx[1], Nx[0], CUFFT_C2C);
	else if (Dim == 3)
		cufftPlan3d(&fftplan, Nx[2], Nx[1], Nx[0], CUFFT_C2C);

}

void allocate_device_potentials(void) {

	for (Gaussian& Iter : Gausses) {
		Iter.Initialize();
	}


	for (FieldPhase& Iter : Fields) {
		Iter.Initialize();
	}

	for (Erf& Iter : Erfs) {
		Iter.Initialize();
	}

	for (GaussianErf& Iter : GaussErfs) {
		Iter.Initialize();
	}
	
	Charge_Pairstyle = new Charges;
}


void allocate_grid_memory(void) {
	tmp = (float*)calloc(M, sizeof(float*));
	tmp2 = (float*)calloc(M, sizeof(float*));
	cpx1 = (cufftComplex*)calloc(M, sizeof(cufftComplex));
	cpx2 = (cufftComplex*)calloc(M, sizeof(cufftComplex));
	k_tmp = (complex<float>*) calloc(M, sizeof(complex<float>));

	if (do_charges == 1) {
		electrostatic_potential = (float*)calloc(M, sizeof(float));
		charge_density_field = (float*)calloc(M, sizeof(float));
		electric_field = (float*)calloc(M * Dim, sizeof(float));
	}

	all_rho = (float*)calloc(M * ntypes, sizeof(float));
	
	// Allocate memory for structure factor calculation
	if (struc_freq > 0) {
		avg_sk = (float**)calloc(ntypes, sizeof(float*));
		for (int i = 0; i < ntypes; i++) {
			avg_sk[i] = (float*)calloc(M, sizeof(float));
			for (int j = 0; j < M; j++) {
				avg_sk[i][j] = 0;
			}
		}
	}

	Components = new FieldComponent[ntypes];
	for (int i = 0; i < ntypes; i++)
		Components[i].Initialize(M);

  for ( int i=0 ; i < n_computes ; i++ ) 
    Computes[i].allocStorage();

}



void allocate_host_particles() {

	x = (float**)calloc(ns, sizeof(float*));
	xo = (float**)calloc(ns, sizeof(float*));
	v = (float**)calloc(ns, sizeof(float*));
	f = (float**)calloc(ns, sizeof(float*));
	if (do_charges == 1) charges = (float*)calloc(ns, sizeof(float));

	for (int i = 0; i < ns; i++) {
		x[i] = (float*)calloc(Dim, sizeof(float));
		xo[i] = (float*)calloc(Dim, sizeof(float));
		v[i] = (float*)calloc(Dim, sizeof(float));
		f[i] = (float*)calloc(Dim, sizeof(float));
	}

	h_ns_float = (float*)calloc(ns * Dim, sizeof(float));

	partic_bondE = (float*)calloc(ns, sizeof(float));
	partic_bondVir = (float*)calloc(ns * n_P_comps, sizeof(float));
	bondVir = (float*)calloc(n_P_comps, sizeof(float));

	tp = (int*)calloc(ns, sizeof(int));
	molecID = (int*)calloc(ns, sizeof(int));

	mass = (float*)calloc(ntypes, sizeof(float));
	Diff = (float*)calloc(ntypes, sizeof(float));

	// Set default diffusivities
	for (int i = 0; i < ntypes; i++)
		Diff[i] = 1.0;

	Ptens = (float*)calloc(n_P_comps, sizeof(float));

	// NOTE: Assumes that a particle is bonded to a maximum 
	// of MAX_BONDS particles
	n_bonds = (int*)calloc(ns, sizeof(int));
	n_angles = (int*)calloc(ns, sizeof(int));
	bonded_to = (int**)calloc(ns, sizeof(int*));
	bond_type = (int**)calloc(ns, sizeof(int*));
	angle_first = (int**)calloc(ns, sizeof(int*));
	angle_mid = (int**)calloc(ns, sizeof(int*));
	angle_end = (int**)calloc(ns, sizeof(int*));
	angle_type = (int**)calloc(ns, sizeof(int*));

	for (int i = 0; i < ns; i++) {
		bonded_to[i] = (int*)calloc(MAX_BONDS, sizeof(int));
		bond_type[i] = (int*)calloc(MAX_BONDS, sizeof(int));

		angle_first[i] = (int*)calloc(MAX_ANGLES, sizeof(int));
		angle_mid[i] = (int*)calloc(MAX_ANGLES, sizeof(int));
		angle_end[i] = (int*)calloc(MAX_ANGLES, sizeof(int));
		angle_type[i] = (int*)calloc(MAX_ANGLES, sizeof(int));
	}

	mem_use += ns * (2 * MAX_BONDS + 4 * MAX_ANGLES + 2) * sizeof(int);

	bond_k = (float*)calloc(nbond_types, sizeof(float));
	bond_req = (float*)calloc(nbond_types, sizeof(float));
	angle_k = (float*)calloc(nangle_types, sizeof(float));
	angle_theta_eq = (float*)calloc(nangle_types, sizeof(float));

	if (do_charges == 1) {
		electrostatic_energy = (float*)calloc(1, sizeof(float));
		electrostatic_energy_direct_computation = (float*)calloc(1, sizeof(float));
		//charges_over_M = (float*)calloc(ns, sizeof(float));
	}

}


//__global__  void nan_kernel_cpx(cufftComplex* d, int len, bool* res) {
//	int idx = threadIdx.x + blockDim.x * blockIdx.x;
//	if (idx < len)
//		if (isnan(d[idx].x) || isnan(d[idx].y)) *res = true;
//}

//bool checknan_cpx(cufftComplex* d, int len) {
//	bool h_nan = false;
//	cudaMemcpy(d_nan, &h_nan, sizeof(bool), cudaMemcpyHostToDevice);
//	nan_kernel_cpx <<<M_Grid, M_Block>>> (d, len, d_nan);
//	cudaMemcpy(&h_nan, d_nan, sizeof(bool), cudaMemcpyDeviceToHost);
//	return h_nan;
//}
