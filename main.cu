#define MAIN
#include "globals.h"
#include "timing.h"
#include <vector>
#include <fstream>

void forces(void);
void update_pairstyles(void);
void calc_properties(void);
void initialize(void);
void write_lammps_traj(void);
void cuda_collect_x(void);
void cuda_collect_f(void);
void cuda_collect_rho(void);
void write_binary(void);
void write_struc_fac(void);
void write_grid_data(const char*, float*);
void write_kspace_data(const char*, complex<float>*);
void write_kspace_cudaComplex(const char*, cufftComplex*);
__global__ void d_prepareDensity(int, float*, cufftComplex*, int);
int print_timestep();
ofstream dout;
void unstack_like_device(int id, int* nn);

__global__ void cu_random_posits(float*, float*, int, int, curandState*);


__global__ void d_real2complex(float*, cufftComplex*, int);
__global__ void d_complex2real(cufftComplex*, float*, float*, int);
__global__ void d_make_step(cufftComplex* , float*, float*, int*, int, int);
__global__ void d_multiplyComplex(cufftComplex*, cufftComplex*,
	cufftComplex*, int);



int main(int argc, char** argv)
{
	main_t_in = int(time(0));
	init_t_in = main_t_in;
	std::vector<std::string> string_vec;
	for (int i = 0; i < argc; i++)
	{
		std::string arg = argv[i];
		string_vec.push_back(arg);
	}

	for (size_t i = 1; i < string_vec.size(); ++i) {
		if (string_vec[i] == "-in") {
			input_file = string_vec[++i];
		}
	}

	initialize();


	// Write initial positions to lammpstrj file
	cuda_collect_x();
	write_lammps_traj();

	dout.open("data.dat");
	dout << "# step Upe Ubond ";
	if (Dim == 2)
		dout << "Pxx Pyy Pxy ";
	else if (Dim == 3)
		dout << " Pxx Pyy Pzz Pxy Pxz Pyz ";
	for (int i = 0; i < Gausses.size(); i++) 
		dout << "Ugauss[" << i << "] ";
	for (int i = 0; i < Erfs.size(); i++) 
		dout << "Uerf[" << i << "] ";
	for (int i = 0; i < GaussErfs.size(); i++) 
		dout << "UGuassErf[" << i << "] ";
	
	dout << endl;
	
	forces();
	cudaDeviceSynchronize();

	cuda_collect_rho();
	cuda_collect_x();
	calc_properties();
	if (grid_freq > 0) {
		print_t_in = int(time(0));
		cudaDeviceSynchronize();

		for (int i = 0; i < ntypes; i++) {
			char nm[30];
			sprintf(nm, "rho%d.dat", i);
			write_grid_data(nm, Components[i].rho);
		}
			
		print_t_out = int(time(0));
		print_tot_time += print_t_out - print_t_in;
	}

	int die_flag = 0;
	die_flag = print_timestep();

	init_t_out = int(time(0));



	///////////////////////////////////////
	// BEGINNING OF MAIN SIMULATION LOOP //
	///////////////////////////////////////

	for (step = 1; step <= max_steps; step++) {
		
		for (int i = 0; i < n_integrators; i++)
			Integrators[i].Integrate_1();
		check_cudaError("Integrator step 1");


		forces();


		for (int i = 0; i < n_integrators; i++)
		  Integrators[i].Integrate_2();
	
		check_cudaError("Integrator step 2");


		
		// Run computes
		for (int i = 0; i < n_computes; i++) {
			if (step > Computes[i].compute_wait && step % Computes[i].compute_freq == 0) {
				
				//cout << "entering compute " << i;
				Computes[i].doCompute();
				//cout << " done!" << endl;
				check_cudaError("Compute");
			}
		}

		// I/O blocks //
		if (traj_freq > 0 && step % traj_freq == 0) {
			print_t_in = int(time(0));
			cudaDeviceSynchronize();

			cuda_collect_x();
			write_lammps_traj();
			print_t_out = int(time(0));
			print_tot_time += print_t_out - print_t_in;
		}

		if (grid_freq > 0 && step % grid_freq == 0) {
			print_t_in = int(time(0));
			cudaDeviceSynchronize();

			cuda_collect_rho();
			for (int i = 0; i < ntypes; i++) {
				char nm[30];
				sprintf(nm, "rho%d.dat", i);
				write_grid_data(nm, Components[i].rho);
			}

			print_t_out = int(time(0));
			print_tot_time += print_t_out - print_t_in;
		}

		if (bin_freq > 0 && step % bin_freq == 0) {
			print_t_in = int(time(0));
			cudaDeviceSynchronize();

			cuda_collect_rho();
			cuda_collect_x();

			write_binary();
			print_t_out = int(time(0));
			print_tot_time += print_t_out - print_t_in;
		}


		// Write to log file, write compute results
		if (step % log_freq == 0) {
			print_t_in = int(time(0));
			cudaDeviceSynchronize();

			calc_properties();

			die_flag = print_timestep();

			for (int i = 0; i < n_computes; i++)
				if (step > Computes[i].compute_wait)
					Computes[i].writeResults(i);

			print_t_out = int(time(0));
			print_tot_time += print_t_out - print_t_in;

			if (die_flag) {
				break;
			}


		}
		
	

		// Finalize time step //
		update_pairstyles();

	}// main loop over steps



	// Write resume frame and finish //
	cuda_collect_x();
	write_lammps_traj();

	main_t_out = int(time(0));
	int dt = main_t_out - main_t_in;
	cout << "Total run time: " << dt / 60 << "m" << dt % 60 << "sec" << endl;
	
	dt = init_t_out - init_t_in;
	cout << "Total init time: " << dt / 60 << "m" << dt % 60 << "sec" << endl;
	
	dt = bond_tot_time;
	cout << "Bond E, P props on host: " << dt / 60 << "m" << dt % 60 << "sec" << endl;

	dt = print_tot_time;
	cout << "I/O + Comm time: " << dt / 60 << "m" << dt % 60 << "sec" << endl;

	dt = compute_tot_time;
	cout << "Computes time: " << dt / 60 << "m" << dt % 60 << "sec" << endl;

	return 0;

}


int print_timestep() {
	int die_flag = 0;
	cout << "Step " << step << " of " << max_steps << " ";

	if (do_charges == 1) {
		cout << " Electrostatic Energy: " << *electrostatic_energy;
	}

	cout << " U/V: " << Upe / V << \
		" Ubond: " << Ubond << \
		" Pdiags: " << Ptens[0] << " " << Ptens[1] << " ";

	if (Dim == 3)
		cout << Ptens[2] << " ";

	dout << step << " " << Upe << " " << Ubond << " ";
	for (int i = 0; i < n_P_comps; i++)
		dout << Ptens[i] << " ";

	
	if (!Gausses.empty()) {
		cout << "Ugauss: ";
			for (Gaussian& Iter : Gausses){
			cout << Iter.energy << " ";
			dout << Iter.energy  << " ";
			if (std::isnan(Iter.energy))
				die_flag = 1;
		}
	}

	if (!Erfs.empty()) {
	cout << "Uerf: ";
		for (Erf& Iter: Erfs) {
			cout << Iter.energy << " ";
			dout << Iter.energy << " ";
			if (std::isnan(Iter.energy))
				die_flag = 1;
		}
	}

	if (!GaussErfs.empty()) {
		cout << "UGaussErf: ";
		for (GaussianErf& Iter : GaussErfs) {
			cout << Iter.energy << " ";
			dout << Iter.energy << " ";
			if (std::isnan(Iter.energy))
				die_flag = 1;
		}
	}

	if (!Fields.empty()) {
		cout << "UFieldPhase: ";
		for (FieldPhase& Iter : Fields) {
			cout << Iter.energy << " ";
			dout << Iter.energy << " ";
			if (std::isnan(Iter.energy))
				die_flag = 1;
		}
	}

    dout << endl;
	cout<<endl;
	return die_flag;
}
