#include "globals.h"
#include "timing.h"
#include <cmath>


__global__ void d_bondStressEnergy(int*, int*, int*, float*,
	float*, float*, float*, float*, float*, float*, int, int, int, int);

void calc_nbEnergy(void);
void calc_nbVirial(void);
void calc_bondedProps(void);
void calc_electrostatic_energy_directly(void);
void host_bonds(void);
void host_angles(void);
void cuda_collect_x(void);
void cuda_collect_charge_density_field(void);
void cuda_collect_electrostatic_potential(void);
void calc_electrostatic_energy(void);
void cuda_collect_electric_field(void);

void calc_properties() {
	// Upe is zeroed in nbEnergy
	calc_nbEnergy();

	// Ptens is zeroed in nbVirial
	calc_nbVirial();

	if (do_charges == 1) {
		// electrostatic_energy zeroed in calc_electrostatic_energy
		calc_electrostatic_energy();


		cuda_collect_electric_field();

		Upe += *electrostatic_energy;

		//calc_electrostatic_energy_directly();//comment out on systems with more than 2 particles
	}


	bond_t_in = int(time(0));
	calc_bondedProps();
	bond_t_out = int(time(0));
	bond_tot_time += bond_t_out - bond_t_in;

	//cout << "Bond Pxx: " << bondVir[0] << " tin, tout: " << bond_t_in \
	//	<< ", " << bond_t_out << endl;
}

void calc_bondedProps() {

	cuda_collect_x();
	host_bonds();

	Upe += Ubond;
	for (int i = 0; i < n_P_comps; i++) {
		Ptens[i] += bondVir[i];
	}



	host_angles();

	Upe += Uangle;
	for (int i = 0; i < n_P_comps; i++) {
		Ptens[i] += angleVir[i];
	}

	
}


void calc_nbEnergy() {
	Upe = 0.0f;

	for (Gaussian& Iter : Gausses) {
		Upe += Iter.CalcEnergy();
	}

	for (Erf& Iter : Erfs) {
		Upe += Iter.CalcEnergy();
	}

	for (GaussianErf& Iter : GaussErfs) {
		Upe += Iter.CalcEnergy();
	}

}

void calc_nbVirial() {
	for (int i = 0; i < n_P_comps; i++)
		Ptens[i] = 0.0f;

	for (Gaussian& Iter : Gausses) {
		Iter.CalcVirial();

		for (int j = 0; j < n_P_comps; j++)
			Ptens[j] += Iter.total_vir[j];
	}

	for (Erf& Iter : Erfs) {
		Iter.CalcVirial();

		for (int j = 0; j < n_P_comps; j++)
			Ptens[j] += Iter.total_vir[j];
	}


	for (GaussianErf& Iter : GaussErfs) {
		Iter.CalcVirial();

		for (int j = 0; j < n_P_comps; j++)
			Ptens[j] += Iter.total_vir[j];
	}
}

void calc_electrostatic_energy() {
	*electrostatic_energy = 0.0f;

	cuda_collect_charge_density_field();
	cuda_collect_electrostatic_potential();
	

	for (int i = 0; i < M; i++) {
		*electrostatic_energy += electrostatic_potential[i] * charge_density_field[i]; // *M; //* d_x[i];
	}
}

void calc_electrostatic_energy_directly() {
	*electrostatic_energy_direct_computation = 0.0f;

	float distance = sqrt(((x[0][0] - x[1][0]) * (x[0][0] - x[1][0])) + ((x[0][1] - x[1][1]) * (x[1][0] - x[1][1])));

	*electrostatic_energy_direct_computation += (charges[0] * charges[1] * charge_bjerrum_length) / distance;
}