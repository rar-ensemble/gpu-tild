#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <complex>
#include <vector>
#include "field_component.h"
#include "pair_style_gaussian.h"
#include "pair_style_fieldphases.h"
#include "pair_style_erf.h"
#include "pair_style_gaussian_erf.h"
#include "pair_style_charges.h"
#include "group.h"
#include "integrator.h"
#include "nlist.h"
#include "Compute.h"

#include <curand.h>
#include <curand_kernel.h>
#include <cuda_runtime.h>
#include <cufft.h>
#include <cufftXt.h>

using namespace std;

#define MAX_BONDS 3
#define MAX_ANGLES 3

#define PI 3.141592654f
#define PI2 6.2831853071f
#define PI4 12.56637061435917295384f

/// Variables for s(k) calculations
#ifndef MAIN
extern
#endif
float ** avg_sk;

#ifndef MAIN
extern
#endif
int n_avg_calc;

#ifndef MAIN
extern
#endif
int ns, Dim, ntypes, step, max_steps, * tp, * molecID, grid_freq,
traj_freq, log_freq, struc_freq, skip_steps, mem_use, device_mem_use, RAND_SEED, read_rand_seed,
Nx[3], pmeorder, M, grid_per_partic, bin_freq,
n_total_bonds, n_total_angles, nbond_types, nangle_types,
* n_bonds, * n_angles, ** bonded_to, ** bond_type,
* d_n_bonds, * d_n_angles, * d_bonded_to, * d_bond_type,
** angle_first, ** angle_mid, ** angle_end, ** angle_type,
n_gaussian_pairstyles, ** gaussian_types,
n_fieldphase_pairstyles, ** fieldphase_types, * fieldphase_dir,
* fieldphase_phase, * fieldphase_n_periods, threads, n_P_comps,
n_erf_pairstyles, ** erf_types, 
n_groups,           // Total number of group definitions.
n_integrators,      // total number of integrators
using_GJF,          // Flag to know to allocate GJF memory
n_neighbor_lists,   // Total number of neighbor lists allocated
allocate_velocities,// Flag to allocate particle velocities
n_computes,         // Total number of computes
do_charges;

#ifndef MAIN
extern
#endif
float L[3], Lh[3], V, ** x, ** xo, ** f, ** v, * mass, * Diff, * h_ns_float,
delt, * bond_k, * bond_req, Ubond, * angle_k, * angle_theta_eq, Uangle,
dx[3], * tmp, * tmp2, * all_rho, gvol, noise_mag,
Upe, * Ptens, * partic_bondE, * partic_bondVir, * bondVir,
* gaussian_prefactor, * gaussian_sigma, * U_gaussian,
* gaussian_final_prefactor,
* fieldphase_prefactor, * U_fieldphase,
* erf_Rp, * erf_Ao, * erf_Xi, * U_erf,
* charges, charge_bjerrum_length, charge_smearing_length,
* charge_density_field, * electrostatic_energy, 
* electrostatic_potential, * electric_field,
* electrostatic_energy_direct_computation;

#ifndef MAIN
extern
#endif
string dump_name, input_file;

#ifndef MAIN
extern
#endif
complex<float> I, * k_tmp ;


#ifndef MAIN
extern
#endif
float* d_x, * d_v, * d_f, * d_L, * d_Lh,
* d_bond_req, * d_bond_k, * d_3n_tmp,
* d_tmp, * d_tmp2, * d_dx, * d_all_rho,
* d_grid_W, * d_all_fx, * d_all_fy, * d_all_fz,
* d_all_fx_charges, * d_all_fy_charges, * d_all_fz_charges,
* d_bondE, * d_bondVir,
* d_xo, * d_prev_noise, * d_mass, * d_Diff,
* d_charges, * d_charge_density_field, * d_electrostatic_potential, * d_electric_field;

#ifndef MAIN
extern
#endif
cufftComplex* d_cpx1, * d_cpx2, * cpx1, * cpx2;

#ifndef MAIN
extern
#endif
bool* d_nan;

#ifndef MAIN
extern
#endif
cufftHandle fftplan;

#ifndef MAIN
extern
#endif
int d_ns, * d_typ, * d_Nx, ns_Block, ns_Grid, M_Block, M_Grid,
* d_grid_inds;

#ifndef MAIN
extern
#endif
cudaError_t cudaReturn;

#ifndef MAIN
extern
#endif
FieldComponent* Components;

#ifndef MAIN
extern
#endif
vector<Gaussian> Gausses;

#ifndef MAIN
extern
#endif
vector<FieldPhase> Fields;

#ifndef MAIN
extern
#endif
vector<Erf> Erfs;

#ifndef MAIN
extern
#endif
vector<Group> Groups;

#ifndef MAIN
extern
#endif
vector<Integrator> Integrators;

#ifndef MAIN
extern
#endif
vector<NList> NLists;

#ifndef MAIN
extern
#endif
vector<GaussianErf> GaussErfs;

#ifndef MAIN
extern
#endif
Charges* Charge_Pairstyle;

#ifndef MAIN
extern
#endif
vector<Compute> Computes;

#ifndef MAIN
extern
#endif
curandState* d_states;

void die(const char*);
void check_cudaError(const char*);
float pbc_mdr2(float*, float*, float*);
float get_k(int, float*, int);
void get_r(int, float*);
float integ_trapPBC(float*);
void write_grid_data(const char*, float*);
