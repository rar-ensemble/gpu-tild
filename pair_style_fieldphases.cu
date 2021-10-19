#include "globals.h"
#include "pair_style_fieldphases.h"
#include "device_utils.cuh"

/*  Field Phase pairstyle
    The point of this type of interaction is to seed 
    specific fields and bias the system to form a target
    structure, such as lamellae, cylinders, etc. The idea 
    is to replace the chi-like interactions between two species
    A and B with the a static field w(r), where A experiences
    +w(r) and B experiences -w(r). This combined with
    kappa-like self repulsions should help the system form
    a target phase. */


__global__ void init_device_fieldphase(float*, float*,
    const int, const float, const int,
    const int, const float*, const float*,
    const int, const int*, const int);
__global__ void d_real2complex(float* , cufftComplex* , int);
__global__ void d_divideByDimension(cufftComplex*, int);
__global__ void d_prepareFieldForce(cufftComplex*, cufftComplex*,
    const float*, const float*, const int, const int, const int,
    const int, const int);
__global__ void d_accumulateGridForce(cufftComplex*, float*, float*,
    const int, const int);

void FieldPhase::Initialize() {
	Initialize_FieldPhase(initial_prefactor, phase, dir, n_periods, M, type1, type2);
}

void FieldPhase::Initialize_FieldPhase(float Ao, int phase,
    int dir, int n_periods, int alloc_size, int typ_A, int typ_B) {

    Initialize_PairStyle(alloc_size, typ_A, typ_B);

    printf("Setting up FieldPhase pair style..."); fflush(stdout);


    init_device_fieldphase<<<M_Grid, M_Block>>>(this->d_u, this->d_f,
        phase, Ao, dir, n_periods, d_L, d_dx, M, d_Nx, Dim);

    cudaMemcpy(this->u, this->d_u, M * sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(this->u, this->d_u, M * sizeof(float), cudaMemcpyDeviceToHost); 

    write_grid_data("fieldphase_map.dat", this->u);

    printf("done!\n"); fflush(stdout);
}

void FieldPhase::CalcForces() {


    // X-component of the force
    d_prepareFieldForce<<<M_Grid, M_Block>>>(d_cpx1, d_cpx2,
        d_all_rho, this->d_f, type1, type2, 0, Dim, M);
    d_accumulateGridForce<<<M_Grid, M_Block>>>(d_cpx1,
        d_all_rho, d_all_fx, type1, M);
    d_accumulateGridForce<<<M_Grid, M_Block>>>(d_cpx2,
        d_all_rho, d_all_fx, type2, M);




    // Y-component of the force
    d_prepareFieldForce<<<M_Grid, M_Block>>>(d_cpx1, d_cpx2,
        d_all_rho, this->d_f, type1, type2, 1, Dim, M);
    d_accumulateGridForce<<<M_Grid, M_Block>>>(d_cpx1,
        d_all_rho, d_all_fy, type1, M);
    d_accumulateGridForce<<<M_Grid, M_Block>>>(d_cpx2,
        d_all_rho, d_all_fy, type2, M);


    if (Dim == 3) {
        // Z-component of the force
        d_prepareFieldForce<<<M_Grid, M_Block>>>(d_cpx1, d_cpx2,
            d_all_rho, this->d_f, type1, type2, 2, Dim, M);
        d_accumulateGridForce<<<M_Grid, M_Block>>>(d_cpx1,
            d_all_rho, d_all_fz, type1, M);
        d_accumulateGridForce<<<M_Grid, M_Block>>>(d_cpx2,
            d_all_rho, d_all_fz, type2, M);
    }


}

FieldPhase::FieldPhase() {

}

FieldPhase::~FieldPhase() {

}



// phase = 0: Lamellar
// phase = 1: BCC
// phase = 2: CYL
// phase = 3: GYR
__global__ void init_device_fieldphase(float* ur, float* fr,
	const int phase, const float Ao, const int dir,
	const int n_periods, const float* dL, const float* dx,
	const int M, const int* Nx, const int Dim) {


	const int ind = blockIdx.x * blockDim.x + threadIdx.x;
	if (ind >= M)
		return;

	float r[3];

	d_get_r(ind, r, Nx, dx, Dim);

	// Initial lamellar fields, forces
	if (phase == 0) {
		float sin_arg = 2.0f * PI * float(n_periods) / dL[dir];
		ur[ind] = Ao * sin(sin_arg * r[dir]);

		for (int j = 0; j < Dim; j++) {
			if (j == dir)
				fr[ind * Dim + j] = Ao * sin_arg * cos(sin_arg * r[dir]);
			else
				fr[ind * Dim + j] = 0.f;
		}
	}

	// CYL phase
	// Assumes the same number of periods in both directions of the 
	// hexagonal plane of the cylinders. 
	// dir is interpretted as the direction of the cylinders
	// For 2D, this should be set to 2
	else if (phase == 2) {
		float dim1_arg, dim2_arg;
		int dim1 = 0, dim2 = 1;

		if (dir == dim1)
			dim1 = 2;
		else if (dir == dim2)
			dim2 = 2;

		dim1_arg = 2.0f * PI * float(n_periods) / dL[dim1];
		dim2_arg = 2.0f * PI * float(n_periods) / dL[dim2];

		ur[ind] = Ao * cos(dim1_arg * r[dim1]) * cos(dim2_arg * r[dim2]);

		if (Dim == 3)
			fr[ind * Dim + dir] = 0.0f;

		fr[ind * Dim + dim1] = Ao * dim1_arg * sin(dim1_arg * r[dim1]) * cos(dim2_arg * r[dim2]);
		fr[ind * Dim + dim2] = Ao * dim2_arg * sin(dim2_arg * r[dim2]) * cos(dim1_arg * r[dim1]);

	}

	// GYR phase
	// Assumes same number of periods in each direction
	else if (phase == 3) {
		if (Dim != 3) {
			return;
		}

		float args[3], cos_dir[3], sin_dir[3];
		for (int j = 0; j < Dim; j++) {
			args[j] = 2.0f * PI * float(n_periods) / dL[j];

			cos_dir[j] = cos(args[j] * r[j]);
			sin_dir[j] = sin(args[j] * r[j]);
		}

		float e_term = sin_dir[0] * cos_dir[1] + sin_dir[1] * cos_dir[2]
			+ sin_dir[2] * cos_dir[0];

		ur[ind] = Ao * (e_term * e_term - 1.44);

		fr[ind * Dim + 0] = -Ao * e_term * args[0] *
			(cos_dir[0] * cos_dir[1] - sin_dir[2] * sin_dir[0]);

		fr[ind * Dim + 1] = -Ao * e_term * args[1] *
			(cos_dir[1] * cos_dir[2] - sin_dir[0] * sin_dir[1]);

		fr[ind * Dim + 2] = -Ao * e_term * args[2] *
			(cos_dir[2] * cos_dir[0] - sin_dir[1] * sin_dir[2]);
	}

}
