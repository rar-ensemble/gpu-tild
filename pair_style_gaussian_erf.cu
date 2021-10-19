#include "globals.h"
#include "pair_style_gaussian_erf.h"
#include "device_utils.cuh"
#include <iostream>
#include <fstream>
using namespace std;

__global__ void init_device_gaussian_erf(cufftComplex*, cufftComplex*, 
    float, float, float, const float*,
    const int, const int*, const int);
__global__ void d_multiply_cufftCpx_scalar(cufftComplex*, float, int);
__global__ void d_complex2real(cufftComplex*, float*, int);
__global__ void d_extractForceComp(cufftComplex*, cufftComplex*,
    const int, const int, const int);
__global__ void d_insertForceCompC2R(float*, cufftComplex*, const int,
    const int, const int);

void GaussianErf::Initialize() {
    this->Initialize_Gaussian_Erf(initial_prefactor, sigma_squared, Rp, M, type1, type2);
}

void GaussianErf::Initialize_Gaussian_Erf(float Ao, float sigma2, 
    float Rp, int alloc_size, int typ_A, int typ_B) {
    Initialize_PairStyle(alloc_size, typ_A, typ_B);

    printf("Setting up Gaussian-Erf pair style..."); fflush(stdout);
 
    init_device_gaussian_erf<<<M_Grid, M_Block>>>(this->d_u_k, this->d_f_k,
        Ao, sigma2, Rp, d_L, M, d_Nx, Dim);

    cufftExecC2C(fftplan, this->d_u_k, d_cpx1, CUFFT_INVERSE);
    d_complex2real<<<M_Grid, M_Block>>>(d_cpx1, this->d_u, M);

    for (int j = 0; j < Dim; j++) {
        d_extractForceComp<<<M_Grid, M_Block>>>(d_cpx1, this->d_f_k, j, Dim, M);
        cufftExecC2C(fftplan, d_cpx1, d_cpx1, CUFFT_INVERSE);
        d_insertForceCompC2R<<<M_Grid, M_Block>>>(this->d_f, d_cpx1, j, Dim, M);
    }

    float k2, kv[3], k;

    // Define the potential and the force in k-space
    for (int i = 0; i < alloc_size; i++) {
        k2 = get_k(i, kv, Dim);
		k = sqrt(k2);

		if (k2 == 0) {
			this->u_k[i] = Ao * // prefactor
				// exp(-k2 * sigma2 * 0.5f) * //Gaussian contribution  = 1
				PI4 * Rp * Rp * Rp / 3;   // step function contribution
		}
		else
		{
			this->u_k[i] = Ao * // prefactor
				exp(-k2 * sigma2 * 0.5f) * //Gaussian contribution of both
				PI4 * (sin(Rp * k) - Rp * k * cos(Rp * k)) / (k2 * k); // step function contribution
		}

        for (int j = 0; j < Dim; j++) {
            this->f_k[j][i] = -I * kv[j] * this->u_k[i];
        }
    }
    
    
    InitializeVirial();
    
    printf("done!\n"); fflush(stdout);
}




GaussianErf::GaussianErf() {

}

GaussianErf::~GaussianErf() {

}


__global__ void init_device_gaussian_erf(cufftComplex* uk, cufftComplex* fk,
    float Ao, float sigma2, float Rp,
    const float* dL, const int M, const int* Nx, const int Dim) {


	const int ind = blockIdx.x * blockDim.x + threadIdx.x;
	if (ind >= M)
		return;

	float k2, kv[3], k;

	k2 = d_get_k(ind, kv, dL, Nx, Dim);
	k = sqrt(k2);

	if (k2 == 0) {
		uk[ind].x = Ao *				// prefactor
			//exp(-k2 * sigma2 * 0.5f) * //Gaussian contribution = 1
			PI4 * Rp * Rp * Rp / 3;   // step function contribution of erfc
	}
	else
	{

		uk[ind].x = Ao *				//prefactor
			exp(-k2 * sigma2 * 0.5f) * //Gaussian contribution of both
			PI4 * (sin(Rp * k) - Rp * k * cos(Rp * k)) / (k2 * k);
		// step function for erfc only
	}
	uk[ind].y = 0.f;
	for (int j = 0; j < Dim; j++) {
		fk[ind * Dim + j].x = 0.f;
		fk[ind * Dim + j].y = -kv[j] * uk[ind].x;
	}

}
