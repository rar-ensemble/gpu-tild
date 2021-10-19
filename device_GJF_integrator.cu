#include <curand_kernel.h>
#include <curand.h>

__global__ void d_GJF_integrator(float* x, float* xo, 
	float* f, float* old_noise,
	float *mass, float *diff, int *typ, float delt, float noise_mag, 
	float* L, float* Lh, 
	int *site_list,
	int ns, 
	int Dim, 
	curandState* d_states ) {

	int list_ind = blockIdx.x * blockDim.x + threadIdx.x;
	if (list_ind >= ns)
		return;

	int ind = site_list[list_ind];

	curandState l_state;

	l_state = d_states[ind];

	int itype = typ[ind];

	float b = mass[itype] / (mass[itype] + delt / (2.0f * diff[itype]));
	
	float a = (1.0f - delt / (2.0f * diff[itype] * mass[itype])) /
		(1.0f + delt / (2.0f * diff[itype] * mass[itype]));
	
	float delt2 = delt * delt;

	for (int j = 0; j < Dim; j++) {
		int aind = ind * Dim + j;

		// Not sure this should be divided by diff[itype]
		float new_noise = noise_mag / diff[itype] * curand_normal(&l_state);
		
		float xtmp = x[aind];

		if (xo[aind] - x[aind] > L[j] / 2.0f)
			xo[aind] -= L[j];
		else if (xo[aind] - x[aind] < -L[j] / 2.0f)
			xo[aind] += L[j];

		
		x[aind] = 2.0f * b * x[aind] - a * xo[aind]
			+ b * delt2 / mass[itype] * f[aind]
			+ b * delt / ( 2.0f * mass[itype] ) * (new_noise + old_noise[aind]);


		xo[aind] = xtmp;
		old_noise[aind] = new_noise;


		if (x[ind * Dim + j] > L[j])
			x[ind * Dim + j] -= L[j];
		else if (x[ind * Dim + j] < 0.0f)
			x[ind * Dim + j] += L[j];
		
	}

	d_states[ind] = l_state;
}