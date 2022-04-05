#include "globals.h"
#include "pair_style.h"

__global__ void d_prepareDensity(int, float*, cufftComplex*, int);
__global__ void d_prepareChargeDensity(float*, cufftComplex*, int);
__global__ void d_prepareForceKSpace(cufftComplex*, cufftComplex*,
    cufftComplex*, const int, const int, const int);
__global__ void d_accumulateGridForce(cufftComplex*, float*, float*,
    const int, const int);
__global__ void d_multiplyComplex(cufftComplex*, cufftComplex*, 
    cufftComplex*, int);
__global__ void d_prepareIntegrand(cufftComplex*, int, float*, float*, int);
__global__ void d_extractForceComp(cufftComplex*, cufftComplex*, const int, const int, const int);
__global__ void d_initVirial(float*, const float*, const float*, const float*, 
    const float*, const int, const int, const int*, const int);
__global__ void d_extractVirialCompR2C(cufftComplex*, const float*, const int,
    const int, const int);
__global__ void d_insertVirialCompC2C(cufftComplex*, const cufftComplex*, const int,
    const int, const int);
__global__ void d_prepareVirial(const cufftComplex*, const cufftComplex*,
    cufftComplex*, const int, const int, const int, const int);
__global__ void d_divideByDimension(cufftComplex*, const int);
__global__ void d_multiply_cufftCpx_scalar(cufftComplex*, float, int);

float reduce_device_float(float*, const int, const int);
//bool checknan_cpx(cufftComplex*, int);
void write_lammps_traj(void);
void cuda_collect_x(void);
void cuda_collect_rho(void);
void write_kspace_cudaComplex(const char*, cufftComplex*);


// Constructor, allocates memory for:
// potential, forces, virial contribution
// in both r- and k-space
void PairStyle::Initialize_PairStyle(int alloc_size, int typ_A, int typ_B) {
    allocated = true;
    int n_off_diag = Dim + (Dim * Dim - Dim) / 2;
        
    size = alloc_size;

    this->u = (float*)calloc(alloc_size, sizeof(float));
    this->f = (float**)calloc(Dim, sizeof(float*));
    this->vir = (float**)calloc(n_off_diag, sizeof(float*));

    this->u_k = (complex<float>*) calloc(alloc_size, sizeof(complex<float>));
    this->f_k = (complex<float>**) calloc(Dim, sizeof(complex<float>*));
    this->vir_k = (complex<float>**) calloc(n_off_diag, sizeof(complex<float>*));

    for (int i = 0; i < Dim; i++) {
        this->f[i] = (float*)calloc(alloc_size, sizeof(float));
        this->f_k[i] = (complex<float>*) calloc(alloc_size, sizeof(complex<float>));
    }

    total_vir = (float*)calloc(n_off_diag, sizeof(float));

    // Vir stores the diagonal plus the off-diagonal terms
    // The term in parenthesis (Dim*Dim-Dim) will always be even
    for (int i = 0; i < n_off_diag; i++) {
        this->vir[i] = (float*)calloc(alloc_size, sizeof(float));
        this->vir_k[i] = (complex<float>*) calloc(alloc_size, sizeof(complex<float>));
    }

    cudaMalloc(&this->d_u, M * sizeof(float));
    cudaMalloc(&this->d_f, Dim * M * sizeof(float));
    cudaMalloc(&this->d_vir, n_P_comps * M * sizeof(float));
    device_mem_use += M * (Dim + 1 + n_P_comps) * sizeof(float);

    cudaMalloc(&this->d_u_k, M * sizeof(cufftComplex));
    cudaMalloc(&this->d_f_k, Dim * M * sizeof(cufftComplex));
    cudaMalloc(&this->d_vir_k, n_P_comps * M * sizeof(cufftComplex));
    device_mem_use += M * (Dim + 1 + n_P_comps) * sizeof(cufftComplex);

    type1 = typ_A;
    type2 = typ_B;

}


// Calculates forces on rho1, rho2 for this pairstyle
void PairStyle::CalcForces() {

    /////////////////////////
    // rho2 acting on rho1 //
    /////////////////////////

    // fft rho2
    d_prepareDensity<<<M_Grid, M_Block>>>(type2, d_all_rho, d_cpx1, M);
  	cudaReturn = cudaGetLastError();
  	if (cudaReturn != cudaSuccess) {
  		char cherror[90];
  		sprintf(cherror, "Cuda failed with error \"%s\" while d_prepareDensity ran\n", 
          cudaGetErrorString(cudaReturn));
  		die(cherror);
  	}

    cufftExecC2C(fftplan, d_cpx1, d_cpx2, CUFFT_FORWARD);
  	cudaReturn = cudaGetLastError();
  	if (cudaReturn != cudaSuccess) {
  		char cherror[90];
  		sprintf(cherror, "Cuda failed with error \"%s\" while cufftExec1 ran\n", cudaGetErrorString(cudaReturn));
  		die(cherror);
  	}

    
    for (int j = 0; j < Dim; j++) {
        // d_cpx1 = d_cpx2 * d_f_k
        d_prepareForceKSpace<<<M_Grid, M_Block>>>(this->d_f_k, 
            d_cpx2, d_cpx1, j, Dim, M);
     	 	cudaReturn = cudaGetLastError();
     	 	if (cudaReturn != cudaSuccess) {
     	 		char cherror[90];
     	 		sprintf(cherror, "Cuda failed with error \"%s\" while d_prepareForceKSpace in j=%d ran\n", cudaGetErrorString(cudaReturn), j);
     	 		die(cherror);
     	 	}
        
        // back to real space, in-place transform
        cufftExecC2C(fftplan, d_cpx1, d_cpx1, CUFFT_INVERSE);
     	 	cudaReturn = cudaGetLastError();
     	 	if (cudaReturn != cudaSuccess) {
     	 		char cherror[90];
     	 		sprintf(cherror, "Cuda failed with error \"%s\" while cufft2 in j=%d ran\n", cudaGetErrorString(cudaReturn), j);
     	 		die(cherror);
     	 	}


        // Accumulate the forces on type 1
        if (j == 0)
            d_accumulateGridForce<<<M_Grid, M_Block>>>(d_cpx1, 
                d_all_rho, d_all_fx, type1, M);
        if (j == 1)
            d_accumulateGridForce<<<M_Grid, M_Block>>>(d_cpx1,
                d_all_rho, d_all_fy, type1, M);
        if (j == 2)
            d_accumulateGridForce<<<M_Grid, M_Block>>>(d_cpx1,
                d_all_rho, d_all_fz, type1, M);
     	 	
        cudaReturn = cudaGetLastError();
     	 	if (cudaReturn != cudaSuccess) {
     	 		char cherror[90];
     	 		sprintf(cherror, "Cuda failed with error \"%s\" while d_accumulateGridForce in j=%d ran\n", 
              cudaGetErrorString(cudaReturn), j);
     	 		die(cherror);
     	 	}
        
    }


    // fft rho1
    d_prepareDensity<<<M_Grid, M_Block>>> (type1, d_all_rho, d_cpx1, M);
    cufftExecC2C(fftplan, d_cpx1, d_cpx2, CUFFT_FORWARD);

  	cudaReturn = cudaGetLastError();
  	if (cudaReturn != cudaSuccess) {
  		char cherror[90];
  		sprintf(cherror, "Cuda failed with error \"%s\" while cufftExec ran\n", cudaGetErrorString(cudaReturn));
  		die(cherror);
  	}


    
    for (int j = 0; j < Dim; j++) {
        // d_cpx1 = d_cpx2 * d_f_k
        d_prepareForceKSpace<<<M_Grid, M_Block>>>(this->d_f_k,
            d_cpx2, d_cpx1, j, Dim, M);

        // back to real space, in-place transform
        cufftExecC2C(fftplan, d_cpx1, d_cpx1, CUFFT_INVERSE);


        // Accumulate the forces on type 2
        if (j == 0)
            d_accumulateGridForce<<<M_Grid, M_Block>>>(d_cpx1,
                d_all_rho, d_all_fx, type2, M);
        if (j == 1)
            d_accumulateGridForce<<<M_Grid, M_Block>>>(d_cpx1,
                d_all_rho, d_all_fy, type2, M);
        if (j == 2)
            d_accumulateGridForce<<<M_Grid, M_Block>>>(d_cpx1,
                d_all_rho, d_all_fz, type2, M);
        
    }// j=0:Dim
    
    /*cudaMemcpy(tmp, d_all_fx, M * sizeof(float), cudaMemcpyDeviceToHost);
    write_grid_data("fx0.dat", tmp);
    cudaMemcpy(tmp, d_all_fy, M * sizeof(float), cudaMemcpyDeviceToHost);
    write_grid_data("fy0.dat", tmp);
    cuda_collect_rho();
    for (int i = 0; i < ntypes; i++) {
        char nm[30];
        sprintf(nm, "rho%d.dat", i);
        write_grid_data(nm, Components[i].rho);
    }
    if ( type1 == 0 && type2 == 1 ) 
        exit(1);*/
}



// Calculates the energy involved in this potential as
// energy = \int dr rho1(r) \int dr' u(r-r') rho2(r')
// The convolution theorem is used to efficiently evaluate
// the integral over r'
float PairStyle::CalcEnergy() {

    // fft rho2
    d_prepareDensity<<<M_Grid, M_Block>>>(type2, d_all_rho, d_cpx1, M);
    cufftExecC2C(fftplan, d_cpx1, d_cpx2, CUFFT_FORWARD);

    // Multiply by d_u_k
    d_multiplyComplex<<<M_Grid, M_Block>>>(this->d_u_k, d_cpx2,
        d_cpx1, M);

    // Back to real space
    cufftExecC2C(fftplan, d_cpx1, d_cpx1, CUFFT_INVERSE);

    d_prepareIntegrand<<<M_Grid, M_Block>>>(d_cpx1, type1, d_all_rho,
        d_tmp, M);

    // Copy the integrand back to host for integration
    // This should be replaced with on-device integration at some point
    cudaMemcpy(tmp, d_tmp, M * sizeof(float), cudaMemcpyDeviceToHost);
    this->energy = integ_trapPBC(tmp);

    //float temp = reduce_device_float(d_tmp, threads, M);
    //cout << "host: " << this->energy << " device: " << temp << endl;
    //die("here!");



    return this->energy;
    
}

void PairStyle::CalcVirial()
{
    // fft rho2
    d_prepareDensity<<<M_Grid, M_Block>>>(this->type2, d_all_rho, d_cpx1, M);
    cufftExecC2C(fftplan, d_cpx1, d_cpx1, CUFFT_FORWARD);
    d_divideByDimension<<<M_Grid, M_Block>>>(d_cpx1, M);


    for (int j = 0; j < n_P_comps; j++) {
        d_prepareVirial<<<M_Grid, M_Block>>>(this->d_vir_k, d_cpx1, d_cpx2,
            j, Dim, n_P_comps, M);

        // back to real space, in-place transform
        cufftExecC2C(fftplan, d_cpx2, d_cpx2, CUFFT_INVERSE);

        d_prepareIntegrand<<<M_Grid, M_Block>>>(d_cpx2, this->type1, d_all_rho,
            d_tmp, M);

        // Copy the integrand back to host for integration
        // This should be replaced with on-device integration at some point
        cudaMemcpy(tmp, d_tmp, M * sizeof(float), cudaMemcpyDeviceToHost);

        this->total_vir[j] = integ_trapPBC(tmp) / V / float(Dim);
    }

}

void PairStyle::InitializeVirial()
{
    
    d_initVirial<<<M_Grid, M_Block>>>(this->d_vir, this->d_f,
        d_L, d_Lh, d_dx, Dim, n_P_comps, d_Nx, M);


    for (int j = 0; j < n_P_comps; j++) {
        d_extractVirialCompR2C<<<M_Grid, M_Block>>>(d_cpx1, this->d_vir, j, n_P_comps, M);
        cufftExecC2C(fftplan, d_cpx1, d_cpx1, CUFFT_FORWARD);
        d_divideByDimension<<<M_Grid, M_Block>>>(d_cpx1, M);
        d_insertVirialCompC2C<<<M_Grid, M_Block>>>(this->d_vir_k, d_cpx1, j, n_P_comps, M);
    }

}



PairStyle::~PairStyle() {
    //printf("PS here for some reason!\n"); fflush(stdout);
    if (allocated){
        int n_off_diag = Dim + (Dim * Dim - Dim) / 2;
        free(this->u);
        free(this->u_k);
    
        for (int i = 0; i < Dim; i++) {
            free(this->f[i]);
            free(this->f_k[i]);
        }

        for (int i = 0; i < n_off_diag; i++) {
            free(this->vir[i]);
            free(this->vir_k[i]);
	}
        free(this->f);
        free(this->f_k);
        free(this->vir);
        free(this->vir_k);
    }

}

PairStyle::PairStyle() {

}

void PairStyle::Update() {
    if (!ramp) return;

    // The factor of 1000x is to attempt to cancel some float
    // precision errors if chi is adjusted slowly over the course
    // of a long simulation. It cancels since the ratio is passed
    // as the argument
    float delta_A = (1000.0f * (final_prefactor - initial_prefactor)) / float(max_steps);
    float cur_A = 1000.0f * initial_prefactor + delta_A * float(step);
    float new_A = cur_A + delta_A;

    cufftComplex temp_cpx;
    cudaMemcpy(&temp_cpx, this->d_u_k, sizeof(cufftComplex), cudaMemcpyDeviceToHost);
    cur_A = temp_cpx.x * 1000.0f;



    d_multiply_cufftCpx_scalar << <M_Grid, M_Block >> > (this->d_u_k, new_A / cur_A, M);
    d_multiply_cufftCpx_scalar << <M_Grid * Dim, M_Block >> > (this->d_f_k,
        new_A / cur_A, M * Dim);
}
