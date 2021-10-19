#include "globals.h"
#include "pair_style_charges.h"

__global__ void d_prepareChargeDensity(float*, cufftComplex*, int);
__global__ void d_prepareElectrostaticPotential(cufftComplex*, cufftComplex*,
    float, float, const int, const int, const float*, const int*);
__global__ void d_prepareElectricField(cufftComplex*, cufftComplex*,
    float, const int, const int, const float*, const int*, const int);
__global__ void d_divideByDimension(cufftComplex*, const int);
__global__ void d_accumulateGridForceWithCharges(cufftComplex*,
    float*, float*, const int);
__global__ void d_setElectrostaticPotential(cufftComplex*,
    float*, const int);
__global__ void d_resetComplexes(cufftComplex*, cufftComplex*, const int);
__global__ void d_setElectricField(cufftComplex*, float*, int, const int);

void Charges::Initialize() {

}

void Charges::CalcCharges() {
    //zero d_cpx1 and d_cpx2
    d_resetComplexes<<<M_Grid, M_Block>>>(d_cpx1, d_cpx2, M);

    //fft charge density
    d_prepareChargeDensity<<<M_Grid, M_Block>>>(d_charge_density_field, d_cpx1, M);

    cufftExecC2C(fftplan, d_cpx1, d_cpx2, CUFFT_FORWARD);//now fourier transformed density data in cpx2

    cudaReturn = cudaGetLastError();
    if (cudaReturn != cudaSuccess) {
        char cherror[90];
        sprintf(cherror, "Cuda failed with error \"%s\" while cufftExec1 ran\n", cudaGetErrorString(cudaReturn));
        die(cherror);
    }

    d_divideByDimension<<<M_Grid, M_Block>>>(d_cpx2, M);//normalizes the charge density field

    //electric potential in cpx1
    d_prepareElectrostaticPotential<<<M_Grid, M_Block>>>(d_cpx2, d_cpx1, charge_bjerrum_length, charge_smearing_length,
        M, Dim, d_L, d_Nx); 

    cudaReturn = cudaGetLastError();
    if (cudaReturn != cudaSuccess) {
        char cherror[90];
        sprintf(cherror, "Cuda failed with error \"%s\" while d_prepareElectrostaticPotential ran\n",
            cudaGetErrorString(cudaReturn));
        die(cherror);
    }

    for (int j = 0; j < Dim; j++) {
        d_prepareElectricField<<<M_Grid, M_Block>>>(d_cpx2, d_cpx1, charge_smearing_length, M, Dim, d_L, d_Nx, j);//new data for electric field in cpx2

        cudaReturn = cudaGetLastError();
        if (cudaReturn != cudaSuccess) {
            char cherror[90];
            sprintf(cherror, "Cuda failed with error \"%s\" while d_prepareElectricField ran\n",
                cudaGetErrorString(cudaReturn));
            die(cherror);
        }

        cufftExecC2C(fftplan, d_cpx2, d_cpx2, CUFFT_INVERSE); //d_cpx2 now holds the electric field, in place transform

        cudaReturn = cudaGetLastError();
        if (cudaReturn != cudaSuccess) {
            char cherror[90];
            sprintf(cherror, "Cuda failed with error \"%s\" while cufftExec2 ran\n", cudaGetErrorString(cudaReturn));
            die(cherror);
        }

        if (j == 0)
            d_accumulateGridForceWithCharges<<<M_Grid, M_Block>>>(d_cpx2,
                d_charge_density_field, d_all_fx_charges, M);
        if (j == 1)
            d_accumulateGridForceWithCharges<<<M_Grid, M_Block>>>(d_cpx2,
                d_charge_density_field, d_all_fy_charges, M);
        if (j == 2)
            d_accumulateGridForceWithCharges<<<M_Grid, M_Block>>>(d_cpx2,
                d_charge_density_field, d_all_fz_charges, M);

        cudaReturn = cudaGetLastError();
        if (cudaReturn != cudaSuccess) {
            char cherror[90];
            sprintf(cherror, "Cuda failed with error \"%s\" while d_accumulateGridForce in j=%d ran\n",
                cudaGetErrorString(cudaReturn), j);
            die(cherror);
        }

        d_setElectricField<<<M_Grid, M_Block>>>(d_cpx2, d_electric_field, j, M);
    }

    //prepares d_electrostatic_potential to be copied onto host
    cufftExecC2C(fftplan, d_cpx1, d_cpx1, CUFFT_INVERSE);

    cudaReturn = cudaGetLastError();
    if (cudaReturn != cudaSuccess) {
        char cherror[90];
        sprintf(cherror, "Cuda failed with error \"%s\" while cufftExec3 ran\n", cudaGetErrorString(cudaReturn));
        die(cherror);
    }

    d_setElectrostaticPotential<<<M_Grid, M_Block>>>(d_cpx1, d_electrostatic_potential, M);

    cudaReturn = cudaGetLastError();
    if (cudaReturn != cudaSuccess) {
        char cherror[90];
        sprintf(cherror, "Cuda failed with error \"%s\" while d_setElectrostaticPotential ran\n",
            cudaGetErrorString(cudaReturn));
        die(cherror);
    }
}