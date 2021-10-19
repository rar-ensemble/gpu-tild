#ifndef PAIR_STYLE
#define PAIR_STYLE
#include <complex>
#include <cufft.h>
#include <cufftXt.h>
#include "field_component.h"
using namespace std;

class PairStyle {
private:
    int size;
public:
    int type1, type2;
    float* u, ** f, ** vir, energy, * total_vir;
    float* rho1, * rho2, ** force1, ** force2;
    float* d_rho1, * d_rho2;
    complex<float>* u_k, ** f_k, ** vir_k;
    
    float* d_u, * d_f, * d_vir;
    cufftComplex* d_u_k, * d_f_k, * d_vir_k;

    double sigma_squared, Rp, U, initial_prefactor, final_prefactor;
    bool ramp = false;
    bool allocated = false;

    PairStyle();
    ~PairStyle();

    void Initialize_PairStyle(int, int, int);
    float CalcEnergy();
    void CalcVirial();
    void InitializeVirial();
    virtual void CalcForces();

	virtual void Initialize() = 0;
    virtual void Update();
};

#endif
