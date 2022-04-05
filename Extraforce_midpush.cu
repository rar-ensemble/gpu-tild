#include <curand_kernel.h>
#include <curand.h>


// Adds Langevin stochastic and friction forces to 
// allow constant T simulations with VV algorithm
__global__ void d_ExtraForce_midpush(
    float* f,               // [ns*Dim], particle forces
    const float* x,         // [ns*Dim], particle positions
    const float fmag,       // Force magnitude
    const int* site_list,   // List of sites in the group
    const float* Lh,        // [Dim] Half box length
    const int ns,           // Number of sites in the list
    const int D) {            // Dimensionality of the simulation
    

    int list_ind = blockIdx.x * blockDim.x + threadIdx.x;
    if (list_ind >= ns)
        return;

    int ind = site_list[list_ind];

    // Operates on (Dim-1) direction
    int pi = ind * D + (D-1);

    if ( x[pi] > Lh[D-1] ) 
        f[pi] -= fmag;
    else
        f[pi] += fmag;

}
