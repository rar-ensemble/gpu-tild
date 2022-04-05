#include <curand_kernel.h>
#include <curand.h>


// Adds Langevin stochastic and friction forces to 
// allow constant T simulations with VV algorithm
__global__ void d_ExtraForce_Langevin(
    float* f,               // [ns*Dim], particle forces
    const float* v,         // [ns*Dim], particle velocities
    const float noise_mag,  // magnitude of the noise, should be sqrt(2.*gamma)
    const float gamma,      // Friction force
    const int* site_list,   // List of sites in the group
    const int ns,           // Number of sites in the list
    const int D,            // Dimensionality of the simulation
    curandState* d_states) {// Status of CUDA rng

    int list_ind = blockIdx.x * blockDim.x + threadIdx.x;
    if (list_ind >= ns)
        return;

    int ind = site_list[list_ind];

    curandState l_state;

    l_state = d_states[ind];

    for (int j = 0; j < D; j++)
        f[ind * D + j] += -gamma * v[ind * D + j] + noise_mag * curand_normal(&l_state);

    d_states[ind] = l_state;

}
