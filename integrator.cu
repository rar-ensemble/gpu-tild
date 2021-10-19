#include "globals.h"
#include "integrator.h"
__global__ void d_EM_integrator(float*, float*, float, float,
    float*, float*, int*, int, int, curandState*);
__global__ void d_GJF_integrator(float*, float*, float*, float*,
    float*, float*, int*, float, float, float*, float*, 
    int*, int, int, curandState*);


// Calls first step of VV integrator
void Integrator::Integrate_1() {    
    
    if (this->name == "VV") {
        int grp_id = this->group_index;
        int GRID = Groups[grp_id].grid;
        int BLOCK = Groups[grp_id].block;
    }
    
}


// Calls the EM or GJF integrators
// Or second step of VV integrator
void Integrator::Integrate_2() {
    int grp_id = this->group_index;
    int GRID = Groups[grp_id].grid;
    int BLOCK = Groups[grp_id].block;

    if ( this->name == "EM" ) 
        d_EM_integrator<<<GRID, BLOCK>>>(d_x, d_f,
            delt, noise_mag, d_L, d_Lh, Groups[grp_id].d_index,
            Groups[grp_id].nsites, Dim, d_states);

    else if ( this->name == "GJF")
        d_GJF_integrator<<<GRID, BLOCK>>>(d_x, d_xo, d_f,
            d_prev_noise, d_mass, d_Diff, d_typ,
            delt, noise_mag, d_L, d_Lh, Groups.at(grp_id).d_index,
            Groups.at(grp_id).nsites, Dim, d_states);
}


// int_name: name of the integrator to be used (EM, GJF, eventually DPD)
// grp_name: group to integrate using int_name
void Integrator::Initialize(string int_name, string grp_name) {

    if (int_name != "EM" && int_name != "GJF") {
        char death[120];
        sprintf(death, "%s is not a valid integrator!", int_name.c_str());
        die(death);
    }
    
    this->group_index = -1;
    for (int i = 0; i < n_groups; i++) {
        if (Groups.at(i).name == grp_name)
            this->group_index = i;
    }

    if (this->group_index == -1) {
        char death[120];
        sprintf(death, "Group %s not found to integrate!", grp_name.c_str());
        die(death);
    }

    this->group_name = grp_name;
    this->name = int_name;

    if (int_name == "GJF")
        using_GJF = 1;

}

Integrator::Integrator() {

}

Integrator::~Integrator() {

}