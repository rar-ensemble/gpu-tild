#include "globals.h"
#include "group.h"


// Allocate *index and *d_index for
// max possible sites since that is biggest group 
// could be.  
void Group::Initialize() {

    // Allocate memory for group indices on host
    this->nsites = ns;
    this->index = (int*)calloc(ns, sizeof(int));

    // Allocate memory for group indices on device
    cudaMalloc(&this->d_index, ns * sizeof(int));
    device_mem_use += ns * sizeof(int);

    // Initialize GPU grid, block size to ns sizes
    this->grid = ns_Grid;
    this->block = ns_Block;
}

Group::Group() {
}

Group::~Group() {
}


void make_group_all() {
    for (int i = 0; i < ns; i++) {
        Groups.at(0).index[i] = i;
    }
    Groups.at(0).nsites = ns;
    Groups.at(0).block = threads;
    Groups.at(0).grid = (int)ceil((float)(ns) / threads);

    cudaMemcpy(Groups.at(0).d_index, Groups.at(0).index, Groups.at(0).nsites * sizeof(int), cudaMemcpyHostToDevice);

    Groups.at(0).name = "all";
    Groups.at(0).command_line = "NONE";
}

void make_group_type(int grp_index, int type_id) {
    int ntype = 0;
    for (int i = 0; i < ns; i++) {
        if (tp[i] == type_id) {
            Groups[grp_index].index[ntype] = i;
            ntype++;
        }
    }

    Groups[grp_index].nsites = ntype;
    Groups[grp_index].block = threads;
    Groups[grp_index].grid = (int)ceil((float)(ntype) / threads);

    int gid = grp_index;
    cudaMemcpy(Groups[gid].d_index, Groups[gid].index, Groups[gid].nsites * sizeof(int), cudaMemcpyHostToDevice);
    
}


